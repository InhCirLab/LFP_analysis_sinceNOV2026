%% ========================================================================
%  LFP ANALYSIS PIPELINE — Post-HippoCookbook preprocessing (REORDENADO)
%  ========================================================================
%
%  PROBLEMA ORIGINAL
%    intervals_nonTheta_nonRipples (usado por Parte I nonrun y por
%    ps1to120 -> allrec / loc / nonlocnonrip) se calcula en detectThetaEpochs
%    LEYENDO el <rec>.ripples.events.mat que exista en ese momento. En el
%    orden antiguo (I -> III -> II) los ripples nuevos de la Parte III se
%    generaban DESPUES de guardar .thetaEpochs6to12.states.mat, asi que
%    Parte I (nonrun) y Parte II usaban los ripples VIEJOS.
%
%  NUEVO ORDEN DE EJECUCION
%    1. P1  PASADA 1  : theta + run (genera thetaRun.newintervals que
%                       necesita la CNN para excluir locomocion). Sin Excel.
%    2. P3  RIPPLES   : deteccion CNN -> guarda <rec>.ripples.events.mat NUEVO
%    3. P1R REFRESH   : vuelve a correr detectThetaEpochs (ya lee los ripples
%                       nuevos), recalcula nonrun, reescribe
%                       .thetaEpochs6to12.states.mat y guarda Excel/figuras.
%                       Imprime por sesion si intervals_nonTheta_nonRipples
%                       cambio respecto a la version vieja.
%    4. P2  POWER SPECTRUM : ps1to120_Alberto lee el .states.mat ya refrescado.
%
%  ps1to120_Alberto.m NO necesita cambios.
%
%  Las funciones auxiliares (extractEpochPower, storeRunPower, buildRow,
%  saveFigSession, saveGroupFig, savePSFig, robust_zscore_matrix, ...) siguen
%  siendo las de siempre, en el path.
% =========================================================================

clearvars; close all; clc;

%% ── POOL PARALELO ────────────────────────────────────────────────────────
%  Se usa en Parte I (un worker por registro) y en ps1to120_Alberto (un worker
%  por sesion). La CNN de la Parte III se deja en serie (Python/GPU + RAM).
N_WORKERS = 8;   % ajusta segun RAM: cada worker carga LFP + wavelets de UNA sesion
if license('test','Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
    parpool('local', min(N_WORKERS, maxNumCompThreads));
end

%% ========================================================================
%  PARTE I — THETA & GAMMA DURANTE LOCOMOTION
%% ========================================================================

%% ── P1.0  PARAMETROS ─────────────────────────────────────────────────────

P.THETA_BAND  = [4  12];
P.LGAMMA_BAND = [20 60];
P.HGAMMA_BAND = [60 100];
P.N_CYCLES    = 8;
P.MIN_RUN_DUR = 1;   % segundos — epocas de carrera mas cortas se descartan

P.DETECT_OPTS = {'saveMat', false, 'force', true, ...
                 'bandpass', [6 12], 'delta_bandpass', [1 3.5], ...
                 'powerThreshold', 1};

P.powerVarNames = ["events (n)", "totaldur (s)", "meandur (s)", ...
    "medianTheta (uv)",  "medianLgamma (uv)",  "medianHgamma (uv)", ...
    "medianTheta (db)",  "medianLgamma (db)",  "medianHgamma (db)", ...
    "zTheta (au)",       "zLgamma (au)",       "zHgamma (au)",      ...
    "zmodTheta (au)",    "zmodLgamma (au)",    "zmodHgamma (au)"];

%% ── P1.1  CONFIGURACION: GRUPOS Y CARPETAS ───────────────────────────────

nGroupsStr = inputdlg('Cuantos grupos experimentales tienes?', 'Grupos', 1, {'1'});
if isempty(nGroupsStr), error('Abortado.'); end
nGroups = str2double(nGroupsStr{1});
assert(~isnan(nGroups) && nGroups >= 1, 'Numero de grupos no valido.');

groupFolders = cell(nGroups, 1);
groupLabels  = cell(nGroups, 1);

for g = 1:nGroups
    lbl = inputdlg(sprintf('Etiqueta del grupo %d (ej: PV, SST, WT...):', g), ...
                   sprintf('Grupo %d', g), 1, {sprintf('G%d', g)});
    if isempty(lbl), error('Abortado.'); end
    groupLabels{g} = strtrim(lbl{1});

    fprintf('>>> Selecciona la carpeta del grupo "%s":\n', groupLabels{g});
    folder = uigetdir(pwd, sprintf('Carpeta del grupo: %s', groupLabels{g}));
    if folder == 0, error('No se selecciono carpeta. Abortado.'); end
    groupFolders{g} = folder;
end

fprintf('>>> Selecciona la carpeta OUTPUT para los resultados de la Parte I:\n');
outFolder = uigetdir(groupFolders{1}, 'OUTPUT — Parte I');
if outFolder == 0, error('Abortado.'); end

figFolder = fullfile(outFolder, 'figures');
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

%% ── P1.2  PASADA 1: THETA + RUN (solo para generar newintervals) ────────
%  Necesaria porque la CNN (Parte III) excluye thetaRun.newintervals.
%  No guarda Excel ni figuras: eso se hace en P1R con los ripples nuevos.

fprintf('\n##### PASADA 1 (theta + run) #####\n');
runThetaStage(groupFolders, groupLabels, P, figFolder, true);


%% ========================================================================
%  PARTE III — RIPPLE ANALYSIS
%% ========================================================================

%% ── P3.0  PARAMETROS ─────────────────────────────────────────────────────

% Rutas del modelo CNN y ejecutable Python
CNN_MODEL_PATH = 'C:\Users\usuario\cnn-matlab-master\cnn';
CNN_PYTHON_EXE = 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe';

% Parametros de deteccion
RIPPLE_THRESHOLD = [2 5];        % [min max] en SD
RIPPLE_PASSBAND  = [80 240];     % Hz
PEAKNORM_CUTOFF  = 50;           % outliers con PeakNormPow > este valor se descartan

% Analisis robusto de outliers — thresholds
OUTLIER_THRESH_MAD  = 3.5;
OUTLIER_MIN_N_REC   = 20;
OUTLIER_MIN_N_SUB   = 20;
OUTLIER_MIN_N_GROUP = 20;
OUTLIER_VARS        = {'PeakAmp','PeakNormPow','Duration','PeakFq'};
OUTLIER_LOG_VARS    = {'PeakAmp','PeakNormPow','Duration'};

%% ── P3.1  CONFIGURACION: GRUPOS Y CARPETAS ───────────────────────────────

nGroupsStr3 = inputdlg('Cuantos grupos para el analisis de ripples?', 'Ripples — Grupos', 1, {num2str(nGroups)});
if isempty(nGroupsStr3), error('Abortado.'); end
nGroupsR = str2double(nGroupsStr3{1});
assert(~isnan(nGroupsR) && nGroupsR >= 1, 'Numero de grupos no valido.');

ripGroupFolders = cell(nGroupsR, 1);
ripGroupLabels  = cell(nGroupsR, 1);

for g = 1:nGroupsR
    defaultLbl = '';
    if g <= nGroups, defaultLbl = groupLabels{g}; else, defaultLbl = sprintf('G%d', g); end
    lbl3 = inputdlg(sprintf('Etiqueta del grupo %d:', g), ...
                    sprintf('Ripples Grupo %d', g), 1, {defaultLbl});
    if isempty(lbl3), error('Abortado.'); end
    ripGroupLabels{g} = strtrim(lbl3{1});

    fprintf('>>> Selecciona la carpeta del grupo "%s":\n', ripGroupLabels{g});
    folder3 = uigetdir(pwd, sprintf('Carpeta del grupo: %s', ripGroupLabels{g}));
    if folder3 == 0, error('No se selecciono carpeta. Abortado.'); end
    ripGroupFolders{g} = folder3;
end

fprintf('>>> Selecciona la carpeta OUTPUT para los resultados de Ripples:\n');
rippleFolder = uigetdir(ripGroupFolders{1}, 'OUTPUT — Parte III Ripples');
if rippleFolder == 0, error('Abortado.'); end

%% ── P3.2  LISTAR SUBCARPETAS POR GRUPO ───────────────────────────────────

lfpPaths = cell(nGroupsR, 1);

for g = 1:nGroupsR
    folderInfo = dir(ripGroupFolders{g});
    subfolders = folderInfo([folderInfo.isdir]);
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
    lfpPaths{g} = {subfolders.name};
    fprintf('  Grupo "%s": %d sesion(es) encontradas.\n', ripGroupLabels{g}, numel(lfpPaths{g}));
end

%% ── P3.3  DETECCION Y CLASIFICACION DE RIPPLES (POR GRUPO) ──────────────

allRipplesStruct = cell(nGroupsR, 1);
allCSD           = cell(nGroupsR, 1);
allCSDEvents     = cell(nGroupsR, 1);
allCSDSummary    = cell(nGroupsR, 1);
sinkSizes        = struct();

for g = 1:nGroupsR
    grpLabel = ripGroupLabels{g};
    paths    = lfpPaths{g};
    nSess    = numel(paths);

    fprintf('\n==== Procesando grupo "%s": %d sesion(es) ====\n', grpLabel, nSess);

    ripplesCell = cell(nSess, 1);
    csdCell     = cell(nSess, 1);
    csdEvCell   = cell(nSess, 1);
    csdSumCell  = cell(nSess, 1);

    slmSize   = nan(nSess, 1);
    radSize   = nan(nSess, 1);
    otherSize = nan(nSess, 1);

    % Evita "Graphics timeout" con figuras invisibles durante deteccion batch
    set(0, 'DefaultFigureVisible', 'off');

    for ii = 1:nSess
        thisFolderName = fullfile(ripGroupFolders{g}, paths{ii});
        fprintf('  [%d/%d] %s\n', ii, nSess, paths{ii});
        cd(thisFolderName);

        % Cargar canales hippocampales
        layersFile = fullfile(thisFolderName, [paths{ii} '.hippocampalLayers.channelinfo.mat']);
        load(layersFile, 'hippocampalLayers');
        rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
        radCh         = hippocampalLayers.layers{1,1}.radiatum;
        SWChannel     = hippocampalLayers.layers{1,1}.radiatum;
        slmCh         = hippocampalLayers.layers{1,1}.slm;

        % Cargar intervalos theta (generados en la PASADA 1 de la Parte I)
        thetaFile = fullfile(thisFolderName, [paths{ii} '.thetaEpochs6to12.states.mat']);
        load(thetaFile, 'thisthetaEpochs');
        thetaIntervals = thisthetaEpochs.thetaRun.newintervals;

        % Cargar LFP
        lfp = getLFP('all');

        try
            % Deteccion con CNN
            ripples = rippleMasterDetector('rippleChannel', rippleChannel, 'SWChannel', SWChannel, ...
                'detector', 'customcnn', 'exec_env', CNN_PYTHON_EXE, 'model_file', CNN_MODEL_PATH, ...
                'passband', RIPPLE_PASSBAND, 'force', true, 'skipStimulationPeriods', true, ...
                'thresholds', RIPPLE_THRESHOLD, 'eventSpikeThreshold', false, ...
                'excludeIntervals', thetaIntervals, ...
                'rippleStats', true);

            try
                ripples = computeRippleStats_new( ...
                    'ripples',       ripples, ...
                    'rippleChannel', rippleChannel, ...
                    'bandpass',      RIPPLE_PASSBAND, ...
                    'thetaEpochs',   thisthetaEpochs, ...
                    'saveSummary',   false, ...
                    'saveMat',       false, ...
                    'plotOpt',       false);
                fprintf('    computeRippleStats_new OK: %d ripples\n', numel(ripples.peaks));
            catch ME_stats
                warning('computeRippleStats_new fallo en sesion %s:\n  %s\nSesion descartada.', ...
                    paths{ii}, ME_stats.message);
                fprintf('    Stack trace:\n');
                for sti = 1:numel(ME_stats.stack)
                    fprintf('      %s (linea %d)\n', ME_stats.stack(sti).name, ME_stats.stack(sti).line);
                end
                ripplesCell{ii} = []; csdCell{ii} = []; csdEvCell{ii} = []; csdSumCell{ii} = [];
                close all; drawnow;
                continue;
            end

            % Fix valid_events: sincronizar rippleStats.data al numero de peaks
            nRipClean = numel(ripples.peaks);
            nRipStats = size(ripples.rippleStats.data.peakAmplitude, 1);
            if nRipStats ~= nRipClean
                fprintf('    FIX valid_events: rippleStats.data tenia %d filas vs %d peaks. Recortando...\n', nRipStats, nRipClean);
                statFields = fieldnames(ripples.rippleStats.data);
                for sf = 1:numel(statFields)
                    fdata = ripples.rippleStats.data.(statFields{sf});
                    if isnumeric(fdata) && size(fdata,1) == nRipStats
                        ripples.rippleStats.data.(statFields{sf}) = fdata(1:nRipClean, :);
                    end
                end
            end

            % ── GUARDAR <rec>.ripples.events.mat ─────────────────────────
            %    Lo leen detectThetaEpochs (en P1R) y nonlocnonrip.
            %    IMPORTANTE: la variable guardada se llama 'ripples' (convencion
            %    buzcode: <rec>.ripples.events.mat contiene la variable 'ripples').
            %    Antes se guardaba como 'ripplesEvents', que un loader buzcode
            %    puede no reconocer.
            ripplesEvents.timestamps  = ripples.timestamps;
            ripplesEvents.peaks       = ripples.peaks;
            ripplesEvents.amplitude   = ripples.peakNormedPower;
            ripplesEvents.duration    = diff(ripples.timestamps, 1, 2);
            ripplesEvents.detectorinfo.detectorname = 'rippleMasterDetector';
            ripplesEvents.detectorinfo.detectionparms.passband   = RIPPLE_PASSBAND;
            ripplesEvents.detectorinfo.detectionparms.thresholds = RIPPLE_THRESHOLD;
            ripplesEventsFile = fullfile(thisFolderName, [paths{ii} '.ripples.events.mat']);
            S_rip.ripples = ripplesEvents;                 % struct auxiliar -> variable 'ripples'
            save(ripplesEventsFile, '-struct', 'S_rip', '-mat');
            clear S_rip ripplesEvents;
            fprintf('    Guardado: %s.ripples.events.mat\n', paths{ii});

            % CSD por ripple
            [csd, ~, csdEvents, ~] = buzcode_eacheventCSD(lfp, ripples.peaks, ...
                'plotCSD', false, 'plotLFP', false);

            % Resumir CSD por ripple (media temporal centrada)
            nRip = size(csdEvents, 3);
            csdSum.tsRipple    = ripples.timestamps;
            csdSum.durRipple   = diff(csdSum.tsRipple, 1, 2) * 1000;  % ms
            csdSum.rangeRipple = round(csdSum.durRipple / 2);

            srCSD = csd.samplingRate;
            nCSDchan = size(csdEvents, 2);
            csdSum.csdSummary = nan(nRip, nCSDchan);

            for jj = 1:nRip
                halfDur_samp = max(round(csdSum.rangeRipple(jj) / 1000 * srCSD), 1);
                idxWin = csd.timestamps >= -halfDur_samp & csd.timestamps <= halfDur_samp;
                if ~any(idxWin)
                    idxWin = true(size(csd.timestamps));  % fallback: ventana completa
                end
                slice   = csdEvents(idxWin, :, jj);
                rowMean = mean(slice, 1);
                csdSum.csdSummary(jj, :) = rowMean;
            end

            % Categorizar por sink (RAD / SLM / Other)
            [~, sinkIdx] = min(csdSum.csdSummary, [], 2);
            slmPeaks   = (sinkIdx >= radCh + 1);
            radPeaks   = (sinkIdx >= rippleChannel + 1) & (sinkIdx <= radCh);
            otherPeaks = ~radPeaks & ~slmPeaks;

            slmSize(ii)   = sum(slmPeaks);
            radSize(ii)   = sum(radPeaks);
            otherSize(ii) = sum(otherPeaks);

            fprintf('    SLM: %d | RAD: %d | Other: %d\n', slmSize(ii), radSize(ii), otherSize(ii));

            ripples.rippleStats.data.radpeaksidx = radPeaks;
            ripples.rippleStats.data.slmpeaksidx = slmPeaks;
            ripples.rippleStats.data.others      = otherPeaks;

            ripplesCell{ii} = ripples;
            csdCell{ii}     = csd;
            csdEvCell{ii}   = csdEvents;
            csdSumCell{ii}  = csdSum;

        catch ME
            warning('Error en sesion %s: %s', paths{ii}, ME.message);
            ripplesCell{ii} = [];
            csdCell{ii}     = [];
            csdEvCell{ii}   = [];
            csdSumCell{ii}  = [];
        end
        close all; drawnow;
    end

    % Restaurar visibilidad de figuras tras deteccion batch
    set(0, 'DefaultFigureVisible', 'on');

    allRipplesStruct{g} = ripplesCell;
    allCSD{g}           = csdCell;
    allCSDEvents{g}     = csdEvCell;
    allCSDSummary{g}    = csdSumCell;

    sinkSizes.(grpLabel).slm   = slmSize;
    sinkSizes.(grpLabel).rad   = radSize;
    sinkSizes.(grpLabel).other = otherSize;

    fprintf('  Grupo "%s" completado.\n', grpLabel);
end

% Guardar estructuras completas de ripples
allRipplesFile = fullfile(rippleFolder, 'allrippledata.mat');
saveVars = cellfun(@(lbl) [lbl 'ripples'], ripGroupLabels, 'UniformOutput', false);
for g = 1:nGroupsR
    eval(sprintf('%s = allRipplesStruct{%d};', saveVars{g}, g));
end
save(allRipplesFile, saveVars{:}, '-mat');
fprintf('\nEstructuras de ripples guardadas: %s\n', allRipplesFile);

%% ── P3.4  EXTRAER MEDIDAS Y ENSAMBLAR TABLAS ────────────────────────────

allRipplesTables = cell(nGroupsR, 1);
allRipplesRaw    = cell(nGroupsR, 1);
allRipplesFilt   = cell(nGroupsR, 1);

for g = 1:nGroupsR
    grpLabel = ripGroupLabels{g};
    paths    = lfpPaths{g};
    ripplesCell = allRipplesStruct{g};

    measures  = [];
    rawWaves  = [];
    filtWaves = [];

    for ii = 1:numel(ripplesCell)
        if isempty(ripplesCell{ii}), continue; end
        ripples = ripplesCell{ii};

        % Extraer medidas
        peakNormPow = ripples.peakNormedPower;
        peakFq      = ripples.rippleStats.data.peakFrequency;
        peakAmp     = ripples.rippleStats.data.peakAmplitude;
        duration    = ripples.rippleStats.data.duration;
        interRipFq  = ripples.rippleStats.data.interRippleFrequency;
        specEnt     = ripples.rippleStats.data.spectralEntropy(:);
        fastRipIdx  = ripples.rippleStats.data.fastRippleIndex(:);
        multitap    = ripples.rippleStats.data.multiTapperFreq;

        rippleIFstat    = ripples.rippleStats.data.rippleIFstat;
        rippleIF_meanFq = ripples.rippleStats.data.rippleIF.meanFreq;
        rippleIF_sdFq   = ripples.rippleStats.data.rippleIF.stdFreq;
        rippleIF_drift  = ripples.rippleStats.data.rippleIF.drift;

        % Sink type
        sinkType = strings(numel(peakNormPow), 1);
        sinkType(ripples.rippleStats.data.radpeaksidx) = "rad";
        sinkType(ripples.rippleStats.data.slmpeaksidx) = "slm";
        sinkType(ripples.rippleStats.data.others)      = "other";

        thisID  = repmat(string(paths{ii}), numel(peakNormPow), 1);
        expCond = repmat(string(grpLabel),  numel(peakNormPow), 1);

        thisTable = table(thisID, expCond, peakNormPow, peakFq, peakAmp, duration, interRipFq, ...
            specEnt, fastRipIdx, multitap, rippleIFstat, rippleIF_meanFq, rippleIF_sdFq, rippleIF_drift, sinkType, ...
            'VariableNames', ["ID","expcond","PeakNormPow","PeakFq","PeakAmp","Duration","interrippleFq", ...
            "spectralEntropy","FastRippleIdx","multitap","freqstat","meanfq_nonstat","std_nonstat","drift_nonstat","sinktype"]);

        measures = [measures; thisTable];

        rawWaves  = [rawWaves;  array2table(ripples.rippleStats.maps.ripples_raw,      'VariableNames', string(ripples.rippleStats.maps.timestamps))];
        filtWaves = [filtWaves; array2table(ripples.rippleStats.maps.ripples_filtered, 'VariableNames', string(ripples.rippleStats.maps.timestamps))];
    end

    allRipplesTables{g} = measures;
    allRipplesRaw{g}    = rawWaves;
    allRipplesFilt{g}   = filtWaves;
end

% Combinar todos los grupos
allRipples     = vertcat(allRipplesTables{:});
allRipplesRaw  = vertcat(allRipplesRaw{:});
allRipplesFilt = vertcat(allRipplesFilt{:});

% Guardar tabla original
writetable(allRipples, fullfile(rippleFolder, 'all_ripples_original.xlsx'));

% Filtrar outliers por PeakNormPow > PEAKNORM_CUTOFF
keepIdx = allRipples.PeakNormPow < PEAKNORM_CUTOFF;
allRipplesKeep = allRipples(keepIdx, :);

% Guardar por grupo
for g = 1:nGroupsR
    grpIdx = strcmp(allRipples.expcond, ripGroupLabels{g});
    grpKeepIdx = grpIdx & keepIdx;
    writetable(allRipples(grpIdx, :),     fullfile(rippleFolder, sprintf('%s_rippledata_original.xlsx', ripGroupLabels{g})));
    writetable(allRipples(grpKeepIdx, :), fullfile(rippleFolder, sprintf('%s_rippledata_keep.xlsx',     ripGroupLabels{g})));
end

% Guardar tabla combinada filtrada
writetable(allRipplesKeep, fullfile(rippleFolder, 'all_ripples_keep.xlsx'));

% Guardar waveforms
waveformFile = fullfile(rippleFolder, 'ripplewaveform.mat');
save(waveformFile, 'allRipplesRaw', 'allRipplesFilt', '-mat');

fprintf('\nTablas de ripples guardadas en: %s\n', rippleFolder);

%% ── P3.5  AÑADIR COLUMNA SUBJECT ─────────────────────────────────────────

if ~ismember('subject', allRipplesKeep.Properties.VariableNames)
    % Extrae ID de animal del nombre de sesion: "C089_190520_sess1" -> "C089"
    allRipplesKeep.subject = extractBefore(allRipplesKeep.ID, '_');
end

%% ── P3.6  OUTLIER REMOVAL (JERARQUICO ROBUSTO) ──────────────────────────

data = allRipplesKeep;
fprintf('\n=== Removing outliers (hierarchical robust MAD) ===\n');

missingVars = setdiff(OUTLIER_VARS, data.Properties.VariableNames);
if ~isempty(missingVars)
    error('Columnas no encontradas en la tabla: %s\nColumnas disponibles: %s', ...
        strjoin(missingVars, ', '), strjoin(data.Properties.VariableNames, ', '));
end

n = height(data);
flagsMat = false(n, numel(OUTLIER_VARS));

uniqueRecs   = unique(data.ID);
uniqueSubs   = unique(data.subject);
uniqueGroups = unique(data.expcond);

for vi = 1:numel(OUTLIER_VARS)
    var  = OUTLIER_VARS{vi};
    xAll = data.(var);
    flagsVar = false(n, 1);
    isLog = ismember(var, OUTLIER_LOG_VARS);

    for r = 1:numel(uniqueRecs)
        recID  = uniqueRecs(r);
        idxRec = strcmp(data.ID, recID);
        valid  = idxRec & ~isnan(xAll);
        if isLog, valid = valid & (xAll > 0); end
        x = xAll(valid);
        if isempty(x), continue; end

        if numel(x) >= OUTLIER_MIN_N_REC
            xproc = x; validUse = valid;
        else
            subjID  = unique(data.subject(idxRec));
            idxSub  = strcmp(data.subject, subjID);
            validSub = idxSub & ~isnan(xAll);
            if isLog, validSub = validSub & (xAll > 0); end
            xSub = xAll(validSub);
            if numel(xSub) >= OUTLIER_MIN_N_SUB
                xproc = xSub; validUse = validSub;
            else
                grpID   = unique(data.expcond(idxRec));
                idxGrp  = strcmp(data.expcond, grpID);
                validGrp = idxGrp & ~isnan(xAll);
                if isLog, validGrp = validGrp & (xAll > 0); end
                xGrp = xAll(validGrp);
                if numel(xGrp) >= OUTLIER_MIN_N_GROUP
                    xproc = xGrp; validUse = validGrp;
                else
                    validGlob = ~isnan(xAll);
                    if isLog, validGlob = validGlob & (xAll > 0); end
                    xGlob = xAll(validGlob);
                    if isempty(xGlob), continue; end
                    xproc = xGlob; validUse = validGlob;
                end
            end
        end

        if isLog, xproc = log10(xproc(xproc>0)); end
        med = median(xproc);
        madv = median(abs(xproc - med));
        scaledMAD = 1.4826 * madv;
        if scaledMAD == 0, scaledMAD = max(std(xproc), eps); end

        idxValidList = find(validUse);
        xvalid = xAll(validUse);
        if isLog, xvalid = log10(xvalid); end
        z = (xvalid - med) / scaledMAD;
        outIdx = idxValidList(abs(z) > OUTLIER_THRESH_MAD);
        flagsVar(outIdx) = true;
    end

    flagsMat(:, vi) = flagsVar;
    fprintf('  %s: %d outliers\n', var, sum(flagsVar));
end

unionIdx = any(flagsMat, 2);
fprintf('TOTAL outliers (union): %d\n', sum(unionIdx));

data = data(~unionIdx, :);
fprintf('Data after outlier removal: %d rows\n', height(data));

%% ── P3.7  HISTOGRAMAS POR GRUPO ──────────────────────────────────────────

vars = {'PeakFq','PeakAmp','spectralEntropy','Duration'};

if nGroupsR <= 3
    cols = [1.00 0.55 0.75; 0.65 0.85 1.00; 0.05 0.25 0.65];
    cols = cols(1:nGroupsR, :);
else
    cols = lines(nGroupsR);
end

alphaVal = 0.35;
nbins    = 30;

fprintf('\n=== Diagnostico de variables para histogramas ===\n');
for i = 1:numel(vars)
    varname = vars{i};
    if ismember(varname, data.Properties.VariableNames)
        nvalid = sum(~isnan(data.(varname)) & isfinite(data.(varname)));
        fprintf('  %s: %d valores validos\n', varname, nvalid);
    else
        fprintf('  %s: *** COLUMNA NO ENCONTRADA ***\n', varname);
        fprintf('    Columnas disponibles: %s\n', strjoin(data.Properties.VariableNames, ', '));
    end
end
fprintf('=================================================\n\n');

fig = figure('Units','normalized','Position',[0.05 0.15 0.9 0.6],'Color','w');
for i = 1:numel(vars)
    ax = subplot(1, 4, i);
    hold(ax, 'on');
    varname = vars{i};

    if ~ismember(varname, data.Properties.VariableNames)
        title(ax, [strrep(varname,'_',' ') ' (NOT FOUND)']);
        continue;
    end

    combined = data.(varname);
    combined = combined(~isnan(combined) & isfinite(combined));
    if isempty(combined), title(ax, strrep(varname,'_',' ')); continue; end
    edges = linspace(min(combined), max(combined), nbins+1);

    plotted   = gobjects(0);
    legLabels = {};
    for g = 1:nGroupsR
        grp = ripGroupLabels{g};
        idx = strcmp(data.expcond, grp);
        grpdata = data.(varname)(idx);
        grpdata = grpdata(~isnan(grpdata) & isfinite(grpdata));
        if isempty(grpdata), continue; end
        h = histogram(ax, grpdata, edges, 'Normalization','probability', ...
            'FaceAlpha', alphaVal, 'EdgeColor','none', 'FaceColor', cols(g,:));
        plotted(end+1) = h;
        legLabels{end+1} = grp;
    end

    title(ax, strrep(varname,'_',' '), 'Interpreter','none');
    xlabel(ax,'Value');
    if i == 1, ylabel(ax,'Probability'); end
    if ~isempty(plotted), legend(ax, plotted, legLabels, 'Location','best'); end
    set(ax,'FontSize',11,'Box','off');
    hold(ax,'off');
end

saveas(fig, fullfile(rippleFolder, 'all_ripples_histograms.png'));
close(fig);

fprintf('\n[Parte III] Completada. Archivos .ripples.events.mat guardados en cada carpeta de sesion.\n');
fprintf('           Resultados de analisis en: %s\n\n', rippleFolder);


%% ========================================================================
%  PARTE I (REFRESH) — RECALCULO CON RIPPLES NUEVOS
%  Ahora detectThetaEpochs lee los .ripples.events.mat recien generados.
%  Recalcula nonrun, reescribe .thetaEpochs6to12.states.mat (que lee
%  ps1to120) y guarda Excel + figuras de la Parte I.
%% ========================================================================

%% ── P1R.1  REFRESH ───────────────────────────────────────────────────────

fprintf('\n##### REFRESH Parte I (con ripples nuevos) #####\n');
Tout = runThetaStage(groupFolders, groupLabels, P, figFolder, false);

%% ── P1R.2  GUARDAR EXCEL ─────────────────────────────────────────────────

writetable(Tout.allrec,  fullfile(outFolder, 'theta_allrec.xlsx'));
writetable(Tout.run,     fullfile(outFolder, 'theta_run.xlsx'));
writetable(Tout.nonrun,  fullfile(outFolder, 'theta_nonrun.xlsx'));
writetable(Tout.loc,     fullfile(outFolder, 'locevent_table.xlsx'));
fprintf('\n[Parte I] Excel guardados en: %s\n', outFolder);

%% ── P1R.3  FIGURAS GRUPALES ──────────────────────────────────────────────

saveGroupFig(Tout.run, Tout.nonrun, groupLabels, figFolder);
fprintf('[Parte I] Completada.\n\n');


%% ========================================================================
%  PARTE II — POWER SPECTRUM 1-120 Hz
%  (lee .thetaEpochs6to12.states.mat ya refrescado con ripples nuevos)
%% ========================================================================

%% ── P2.0  PARAMETROS DE BANDAS ───────────────────────────────────────────

PS_BANDS = struct( ...
    'delta',       struct('active', true, 'range', [1    3.5]), ...
    'theta',       struct('active', true, 'range', [6    12]),  ...
    'beta',        struct('active', true, 'range', [12   30]),  ...
    'lowgamma',    struct('active', true, 'range', [30   60]),  ...
    'mediumgamma', struct('active', true, 'range', [60   90]),  ...
    'highgamma',   struct('active', true, 'range', [90   120]));

LAYER_NAMES = {'or', 'pyr', 'rad', 'slm'};

PS_CONDITIONS = {
    'allrec',              'AllRec_dB',       true;
    'allrecinuv',          'AllRec_uV',       true;
    'locindb',             'Loc_dB',          true;
    'locinuv',             'Loc_uV',          true;
    'nonlocnonrippleindb', 'NonLocNonRip_dB', true;
    'nonlocnonrippleinuv', 'NonLocNonRip_uV', true;
    'nonlocindb',          'NonLoc_dB',       true;
    'nonlocinuv',          'NonLoc_uV',       true;
};

%% ── P2.1  CONFIGURACION: GRUPOS Y CARPETA OUTPUT ─────────────────────────

fprintf('=== PARTE II: Power Spectrum 1-120 Hz ===\n');

nGrpsStr2 = inputdlg('Cuantos grupos para el Power Spectrum?', 'PS — Grupos', 1, {'2'});
if isempty(nGrpsStr2), error('Abortado.'); end
nGroupsPS = str2double(nGrpsStr2{1});
assert(~isnan(nGroupsPS) && nGroupsPS >= 1, 'Numero de grupos no valido.');

psGroupLabels = cell(nGroupsPS, 1);
psGroupData   = cell(nGroupsPS, 1);

for g = 1:nGroupsPS
    lbl2 = inputdlg(sprintf('Etiqueta del grupo %d (debe coincidir con la que usa ps1to120):', g), ...
                    sprintf('PS Grupo %d', g), 1, {sprintf('G%d', g)});
    if isempty(lbl2), error('Abortado.'); end
    psGroupLabels{g} = strtrim(lbl2{1});

    fprintf('  Cargando datos para grupo "%s" — selecciona su carpeta cuando aparezca el dialogo...\n', psGroupLabels{g});
    psGroupData{g} = ps1to120_Alberto(psGroupLabels{g});
    fprintf('  Grupo "%s" cargado: %d animales.\n', psGroupLabels{g}, height(psGroupData{g}.idtable));
end

fprintf('>>> Selecciona carpeta OUTPUT para los resultados del Power Spectrum:\n');
psOutFolder = uigetdir(pwd, 'OUTPUT — Parte II Power Spectrum');
if psOutFolder == 0, error('Abortado.'); end

psFigFolder = fullfile(psOutFolder, 'figures', 'ps');
if ~exist(psFigFolder, 'dir'), mkdir(psFigFolder); end

freqs = psGroupData{1}.freqs;

%% ── P2.2  BUCLE PRINCIPAL: BANDAS × CONDICIONES ─────────────────────────

bandNames    = fieldnames(PS_BANDS);
groupSummary = table();
nL           = numel(LAYER_NAMES);

for bb = 1:numel(bandNames)

    bandName  = bandNames{bb};
    bandCfg   = PS_BANDS.(bandName);
    if ~bandCfg.active, continue; end

    freqIdx  = freqs >= bandCfg.range(1) & freqs <= bandCfg.range(2);
    freqVec  = freqs(freqIdx);
    freqVec  = freqVec(:);
    fprintf('\n  Banda: %s [%.1f-%.1f Hz] (%d bins)\n', ...
        bandName, bandCfg.range(1), bandCfg.range(2), sum(freqIdx));

    for cc = 1:size(PS_CONDITIONS, 1)

        condField  = PS_CONDITIONS{cc, 1};
        condLabel  = PS_CONDITIONS{cc, 2};
        condIsDB   = PS_CONDITIONS{cc, 3};

        peakFqRows   = [];
        powPeakRows  = [];
        meanPowRows  = [];
        groupCol     = {};
        idCol        = {};

        for g = 1:nGroupsPS

            gData  = psGroupData{g};
            gLabel = psGroupLabels{g};

            % --- 1) IDs de los registros (referencia para el numero de filas) ---
            ids   = getRecIDs(gData.idtable);
            nAnim = numel(ids);
            if nAnim == 0
                fprintf('    [AVISO] Grupo "%s" sin registros en idtable. Saltando.\n', gLabel);
                continue;
            end

            % --- 2) Matrices por capa, alineadas con ids ---
            %     Si la capa o el campo no existe (o esta vacio) -> todo NaN
            layerMats = cell(nL, 1);
            for L = 1:nL
                lname = LAYER_NAMES{L};

                hasData = isfield(gData, lname) && isstruct(gData.(lname)) && ...
                          isfield(gData.(lname), condField) && ...
                          ~isempty(gData.(lname).(condField));

                if hasData
                    matFull = gData.(lname).(condField);
                else
                    fprintf('    [AVISO] Grupo "%s": capa "%s" sin datos para "%s". Se rellena con NaN.\n', ...
                        gLabel, lname, condField);
                    matFull = nan(nAnim, numel(freqs));
                end

                nr = size(matFull, 1);
                if nr ~= nAnim
                    fprintf(['    [AVISO] Grupo "%s", capa "%s", "%s": %d filas vs %d registros en idtable. ' ...
                             'Se rellena/recorta por el final -> REVISA que el orden coincide.\n'], ...
                             gLabel, lname, condField, nr, nAnim);
                    if nr < nAnim
                        matFull(end+1:nAnim, :) = NaN;
                    else
                        matFull = matFull(1:nAnim, :);
                    end
                end

                if ~condIsDB
                    matFull = pow2db(matFull);
                end
                matFull(isinf(matFull)) = NaN;       % pow2db(0) = -Inf -> NaN
                layerMats{L} = matFull(:, freqIdx);
            end

            % --- 3) Metricas (los NaN se mantienen como hueco) ---
            peakFqBlock  = nan(nAnim, nL);
            powPeakBlock = nan(nAnim, nL);
            meanPowBlock = nan(nAnim, nL);

            for L = 1:nL
                mat    = layerMats{L};
                allNaN = all(isnan(mat), 2);

                [maxVals, maxIdx] = max(mat, [], 2, 'omitnan');
                pk = freqVec(maxIdx);
                pk = pk(:);

                maxVals(allNaN) = NaN;
                pk(allNaN)      = NaN;

                powPeakBlock(:, L) = maxVals;
                peakFqBlock(:, L)  = pk;
                meanPowBlock(:, L) = mean(mat, 2, 'omitnan');   % all-NaN -> NaN
            end

            % --- 4) Acumular ---
            peakFqRows  = [peakFqRows;  peakFqBlock];
            powPeakRows = [powPeakRows; powPeakBlock];
            meanPowRows = [meanPowRows; meanPowBlock];
            groupCol    = [groupCol;    repmat({gLabel}, nAnim, 1)];
            idCol       = [idCol;       ids(:)];
        end

        if isempty(peakFqRows), continue; end

        yUnitPow = 'dB';
        if ~condIsDB, yUnitPow = 'uV^2/Hz'; end

        % Tablas con ID (para Excel)
        idGrpTab = table(idCol, groupCol, 'VariableNames', ["ID","Group"]);
        tPeakFq  = [idGrpTab  array2table(peakFqRows,  'VariableNames', LAYER_NAMES)];
        tPowPeak = [idGrpTab  array2table(powPeakRows, 'VariableNames', LAYER_NAMES)];
        tMeanPow = [idGrpTab  array2table(meanPowRows, 'VariableNames', LAYER_NAMES)];

        % writetable deja la celda vacia para cada NaN numerico
        writetable(tPeakFq,  fullfile(psOutFolder, sprintf('ps_peakfq_%s_%s.xlsx',    bandName, condLabel)));
        writetable(tPowPeak, fullfile(psOutFolder, sprintf('ps_powtopeak_%s_%s.xlsx', bandName, condLabel)));
        writetable(tMeanPow, fullfile(psOutFolder, sprintf('ps_meanpow_%s_%s.xlsx',   bandName, condLabel)));

        % Resumen grupal (mean/SD con omitnan; n = todos los registros del grupo)
        for g = 1:nGroupsPS
            gMask = strcmp(tMeanPow.Group, psGroupLabels{g});
            vals  = tMeanPow(gMask, :);
            n     = sum(gMask);
            if n == 0, continue; end

            rowData = table( ...
                {bandName}, {condLabel}, {yUnitPow}, {psGroupLabels{g}}, {n}, ...
                mean(vals.or,  'omitnan'), std(vals.or,  'omitnan'), ...
                mean(vals.pyr, 'omitnan'), std(vals.pyr, 'omitnan'), ...
                mean(vals.rad, 'omitnan'), std(vals.rad, 'omitnan'), ...
                mean(vals.slm, 'omitnan'), std(vals.slm, 'omitnan'), ...
                'VariableNames', ["Band","Condition","Unit","Group","n", ...
                    "mean_or","sd_or","mean_pyr","sd_pyr", ...
                    "mean_rad","sd_rad","mean_slm","sd_slm"]);

            groupSummary = [groupSummary; rowData];
        end

        % Figuras: sin la columna ID (formato Group + capas, como antes)
        savePSFig(removevars(tPeakFq,  'ID'), sprintf('Peak Frequency — %s — %s', bandName, condLabel), ...
                  'Hz',     LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('peakfq_%s_%s',    bandName, condLabel));
        savePSFig(removevars(tPowPeak, 'ID'), sprintf('Power to Peak — %s — %s', bandName, condLabel), ...
                  yUnitPow, LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('powtopeak_%s_%s', bandName, condLabel));
        savePSFig(removevars(tMeanPow, 'ID'), sprintf('Mean Power — %s — %s', bandName, condLabel), ...
                  yUnitPow, LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('meanpow_%s_%s',   bandName, condLabel));

    end   % condiciones

    fprintf('    Banda "%s" completada.\n', bandName);
end   % bandas

if ~isempty(groupSummary)
    writetable(groupSummary, fullfile(psOutFolder, 'ps_group_summary.xlsx'));
    fprintf('  Resumen grupal guardado: ps_group_summary.xlsx\n');
end

fprintf('\n[Parte II] Completada. Resultados en: %s\n', psOutFolder);


%% ========================================================================
%  FUNCIONES LOCALES — ETAPA THETA (se llama 2 veces: pasada 1 y refresh)
%  runOnly = true  -> solo run (genera newintervals), sin nonrun/Excel/figs
%  runOnly = false -> run + nonrun completo, devuelve tablas y guarda figs
%  Los registros se procesan en parfor; figuras y prints van en serie despues.
%% ========================================================================
function T = runThetaStage(groupFolders, groupLabels, P, figFolder, runOnly)

nGroups = numel(groupFolders);
T.allrec = table();  T.run = table();  T.nonrun = table();  T.loc = table();

for g = 1:nGroups

    mainFolder = groupFolders{g};
    expCond    = groupLabels{g};

    entries  = dir(mainFolder);
    entries  = entries([entries.isdir]);
    entries  = entries(~ismember({entries.name}, {'.', '..'}));
    recNames = {entries.name};
    nRec     = numel(recNames);
    assert(nRec > 0, 'No se encontraron subcarpetas en: %s', mainFolder);
    fprintf('\n==== Grupo "%s": %d registro(s) (parfor) ====\n', expCond, nRec);

    outs = cell(nRec, 1);
    parfor ii = 1:nRec
        outs{ii} = processRecording(fullfile(mainFolder, recNames{ii}), recNames{ii}, P, runOnly);
    end

    zAllRec = nan(nRec, 3);   zRun = nan(nRec, 15);   zNonRun = nan(nRec, 15);
    locEvents = nan(nRec,1);  locTotalDur = nan(nRec,1);  locMeanDur = nan(nRec,1);
    locRecLen = nan(nRec,1);  locPercent = nan(nRec,1);

    for ii = 1:nRec
        o = outs{ii};
        fprintf('  [%d/%d] %s\n', ii, nRec, recNames{ii});
        for m = 1:numel(o.msg), fprintf('     %s\n', o.msg{m}); end

        if ~runOnly
            zAllRec(ii,:) = o.zAll;      zRun(ii,:) = o.zRun;      zNonRun(ii,:) = o.zNonRun;
            locEvents(ii) = o.nR;        locTotalDur(ii) = o.totDur;  locMeanDur(ii) = o.meanDur;
            locRecLen(ii) = o.recLen;    locPercent(ii)  = o.pct;
            try
                saveFigSession(recNames{ii}, expCond, o.runPow, o.nonrunPow, figFolder);
            catch ME
                warning('saveFigSession fallo en %s: %s', recNames{ii}, ME.message);
            end
            close all;
        end
    end

    if runOnly, continue; end

    idCol   = recNames(:);
    condCol = repmat({expCond}, nRec, 1);
    idTab   = table(idCol, condCol, 'VariableNames', ["ID","ExpCond"]);

    T.allrec  = [T.allrec;  idTab  array2table(zAllRec, 'VariableNames', ["thetapower (db)","lgammapower (db)","hgammapower (db)"])];
    T.run     = [T.run;     idTab  array2table(zRun,    'VariableNames', P.powerVarNames)];
    T.nonrun  = [T.nonrun;  idTab  array2table(zNonRun, 'VariableNames', P.powerVarNames)];

    locVars = ["ID","ExpCond","events(n)","totaldur (s)","meandur (s)","reclength (s)","timespent (%)"];
    T.loc = [T.loc; table(idCol, condCol, locEvents, locTotalDur, locMeanDur, locRecLen, locPercent, 'VariableNames', locVars)];

    fprintf('  Grupo "%s" completado.\n', expCond);
end
end


function o = processRecording(recPath, recName, P, runOnly)
% Procesa UN registro (corre dentro de un worker). Devuelve struct 'o'.

o = struct('zAll',nan(1,3), 'zRun',nan(1,15), 'zNonRun',nan(1,15), ...
           'nR',NaN, 'totDur',NaN, 'meanDur',NaN, 'recLen',NaN, 'pct',NaN, ...
           'runPow',[], 'nonrunPow',[]);
o.msg = {};

cd(recPath);
statesFile = fullfile(recPath, [recName '.thetaEpochs6to12.states.mat']);

% Version vieja de nonTheta_nonRipples (para comparar tras el refresh)
oldNLNR = [];
if ~runOnly && isfile(statesFile)
    tmpOld = load(statesFile, 'thisthetaEpochs');
    if isfield(tmpOld.thisthetaEpochs, 'intervals_nonTheta_nonRipples')
        oldNLNR = tmpOld.thisthetaEpochs.intervals_nonTheta_nonRipples;
    end
end

te         = detectThetaEpochs(P.DETECT_OPTS{:});
powts      = te.timestamps(:);
run_int    = te.thetaRun.intervals;
nonrun_int = te.intervals_nonTheta_nonRipples;

if ~runOnly
    if isempty(oldNLNR)
        o.msg{end+1} = sprintf('nonTheta_nonRipples: sin version previa para comparar (%d intervalos)', size(nonrun_int,1));
    elseif isequal(size(oldNLNR), size(nonrun_int)) && all(abs(oldNLNR(:) - nonrun_int(:)) < 1e-9)
        o.msg{end+1} = sprintf('nonTheta_nonRipples: IGUAL a la version previa (%d intervalos)', size(nonrun_int,1));
    else
        o.msg{end+1} = sprintf('nonTheta_nonRipples: ACTUALIZADO (%d -> %d intervalos)', size(oldNLNR,1), size(nonrun_int,1));
    end
end

S = load(fullfile(recPath, [recName '.hippocampalLayers.channelinfo.mat']), 'hippocampalLayers');
hippocampalLayers = S.hippocampalLayers;
nCh = hippocampalLayers.channels(end);

% Obtener SLM, fallback si falta o NaN
if isfield(hippocampalLayers.layers{1,1}, 'slm') && ~isempty(hippocampalLayers.layers{1,1}.slm)
    slmCh = hippocampalLayers.layers{1,1}.slm;
    if isnan(slmCh)
        o.msg{end+1} = sprintf('AVISO: SLM = NaN en %s. Usando ultimo canal (%d).', recName, nCh);
        slmCh = nCh;
    end
else
    o.msg{end+1} = sprintf('AVISO: SLM no encontrado en %s. Usando ultimo canal (%d).', recName, nCh);
    slmCh = nCh;
end
lfpT         = getLFP(slmCh, 'noPrompts', true);
samplingRate = lfpT.samplingRate;
lfpData      = double(lfpT.data(:,1));

run_int = run_int((run_int(:,2) - run_int(:,1)) >= P.MIN_RUN_DUR, :);

[tW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, P.THETA_BAND(1),  P.THETA_BAND(2),  P.N_CYCLES, 0);
[lW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, P.LGAMMA_BAND(1), P.LGAMMA_BAND(2), P.N_CYCLES, 0);
[hW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, P.HGAMMA_BAND(1), P.HGAMMA_BAND(2), P.N_CYCLES, 0);

rawT = median(tW,1);  rawL = median(lW,1);  rawH = median(hW,1);
tDB = pow2db(tW);  tDB_med = median(tDB,1);
lDB = pow2db(lW);  lDB_med = median(lDB,1);
hDB = pow2db(hW);  hDB_med = median(hDB,1);
tZ  = mean(zscore(tDB),1);  lZ = mean(zscore(lDB),1);  hZ = mean(zscore(hDB),1);
tZm = median(robust_zscore_matrix(tDB),1);
lZm = median(robust_zscore_matrix(lDB),1);
hZm = median(robust_zscore_matrix(hDB),1);

runPow = extractEpochPower(run_int, powts, rawT,rawL,rawH, tDB_med,lDB_med,hDB_med, tZ,lZ,hZ, tZm,lZm,hZm);
te     = storeRunPower(te, run_int, runPow);

nR   = size(run_int,1);
rDur = run_int(:,2) - run_int(:,1);

if ~runOnly
    nonrunPow = extractEpochPower(nonrun_int, powts, rawT,rawL,rawH, tDB_med,lDB_med,hDB_med, tZ,lZ,hZ, tZm,lZm,hZm);
    te        = storeNonRunPower(te, nonrun_int, nonrunPow);

    nNR  = size(nonrun_int,1);
    nDur = nonrun_int(:,2) - nonrun_int(:,1);

    o.zAll    = [median(tDB_med) median(lDB_med) median(hDB_med)];
    o.zRun    = buildRow(nR,  rDur, runPow);
    o.zNonRun = buildRow(nNR, nDur, nonrunPow);
    o.nR      = nR;
    o.totDur  = sum(rDur);
    o.meanDur = mean(rDur);
    o.recLen  = te.thetaRun.timestamps(end);
    o.pct     = o.totDur / o.recLen;
    o.runPow    = runPow;
    o.nonrunPow = nonrunPow;
end

saveStates(statesFile, te);
o.msg{end+1} = 'OK';
end


function saveStates(statesFile, thisthetaEpochs) %#ok<INUSD>
% save() no se puede llamar directamente dentro de parfor; aqui si.
save(statesFile, 'thisthetaEpochs');
end

function ids = getRecIDs(idtable)
% Devuelve cellstr (n x 1) con el nombre de cada registro a partir de idtable.
    vn  = idtable.Properties.VariableNames;
    pat = {'id','name','rec','sess','folder','animal'};
    col = 1;
    for k = 1:numel(pat)
        hit = find(contains(lower(vn), pat{k}), 1);
        if ~isempty(hit), col = hit; break; end
    end
    x   = idtable.(vn{col});
    ids = cellstr(string(x));
    ids = ids(:);
end

% %% ========================================================================
% %  LFP ANALYSIS PIPELINE — Post-HippoCookbook preprocessing
% %  ========================================================================
% %
% %  ── PARTE I: Theta & Gamma durante locomotion ────────────────────────────
% %    Detecta epocas theta, extrae potencia wavelet (theta/lgamma/hgamma)
% %    durante run y non-run. Guarda .mat por sesion + tablas Excel + figuras.
% %
% %  ── PARTE II: Power Spectrum 1-120 Hz (peak frequency + power to peak) ───
% %    Llama a ps1to120() por grupo para calcular espectros. Luego extrae
% %    peak frequency y potencia al pico para cada banda, capa y condicion.
% %    Guarda tablas Excel + figuras comparativas por banda.
% %    Se puede ejecutar de forma INDEPENDIENTE (Ctrl+Enter en su seccion).
% %
% %  ── PARTE III: Ripple Analysis ───────────────────────────────────────────
% %    Detecta SWRs con CNN (rippleMasterDetector), categoriza por
% %    sink (RAD/SLM/Other) via CSD, extrae medidas (amplitud, frecuencia,
% %    duracion, entropia, etc.), y guarda:
% %      · <rec>.ripples.events.mat  → para usar con nonlocnonrip / ps1to120
% %      · allrippledata.mat         → estructuras completas por grupo
% %      · tablas Excel de medidas filtradas
% %      · waveforms raw + filtered
% %
% %  ── PARTE IV: Análisis de delta (ARCHIVO SEPARADO: analyzedelta.m) ──────
% %    Calcula el Power Spectral Density (PSD) en la banda delta [1-3.5 Hz]
% %    por grupo experimental y capa hippocampal
% %
% %  INSTRUCCIONES:
% %    1. Edita el bloque P3.0 si necesitas cambiar paths CNN / thresholds
% %    2. F5 = corre todo | Ctrl+Enter seccion a seccion
% %    3. Los dialogos de carpetas aparecen al inicio de cada parte
% %
% %  SALIDA Parte I  (carpeta output elegida al inicio):
% %    theta_allrec.xlsx, theta_run.xlsx, theta_nonrun.xlsx
% %    locevent_table.xlsx
% %    figures/  ->  PNG por sesion + group_comparison_dB/zmod.png/.fig
% %
% %  SALIDA Parte II (carpeta output elegida al inicio de la parte):
% %    ps_peakfq_<banda>_<condicion>.xlsx
% %    ps_powtopeak_<banda>_<condicion>.xlsx
% %    ps_meanpow_<banda>_<condicion>.xlsx
% %    ps_group_summary.xlsx
% %    figures/ps/  ->  figura comparativa por banda y condicion
% %
% %  SALIDA Parte III (carpeta output elegida en P3.1):
% %    <sesion>/<sesion>.ripples.events.mat  — eventos por sesion (nonlocnonrip)
% %    allrippledata.mat                     — estructuras completas
% %    <grupo>_rippledata_keep.xlsx          — medidas filtradas por grupo
% %    all_ripples_combined.xlsx             — tabla combinada todos los grupos
% %    ripplewaveform.mat                    — formas de onda raw + filtered
% %    Figuras: histogramas (SVG)
% %
% %  SALIDA Parte IV (ver analyzedelta.m para detalles)
% % =========================================================================
% 
% clearvars; close all; clc;
% 
% %% ========================================================================
% %  PARTE I — THETA & GAMMA DURANTE LOCOMOTION
% %% ========================================================================
% 
% %% ── P1.0  PARAMETROS ─────────────────────────────────────────────────────
% 
% THETA_BAND  = [4  12];
% LGAMMA_BAND = [20 60];
% HGAMMA_BAND = [60 100];
% N_CYCLES    = 8;
% MIN_RUN_DUR = 1;   % segundos — epocas de carrera mas cortas se descartan
% 
% DETECT_OPTS = {'saveMat', false, 'force', true, ...
%                'bandpass', [6 12], 'delta_bandpass', [1 3.5], ...
%                'powerThreshold', 1};
% 
% %% ── P1.1  CONFIGURACION: GRUPOS Y CARPETAS ───────────────────────────────
% 
% nGroupsStr = inputdlg('Cuantos grupos experimentales tienes?', 'Grupos', 1, {'1'});
% if isempty(nGroupsStr), error('Abortado.'); end
% nGroups = str2double(nGroupsStr{1});
% assert(~isnan(nGroups) && nGroups >= 1, 'Numero de grupos no valido.');
% 
% groupFolders = cell(nGroups, 1);
% groupLabels  = cell(nGroups, 1);
% 
% for g = 1:nGroups
%     lbl = inputdlg(sprintf('Etiqueta del grupo %d (ej: PV, SST, WT...):', g), ...
%                    sprintf('Grupo %d', g), 1, {sprintf('G%d', g)});
%     if isempty(lbl), error('Abortado.'); end
%     groupLabels{g} = strtrim(lbl{1});
% 
%     fprintf('>>> Selecciona la carpeta del grupo "%s":\n', groupLabels{g});
%     folder = uigetdir(pwd, sprintf('Carpeta del grupo: %s', groupLabels{g}));
%     if folder == 0, error('No se selecciono carpeta. Abortado.'); end
%     groupFolders{g} = folder;
% end
% 
% fprintf('>>> Selecciona la carpeta OUTPUT para los resultados de la Parte I:\n');
% outFolder = uigetdir(groupFolders{1}, 'OUTPUT — Parte I');
% if outFolder == 0, error('Abortado.'); end
% 
% figFolder = fullfile(outFolder, 'figures');
% if ~exist(figFolder, 'dir'), mkdir(figFolder); end
% 
% %% ── P1.2  NOMBRES DE COLUMNAS ────────────────────────────────────────────
% 
% powerVarNames = ["events (n)", "totaldur (s)", "meandur (s)", ...
%     "medianTheta (uv)",  "medianLgamma (uv)",  "medianHgamma (uv)", ...
%     "medianTheta (db)",  "medianLgamma (db)",  "medianHgamma (db)", ...
%     "zTheta (au)",       "zLgamma (au)",       "zHgamma (au)",      ...
%     "zmodTheta (au)",    "zmodLgamma (au)",    "zmodHgamma (au)"];
% 
% %% ── P1.3  BUCLE POR GRUPO ────────────────────────────────────────────────
% 
% tableAllRec_all = table();
% tableRun_all    = table();
% tableNonRun_all = table();
% tableLoc_all    = table();
% 
% for g = 1:nGroups
% 
%     mainFolder = groupFolders{g};
%     expCond    = groupLabels{g};
% 
%     entries  = dir(mainFolder);
%     entries  = entries([entries.isdir]);
%     entries  = entries(~ismember({entries.name}, {'.', '..'}));
%     recNames = {entries.name};
%     nRec     = numel(recNames);
%     assert(nRec > 0, 'No se encontraron subcarpetas en: %s', mainFolder);
%     fprintf('\n==== Grupo "%s": %d registro(s) ====\n', expCond, nRec);
% 
%     zAllRec     = nan(nRec, 3);
%     zRun        = nan(nRec, 15);
%     zNonRun     = nan(nRec, 15);
%     locEvents   = nan(nRec, 1);
%     locTotalDur = nan(nRec, 1);
%     locMeanDur  = nan(nRec, 1);
%     locRecLen   = nan(nRec, 1);
%     locPercent  = nan(nRec, 1);
%     thetaEpochs = cell(nRec, 1);
% 
%     for ii = 1:nRec
% 
%         recName = recNames{ii};
%         recPath = fullfile(mainFolder, recName);
%         fprintf('  [%d/%d] %s\n', ii, nRec, recName);
% 
%         cd(recPath);
%         thetaEpochs{ii} = detectThetaEpochs(DETECT_OPTS{:});
%         powts           = thetaEpochs{ii}.timestamps(:);
%         run_int         = thetaEpochs{ii}.thetaRun.intervals;
%         nonrun_int      = thetaEpochs{ii}.intervals_nonTheta_nonRipples;
% 
%         layersFile = fullfile(recPath, [recName '.hippocampalLayers.channelinfo.mat']);
%         load(layersFile, 'hippocampalLayers');
%         %slmCh        = hippocampalLayers.layers{1,1}.slm;
%         nCh = hippocampalLayers.channels(end);   % porque es un 1xN double
% 
%         % Obtener SLM, fallback si falta o NaN
%         if isfield(hippocampalLayers.layers{1,1}, 'slm') && ...
%                 ~isempty(hippocampalLayers.layers{1,1}.slm)
% 
%             slmCh = hippocampalLayers.layers{1,1}.slm;
% 
%             if isnan(slmCh)
%                 warning('SLM = NaN en %s. Usando último canal (%d).', recName, nCh);
%                 slmCh = nCh;
%             end
% 
%         else
%             warning('SLM no encontrado en %s. Usando último canal (%d).', recName, nCh);
%             slmCh = nCh;
%         end
%         lfpT         = getLFP(slmCh, 'noPrompts', true);
%         samplingRate = lfpT.samplingRate;
%         lfpData      = double(lfpT.data(:,1));
% 
%         run_int = run_int((run_int(:,2) - run_int(:,1)) >= MIN_RUN_DUR, :);
% 
%         [tW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, THETA_BAND(1),  THETA_BAND(2),  N_CYCLES, 0);
%         [lW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, LGAMMA_BAND(1), LGAMMA_BAND(2), N_CYCLES, 0);
%         [hW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, HGAMMA_BAND(1), HGAMMA_BAND(2), N_CYCLES, 0);
% 
%         rawT = median(tW,1);  rawL = median(lW,1);  rawH = median(hW,1);
%         tDB = pow2db(tW);  tDB_med = median(tDB,1);
%         lDB = pow2db(lW);  lDB_med = median(lDB,1);
%         hDB = pow2db(hW);  hDB_med = median(hDB,1);
%         tZ  = mean(zscore(tDB),1);  lZ = mean(zscore(lDB),1);  hZ = mean(zscore(hDB),1);
%         tZm = median(robust_zscore_matrix(tDB),1);
%         lZm = median(robust_zscore_matrix(lDB),1);
%         hZm = median(robust_zscore_matrix(hDB),1);
% 
%         runPow    = extractEpochPower(run_int,    powts, rawT,rawL,rawH, tDB_med,lDB_med,hDB_med, tZ,lZ,hZ, tZm,lZm,hZm);
%         nonrunPow = extractEpochPower(nonrun_int, powts, rawT,rawL,rawH, tDB_med,lDB_med,hDB_med, tZ,lZ,hZ, tZm,lZm,hZm);
% 
%         thetaEpochs{ii} = storeRunPower(thetaEpochs{ii},    run_int,    runPow);
%         thetaEpochs{ii} = storeNonRunPower(thetaEpochs{ii}, nonrun_int, nonrunPow);
% 
%         nR   = size(run_int,1);    nNR  = size(nonrun_int,1);
%         rDur = run_int(:,2)    - run_int(:,1);
%         nDur = nonrun_int(:,2) - nonrun_int(:,1);
% 
%         zAllRec(ii,:) = [median(tDB_med) median(lDB_med) median(hDB_med)];
%         zRun(ii,:)    = buildRow(nR,  rDur, runPow);
%         zNonRun(ii,:) = buildRow(nNR, nDur, nonrunPow);
% 
%         locEvents(ii)   = nR;
%         locTotalDur(ii) = sum(rDur);
%         locMeanDur(ii)  = mean(rDur);
%         locRecLen(ii)   = thetaEpochs{ii}.thetaRun.timestamps(end);
%         locPercent(ii)  = locTotalDur(ii) / locRecLen(ii);
% 
%         thisthetaEpochs = thetaEpochs{ii};
%         save(fullfile(recPath, [recName '.thetaEpochs6to12.states.mat']), 'thisthetaEpochs');
%         saveFigSession(recName, expCond, runPow, nonrunPow, figFolder);
% 
%         clearvars -except ...
%             g nGroups groupFolders groupLabels ...
%             ii nRec recNames mainFolder outFolder figFolder expCond ...
%             THETA_BAND LGAMMA_BAND HGAMMA_BAND N_CYCLES MIN_RUN_DUR DETECT_OPTS ...
%             powerVarNames thetaEpochs ...
%             zAllRec zRun zNonRun ...
%             locEvents locTotalDur locMeanDur locRecLen locPercent ...
%             tableAllRec_all tableRun_all tableNonRun_all tableLoc_all;
%         close all;
%         fprintf('     OK\n');
%     end
% 
%     idCol   = recNames(:);
%     condCol = repmat({expCond}, nRec, 1);
%     idTab   = table(idCol, condCol, 'VariableNames', ["ID","ExpCond"]);
% 
%     tableAllRec_all = [tableAllRec_all; idTab  array2table(zAllRec,  'VariableNames', ["thetapower (db)","lgammapower (db)","hgammapower (db)"])];
%     tableRun_all    = [tableRun_all;    idTab  array2table(zRun,     'VariableNames', powerVarNames)];
%     tableNonRun_all = [tableNonRun_all; idTab  array2table(zNonRun,  'VariableNames', powerVarNames)];
% 
%     locVars = ["ID","ExpCond","events(n)","totaldur (s)","meandur (s)","reclength (s)","timespent (%)"];
%     tableLoc_all = [tableLoc_all; table(idCol, condCol, locEvents, locTotalDur, locMeanDur, locRecLen, locPercent, 'VariableNames', locVars)];
% 
%     fprintf('  Grupo "%s" completado.\n', expCond);
% end
% 
% %% ── P1.4  GUARDAR EXCEL ──────────────────────────────────────────────────
% 
% writetable(tableAllRec_all, fullfile(outFolder, 'theta_allrec.xlsx'));
% writetable(tableRun_all,    fullfile(outFolder, 'theta_run.xlsx'));
% writetable(tableNonRun_all, fullfile(outFolder, 'theta_nonrun.xlsx'));
% writetable(tableLoc_all,    fullfile(outFolder, 'locevent_table.xlsx'));
% fprintf('\n[Parte I] Excel guardados en: %s\n', outFolder);
% 
% %% ── P1.5  FIGURAS GRUPALES ───────────────────────────────────────────────
% 
% saveGroupFig(tableRun_all, tableNonRun_all, groupLabels, figFolder);
% fprintf('[Parte I] Completada.\n\n');
% 
% 
% 
% %% ========================================================================
% %  PARTE III — RIPPLE ANALYSIS
% %% ========================================================================
% 
% %% ── P3.0  PARAMETROS ─────────────────────────────────────────────────────
% 
% % Rutas del modelo CNN y ejecutable Python
% CNN_MODEL_PATH = 'C:\Users\usuario\cnn-matlab-master\cnn';
% CNN_PYTHON_EXE = 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe';
% 
% % Parametros de deteccion
% RIPPLE_THRESHOLD = [2 5];        % [min max] en SD
% RIPPLE_PASSBAND  = [80 240];     % Hz
% PEAKNORM_CUTOFF  = 50;           % outliers con PeakNormPow > este valor se descartan
% 
% % Analisis robusto de outliers — thresholds
% OUTLIER_THRESH_MAD  = 3.5;
% OUTLIER_MIN_N_REC   = 20;
% OUTLIER_MIN_N_SUB   = 20;
% OUTLIER_MIN_N_GROUP = 20;
% OUTLIER_VARS        = {'PeakAmp','PeakNormPow','Duration','PeakFq'};
% OUTLIER_LOG_VARS    = {'PeakAmp','PeakNormPow','Duration'};
% 
% %% ── P3.1  CONFIGURACION: GRUPOS Y CARPETAS ───────────────────────────────
% %
% %  NOTA: si acabas de ejecutar la Parte I en esta misma sesion de MATLAB,
% %  las variables nGroups / groupFolders / groupLabels ya existen. El dialogo
% %  te permite reutilizar los mismos grupos o cambiarlos si quieres analizar
% %  solo un subconjunto.
% 
% nGroupsStr3 = inputdlg('Cuantos grupos para el analisis de ripples?', 'Ripples — Grupos', 1, {num2str(nGroups)});
% if isempty(nGroupsStr3), error('Abortado.'); end
% nGroupsR = str2double(nGroupsStr3{1});
% assert(~isnan(nGroupsR) && nGroupsR >= 1, 'Numero de grupos no valido.');
% 
% ripGroupFolders = cell(nGroupsR, 1);
% ripGroupLabels  = cell(nGroupsR, 1);
% 
% for g = 1:nGroupsR
%     defaultLbl = '';
%     if g <= nGroups, defaultLbl = groupLabels{g}; else, defaultLbl = sprintf('G%d', g); end
%     lbl3 = inputdlg(sprintf('Etiqueta del grupo %d:', g), ...
%                     sprintf('Ripples Grupo %d', g), 1, {defaultLbl});
%     if isempty(lbl3), error('Abortado.'); end
%     ripGroupLabels{g} = strtrim(lbl3{1});
% 
%     fprintf('>>> Selecciona la carpeta del grupo "%s":\n', ripGroupLabels{g});
%     folder3 = uigetdir(pwd, sprintf('Carpeta del grupo: %s', ripGroupLabels{g}));
%     if folder3 == 0, error('No se selecciono carpeta. Abortado.'); end
%     ripGroupFolders{g} = folder3;
% end
% 
% fprintf('>>> Selecciona la carpeta OUTPUT para los resultados de Ripples:\n');
% rippleFolder = uigetdir(ripGroupFolders{1}, 'OUTPUT — Parte III Ripples');
% if rippleFolder == 0, error('Abortado.'); end
% 
% %% ── P3.2  LISTAR SUBCARPETAS POR GRUPO ───────────────────────────────────
% 
% lfpPaths = cell(nGroupsR, 1);
% 
% for g = 1:nGroupsR
%     folderInfo = dir(ripGroupFolders{g});
%     subfolders = folderInfo([folderInfo.isdir]);
%     subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
%     lfpPaths{g} = {subfolders.name};
%     fprintf('  Grupo "%s": %d sesion(es) encontradas.\n', ripGroupLabels{g}, numel(lfpPaths{g}));
% end
% 
% %% ── P3.3  DETECCION Y CLASIFICACION DE RIPPLES (POR GRUPO) ──────────────
% 
% allRipplesStruct = cell(nGroupsR, 1);
% allCSD           = cell(nGroupsR, 1);
% allCSDEvents     = cell(nGroupsR, 1);
% allCSDSummary    = cell(nGroupsR, 1);
% sinkSizes        = struct();
% 
% for g = 1:nGroupsR
%     grpLabel = ripGroupLabels{g};
%     paths    = lfpPaths{g};
%     nSess    = numel(paths);
% 
%     fprintf('\n==== Procesando grupo "%s": %d sesion(es) ====\n', grpLabel, nSess);
% 
%     ripplesCell = cell(nSess, 1);
%     csdCell     = cell(nSess, 1);
%     csdEvCell   = cell(nSess, 1);
%     csdSumCell  = cell(nSess, 1);
% 
%     slmSize   = nan(nSess, 1);
%     radSize   = nan(nSess, 1);
%     otherSize = nan(nSess, 1);
% 
%     % Evita "Graphics timeout" con figuras invisibles durante deteccion batch
%     set(0, 'DefaultFigureVisible', 'off');
% 
%     for ii = 1:nSess
%         thisFolderName = fullfile(ripGroupFolders{g}, paths{ii});
%         fprintf('  [%d/%d] %s\n', ii, nSess, paths{ii});
%         cd(thisFolderName);
% 
%         % Cargar canales hippocampales
%         layersFile = fullfile(thisFolderName, [paths{ii} '.hippocampalLayers.channelinfo.mat']);
%         load(layersFile, 'hippocampalLayers');
%         rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
%         radCh         = hippocampalLayers.layers{1,1}.radiatum;
%         SWChannel     = hippocampalLayers.layers{1,1}.radiatum;
%         slmCh         = hippocampalLayers.layers{1,1}.slm;
% 
%         % Cargar intervalos theta (generados en Parte I)
%         thetaFile = fullfile(thisFolderName, [paths{ii} '.thetaEpochs6to12.states.mat']);
%         load(thetaFile, 'thisthetaEpochs');
%         thetaIntervals = thisthetaEpochs.thetaRun.newintervals;
% 
%         % Cargar LFP
%         lfp = getLFP('all');
% 
%         try
%             % Deteccion con CNN
%             ripples = rippleMasterDetector('rippleChannel', rippleChannel, 'SWChannel', SWChannel, ...
%                 'detector', 'customcnn', 'exec_env', CNN_PYTHON_EXE, 'model_file', CNN_MODEL_PATH, ...
%                 'passband', RIPPLE_PASSBAND, 'force', true, 'skipStimulationPeriods', true, ...
%                 'thresholds', RIPPLE_THRESHOLD, 'eventSpikeThreshold', false, ...
%                 'excludeIntervals', thetaIntervals, ...
%                 'rippleStats', true);
% 
%             try
%                 ripples = computeRippleStats_new( ...
%                     'ripples',       ripples, ...
%                     'rippleChannel', rippleChannel, ...
%                     'bandpass',      RIPPLE_PASSBAND, ...
%                     'thetaEpochs',   thisthetaEpochs, ...
%                     'saveSummary',   false, ...
%                     'saveMat',       false, ...
%                     'plotOpt',       false);
%                 fprintf('    computeRippleStats_new OK: %d ripples\n', numel(ripples.peaks));
%             catch ME_stats
%                 warning('computeRippleStats_new fallo en sesion %s:\n  %s\nSesion descartada.', ...
%                     paths{ii}, ME_stats.message);
%                 fprintf('    Stack trace:\n');
%                 for sti = 1:numel(ME_stats.stack)
%                     fprintf('      %s (linea %d)\n', ME_stats.stack(sti).name, ME_stats.stack(sti).line);
%                 end
%                 ripplesCell{ii} = []; csdCell{ii} = []; csdEvCell{ii} = []; csdSumCell{ii} = [];
%                 close all; drawnow;
%                 continue;
%             end
% 
%             % Fix valid_events: sincronizar rippleStats.data al numero de peaks
%             nRipClean = numel(ripples.peaks);
%             nRipStats = size(ripples.rippleStats.data.peakAmplitude, 1);
%             if nRipStats ~= nRipClean
%                 fprintf('    FIX valid_events: rippleStats.data tenia %d filas vs %d peaks. Recortando...\n', nRipStats, nRipClean);
%                 statFields = fieldnames(ripples.rippleStats.data);
%                 for sf = 1:numel(statFields)
%                     fdata = ripples.rippleStats.data.(statFields{sf});
%                     if isnumeric(fdata) && size(fdata,1) == nRipStats
%                         ripples.rippleStats.data.(statFields{sf}) = fdata(1:nRipClean, :);
%                     end
%                 end
%             end
% 
%             % ── GUARDAR <rec>.ripples.events.mat ─────────────────────────
%             %    Formato compatible con nonlocnonrip / ps1to120.
%             %    Contiene: timestamps [n x 2], peaks [n x 1], y campo
%             %    detectorinfo para que buzcode lo reconozca como events file.
%             ripplesEvents.timestamps  = ripples.timestamps;
%             ripplesEvents.peaks       = ripples.peaks;
%             ripplesEvents.amplitude   = ripples.peakNormedPower;
%             ripplesEvents.duration    = diff(ripples.timestamps, 1, 2);
%             ripplesEvents.detectorinfo.detectorname = 'rippleMasterDetector';
%             ripplesEvents.detectorinfo.detectionparms.passband   = RIPPLE_PASSBAND;
%             ripplesEvents.detectorinfo.detectionparms.thresholds = RIPPLE_THRESHOLD;
%             ripplesEventsFile = fullfile(thisFolderName, [paths{ii} '.ripples.events.mat']);
%             save(ripplesEventsFile, 'ripplesEvents', '-mat');
%             fprintf('    Guardado: %s.ripples.events.mat\n', paths{ii});
% 
%             % CSD por ripple
%             [csd, ~, csdEvents, ~] = buzcode_eacheventCSD(lfp, ripples.peaks, ...
%                 'plotCSD', false, 'plotLFP', false);
% 
%             % Resumir CSD por ripple (media temporal centrada)
%             nRip = size(csdEvents, 3);
%             csdSum.tsRipple    = ripples.timestamps;
%             csdSum.durRipple   = diff(csdSum.tsRipple, 1, 2) * 1000;  % ms
%             csdSum.rangeRipple = round(csdSum.durRipple / 2);
% 
%             srCSD = csd.samplingRate;
%             nCSDchan = size(csdEvents, 2);
%             csdSum.csdSummary = nan(nRip, nCSDchan);
% 
%             for jj = 1:nRip
%                 halfDur_samp = max(round(csdSum.rangeRipple(jj) / 1000 * srCSD), 1);
%                 idxWin = csd.timestamps >= -halfDur_samp & csd.timestamps <= halfDur_samp;
%                 if ~any(idxWin)
%                     idxWin = true(size(csd.timestamps));  % fallback: ventana completa
%                 end
%                 slice   = csdEvents(idxWin, :, jj);
%                 rowMean = mean(slice, 1);
%                 csdSum.csdSummary(jj, :) = rowMean;
%             end
% 
%             % Categorizar por sink (RAD / SLM / Other)
%             [~, sinkIdx] = min(csdSum.csdSummary, [], 2);
%             slmPeaks   = (sinkIdx >= radCh + 1);
%             radPeaks   = (sinkIdx >= rippleChannel + 1) & (sinkIdx <= radCh);
%             otherPeaks = ~radPeaks & ~slmPeaks;
% 
%             slmSize(ii)   = sum(slmPeaks);
%             radSize(ii)   = sum(radPeaks);
%             otherSize(ii) = sum(otherPeaks);
% 
%             fprintf('    SLM: %d | RAD: %d | Other: %d\n', slmSize(ii), radSize(ii), otherSize(ii));
% 
%             ripples.rippleStats.data.radpeaksidx = radPeaks;
%             ripples.rippleStats.data.slmpeaksidx = slmPeaks;
%             ripples.rippleStats.data.others      = otherPeaks;
% 
%             ripplesCell{ii} = ripples;
%             csdCell{ii}     = csd;
%             csdEvCell{ii}   = csdEvents;
%             csdSumCell{ii}  = csdSum;
% 
%         catch ME
%             warning('Error en sesion %s: %s', paths{ii}, ME.message);
%             ripplesCell{ii} = [];
%             csdCell{ii}     = [];
%             csdEvCell{ii}   = [];
%             csdSumCell{ii}  = [];
%         end
%         close all; drawnow;
%     end
% 
%     % Restaurar visibilidad de figuras tras deteccion batch
%     set(0, 'DefaultFigureVisible', 'on');
% 
%     allRipplesStruct{g} = ripplesCell;
%     allCSD{g}           = csdCell;
%     allCSDEvents{g}     = csdEvCell;
%     allCSDSummary{g}    = csdSumCell;
% 
%     sinkSizes.(grpLabel).slm   = slmSize;
%     sinkSizes.(grpLabel).rad   = radSize;
%     sinkSizes.(grpLabel).other = otherSize;
% 
%     fprintf('  Grupo "%s" completado.\n', grpLabel);
% end
% 
% % Guardar estructuras completas de ripples
% allRipplesFile = fullfile(rippleFolder, 'allrippledata.mat');
% saveVars = cellfun(@(lbl) [lbl 'ripples'], ripGroupLabels, 'UniformOutput', false);
% for g = 1:nGroupsR
%     eval(sprintf('%s = allRipplesStruct{%d};', saveVars{g}, g));
% end
% save(allRipplesFile, saveVars{:}, '-mat');
% fprintf('\nEstructuras de ripples guardadas: %s\n', allRipplesFile);
% 
% %% ── P3.4  EXTRAER MEDIDAS Y ENSAMBLAR TABLAS ────────────────────────────
% 
% allRipplesTables = cell(nGroupsR, 1);
% allRipplesRaw    = cell(nGroupsR, 1);
% allRipplesFilt   = cell(nGroupsR, 1);
% 
% for g = 1:nGroupsR
%     grpLabel = ripGroupLabels{g};
%     paths    = lfpPaths{g};
%     ripplesCell = allRipplesStruct{g};
% 
%     measures  = [];
%     rawWaves  = [];
%     filtWaves = [];
% 
%     for ii = 1:numel(ripplesCell)
%         if isempty(ripplesCell{ii}), continue; end
%         ripples = ripplesCell{ii};
% 
%         % Extraer medidas
%         peakNormPow = ripples.peakNormedPower;
%         peakFq      = ripples.rippleStats.data.peakFrequency;
%         peakAmp     = ripples.rippleStats.data.peakAmplitude;
%         duration    = ripples.rippleStats.data.duration;
%         interRipFq  = ripples.rippleStats.data.interRippleFrequency;
%         specEnt     = ripples.rippleStats.data.spectralEntropy(:);
%         fastRipIdx  = ripples.rippleStats.data.fastRippleIndex(:);
%         multitap    = ripples.rippleStats.data.multiTapperFreq;
% 
%         rippleIFstat    = ripples.rippleStats.data.rippleIFstat;
%         rippleIF_meanFq = ripples.rippleStats.data.rippleIF.meanFreq;
%         rippleIF_sdFq   = ripples.rippleStats.data.rippleIF.stdFreq;
%         rippleIF_drift  = ripples.rippleStats.data.rippleIF.drift;
% 
%         % Sink type
%         sinkType = strings(numel(peakNormPow), 1);
%         sinkType(ripples.rippleStats.data.radpeaksidx) = "rad";
%         sinkType(ripples.rippleStats.data.slmpeaksidx) = "slm";
%         sinkType(ripples.rippleStats.data.others)      = "other";
% 
%         thisID  = repmat(string(paths{ii}), numel(peakNormPow), 1);
%         expCond = repmat(string(grpLabel),  numel(peakNormPow), 1);
% 
%         thisTable = table(thisID, expCond, peakNormPow, peakFq, peakAmp, duration, interRipFq, ...
%             specEnt, fastRipIdx, multitap, rippleIFstat, rippleIF_meanFq, rippleIF_sdFq, rippleIF_drift, sinkType, ...
%             'VariableNames', ["ID","expcond","PeakNormPow","PeakFq","PeakAmp","Duration","interrippleFq", ...
%             "spectralEntropy","FastRippleIdx","multitap","freqstat","meanfq_nonstat","std_nonstat","drift_nonstat","sinktype"]);
% 
%         measures = [measures; thisTable];
% 
%         rawWaves  = [rawWaves;  array2table(ripples.rippleStats.maps.ripples_raw,      'VariableNames', string(ripples.rippleStats.maps.timestamps))];
%         filtWaves = [filtWaves; array2table(ripples.rippleStats.maps.ripples_filtered, 'VariableNames', string(ripples.rippleStats.maps.timestamps))];
%     end
% 
%     allRipplesTables{g} = measures;
%     allRipplesRaw{g}    = rawWaves;
%     allRipplesFilt{g}   = filtWaves;
% end
% 
% % Combinar todos los grupos
% allRipples     = vertcat(allRipplesTables{:});
% allRipplesRaw  = vertcat(allRipplesRaw{:});
% allRipplesFilt = vertcat(allRipplesFilt{:});
% 
% % Guardar tabla original
% writetable(allRipples, fullfile(rippleFolder, 'all_ripples_original.xlsx'));
% 
% % Filtrar outliers por PeakNormPow > PEAKNORM_CUTOFF
% keepIdx = allRipples.PeakNormPow < PEAKNORM_CUTOFF;
% allRipplesKeep = allRipples(keepIdx, :);
% 
% % Guardar por grupo
% for g = 1:nGroupsR
%     grpIdx = strcmp(allRipples.expcond, ripGroupLabels{g});
%     grpKeepIdx = grpIdx & keepIdx;
%     writetable(allRipples(grpIdx, :),     fullfile(rippleFolder, sprintf('%s_rippledata_original.xlsx', ripGroupLabels{g})));
%     writetable(allRipples(grpKeepIdx, :), fullfile(rippleFolder, sprintf('%s_rippledata_keep.xlsx',     ripGroupLabels{g})));
% end
% 
% % Guardar tabla combinada filtrada
% writetable(allRipplesKeep, fullfile(rippleFolder, 'all_ripples_keep.xlsx'));
% 
% % Guardar waveforms
% waveformFile = fullfile(rippleFolder, 'ripplewaveform.mat');
% save(waveformFile, 'allRipplesRaw', 'allRipplesFilt', '-mat');
% 
% fprintf('\nTablas de ripples guardadas en: %s\n', rippleFolder);
% 
% %% ── P3.5  AÑADIR COLUMNA SUBJECT ─────────────────────────────────────────
% 
% if ~ismember('subject', allRipplesKeep.Properties.VariableNames)
%     % Extrae ID de animal del nombre de sesion: "C089_190520_sess1" -> "C089"
%     allRipplesKeep.subject = extractBefore(allRipplesKeep.ID, '_');
%     % Si no hay patron, usar ID directamente:
%     % allRipplesKeep.subject = allRipplesKeep.ID;
% end
% 
% %% ── P3.6  OUTLIER REMOVAL (JERARQUICO ROBUSTO) ──────────────────────────
% 
% data = allRipplesKeep;
% fprintf('\n=== Removing outliers (hierarchical robust MAD) ===\n');
% 
% missingVars = setdiff(OUTLIER_VARS, data.Properties.VariableNames);
% if ~isempty(missingVars)
%     error('Columnas no encontradas en la tabla: %s\nColumnas disponibles: %s', ...
%         strjoin(missingVars, ', '), strjoin(data.Properties.VariableNames, ', '));
% end
% 
% n = height(data);
% flagsMat = false(n, numel(OUTLIER_VARS));
% 
% uniqueRecs   = unique(data.ID);
% uniqueSubs   = unique(data.subject);
% uniqueGroups = unique(data.expcond);
% 
% for vi = 1:numel(OUTLIER_VARS)
%     var  = OUTLIER_VARS{vi};
%     xAll = data.(var);
%     flagsVar = false(n, 1);
%     isLog = ismember(var, OUTLIER_LOG_VARS);
% 
%     for r = 1:numel(uniqueRecs)
%         recID  = uniqueRecs(r);
%         idxRec = strcmp(data.ID, recID);
%         valid  = idxRec & ~isnan(xAll);
%         if isLog, valid = valid & (xAll > 0); end
%         x = xAll(valid);
%         if isempty(x), continue; end
% 
%         if numel(x) >= OUTLIER_MIN_N_REC
%             xproc = x; validUse = valid;
%         else
%             subjID  = unique(data.subject(idxRec));
%             idxSub  = strcmp(data.subject, subjID);
%             validSub = idxSub & ~isnan(xAll);
%             if isLog, validSub = validSub & (xAll > 0); end
%             xSub = xAll(validSub);
%             if numel(xSub) >= OUTLIER_MIN_N_SUB
%                 xproc = xSub; validUse = validSub;
%             else
%                 grpID   = unique(data.expcond(idxRec));
%                 idxGrp  = strcmp(data.expcond, grpID);
%                 validGrp = idxGrp & ~isnan(xAll);
%                 if isLog, validGrp = validGrp & (xAll > 0); end
%                 xGrp = xAll(validGrp);
%                 if numel(xGrp) >= OUTLIER_MIN_N_GROUP
%                     xproc = xGrp; validUse = validGrp;
%                 else
%                     validGlob = ~isnan(xAll);
%                     if isLog, validGlob = validGlob & (xAll > 0); end
%                     xGlob = xAll(validGlob);
%                     if isempty(xGlob), continue; end
%                     xproc = xGlob; validUse = validGlob;
%                 end
%             end
%         end
% 
%         if isLog, xproc = log10(xproc(xproc>0)); end
%         med = median(xproc);
%         madv = median(abs(xproc - med));
%         scaledMAD = 1.4826 * madv;
%         if scaledMAD == 0, scaledMAD = max(std(xproc), eps); end
% 
%         idxValidList = find(validUse);
%         xvalid = xAll(validUse);
%         if isLog, xvalid = log10(xvalid); end
%         z = (xvalid - med) / scaledMAD;
%         outIdx = idxValidList(abs(z) > OUTLIER_THRESH_MAD);
%         flagsVar(outIdx) = true;
%     end
% 
%     flagsMat(:, vi) = flagsVar;
%     fprintf('  %s: %d outliers\n', var, sum(flagsVar));
% end
% 
% unionIdx = any(flagsMat, 2);
% fprintf('TOTAL outliers (union): %d\n', sum(unionIdx));
% 
% data = data(~unionIdx, :);
% fprintf('Data after outlier removal: %d rows\n', height(data));
% 
% %% ── P3.7  HISTOGRAMAS POR GRUPO ──────────────────────────────────────────
% 
% vars = {'PeakFq','PeakAmp','spectralEntropy','Duration'};
% 
% if nGroupsR <= 3
%     cols = [1.00 0.55 0.75; 0.65 0.85 1.00; 0.05 0.25 0.65];
%     cols = cols(1:nGroupsR, :);
% else
%     cols = lines(nGroupsR);
% end
% 
% alphaVal = 0.35;
% nbins    = 30;
% 
% fprintf('\n=== Diagnostico de variables para histogramas ===\n');
% for i = 1:numel(vars)
%     varname = vars{i};
%     if ismember(varname, data.Properties.VariableNames)
%         nvalid = sum(~isnan(data.(varname)) & isfinite(data.(varname)));
%         fprintf('  %s: %d valores validos\n', varname, nvalid);
%     else
%         fprintf('  %s: *** COLUMNA NO ENCONTRADA ***\n', varname);
%         fprintf('    Columnas disponibles: %s\n', strjoin(data.Properties.VariableNames, ', '));
%     end
% end
% fprintf('=================================================\n\n');
% 
% fig = figure('Units','normalized','Position',[0.05 0.15 0.9 0.6],'Color','w');
% for i = 1:numel(vars)
%     ax = subplot(1, 4, i);
%     hold(ax, 'on');
%     varname = vars{i};
% 
%     if ~ismember(varname, data.Properties.VariableNames)
%         title(ax, [strrep(varname,'_',' ') ' (NOT FOUND)']);
%         continue;
%     end
% 
%     combined = data.(varname);
%     combined = combined(~isnan(combined) & isfinite(combined));
%     if isempty(combined), title(ax, strrep(varname,'_',' ')); continue; end
%     edges = linspace(min(combined), max(combined), nbins+1);
% 
%     plotted   = gobjects(0);
%     legLabels = {};
%     for g = 1:nGroupsR
%         grp = ripGroupLabels{g};
%         idx = strcmp(data.expcond, grp);
%         grpdata = data.(varname)(idx);
%         grpdata = grpdata(~isnan(grpdata) & isfinite(grpdata));
%         if isempty(grpdata), continue; end
%         h = histogram(ax, grpdata, edges, 'Normalization','probability', ...
%             'FaceAlpha', alphaVal, 'EdgeColor','none', 'FaceColor', cols(g,:));
%         plotted(end+1) = h;
%         legLabels{end+1} = grp;
%     end
% 
%     title(ax, strrep(varname,'_',' '), 'Interpreter','none');
%     xlabel(ax,'Value');
%     if i == 1, ylabel(ax,'Probability'); end
%     if ~isempty(plotted), legend(ax, plotted, legLabels, 'Location','best'); end
%     set(ax,'FontSize',11,'Box','off');
%     hold(ax,'off');
% end
% 
% saveas(fig, fullfile(rippleFolder, 'all_ripples_histograms.png'));
% close(fig);
% 
% fprintf('\n[Parte III] Completada. Archivos .ripples.events.mat guardados en cada carpeta de sesion.\n');
% fprintf('           Resultados de analisis en: %s\n\n', rippleFolder);
% 
% %% ========================================================================
% %  PARTE II — POWER SPECTRUM 1-120 Hz
% %  (puede ejecutarse de forma independiente con Ctrl+Enter)
% %% ========================================================================
% 
% %% ── P2.0  PARAMETROS DE BANDAS ───────────────────────────────────────────
% %  Desactiva una banda poniendo 'active', false
% 
% PS_BANDS = struct( ...
%     'delta',       struct('active', true, 'range', [1    3.5]), ...
%     'theta',       struct('active', true, 'range', [6    12]),  ...
%     'beta',        struct('active', true, 'range', [12   30]),  ...
%     'lowgamma',    struct('active', true, 'range', [30   60]),  ...
%     'mediumgamma', struct('active', true, 'range', [60   90]),  ...
%     'highgamma',   struct('active', true, 'range', [90   120]));
% 
% % Capas hippocampales — orden fijo que devuelve ps1to120()
% LAYER_NAMES = {'or', 'pyr', 'rad', 'slm'};
% 
% % Condiciones — campos dentro de cada capa en la estructura de ps1to120()
% % Formato: {nombre_campo, etiqueta_archivos, es_db}
% PS_CONDITIONS = {
%     'allrec',              'AllRec_dB',       true;
%     'allrecinuv',          'AllRec_uV',       true;
%     'locindb',             'Loc_dB',          true;
%     'locinuv',             'Loc_uV',          true;
%     'nonlocnonrippleindb', 'NonLocNonRip_dB', true;
%     'nonlocnonrippleinuv', 'NonLocNonRip_uV', true;
%     'nonlocindb',          'NonLoc_dB',       true;
%     'nonlocinuv',          'NonLoc_uV',       true;
% };
% 
% %% ── P2.1  CONFIGURACION: GRUPOS Y CARPETA OUTPUT ─────────────────────────
% 
% fprintf('=== PARTE II: Power Spectrum 1-120 Hz ===\n');
% 
% nGrpsStr2 = inputdlg('Cuantos grupos para el Power Spectrum?', 'PS — Grupos', 1, {'2'});
% if isempty(nGrpsStr2), error('Abortado.'); end
% nGroupsPS = str2double(nGrpsStr2{1});
% assert(~isnan(nGroupsPS) && nGroupsPS >= 1, 'Numero de grupos no valido.');
% 
% psGroupLabels = cell(nGroupsPS, 1);
% psGroupData   = cell(nGroupsPS, 1);
% 
% % ps1to120() abre su propio uigetdir internamente — solo necesitamos la etiqueta
% for g = 1:nGroupsPS
%     lbl2 = inputdlg(sprintf('Etiqueta del grupo %d (debe coincidir con la que usa ps1to120):', g), ...
%                     sprintf('PS Grupo %d', g), 1, {sprintf('G%d', g)});
%     if isempty(lbl2), error('Abortado.'); end
%     psGroupLabels{g} = strtrim(lbl2{1});
% 
%     fprintf('  Cargando datos para grupo "%s" — selecciona su carpeta cuando aparezca el dialogo...\n', psGroupLabels{g});
%     psGroupData{g} = ps1to120_Alberto(psGroupLabels{g});
%     fprintf('  Grupo "%s" cargado: %d animales.\n', psGroupLabels{g}, height(psGroupData{g}.idtable));
% end
% 
% fprintf('>>> Selecciona carpeta OUTPUT para los resultados del Power Spectrum:\n');
% psOutFolder = uigetdir(pwd, 'OUTPUT — Parte II Power Spectrum');
% if psOutFolder == 0, error('Abortado.'); end
% 
% psFigFolder = fullfile(psOutFolder, 'figures', 'ps');
% if ~exist(psFigFolder, 'dir'), mkdir(psFigFolder); end
% 
% freqs = psGroupData{1}.freqs;   % vector de frecuencias (identico entre grupos)
% 
% %% ── P2.2  BUCLE PRINCIPAL: BANDAS × CONDICIONES ─────────────────────────
% 
% bandNames = fieldnames(PS_BANDS);
% 
% % Accumulator for group summary (mean+SD across animals, all bands x conditions)
% groupSummary = table();
% 
% for bb = 1:numel(bandNames)
% 
%     bandName  = bandNames{bb};
%     bandCfg   = PS_BANDS.(bandName);
%     if ~bandCfg.active, continue; end
% 
%     freqIdx  = freqs >= bandCfg.range(1) & freqs <= bandCfg.range(2);
%     freqVec  = freqs(freqIdx);
%     fprintf('\n  Banda: %s [%.1f-%.1f Hz] (%d bins)\n', ...
%         bandName, bandCfg.range(1), bandCfg.range(2), sum(freqIdx));
% 
%     for cc = 1:size(PS_CONDITIONS, 1)
% 
%         condField  = PS_CONDITIONS{cc, 1};   % nombre del campo en la estructura
%         condLabel  = PS_CONDITIONS{cc, 2};   % etiqueta para archivos y titulos
%         condIsDB   = PS_CONDITIONS{cc, 3};   % true si ya esta en dB
% 
%         peakFqRows   = [];
%         powPeakRows  = [];
%         meanPowRows  = [];   % media de potencia sobre la banda (nueva)
%         groupCol     = {};
% 
%         for g = 1:nGroupsPS
% 
%             gData  = psGroupData{g};
%             gLabel = psGroupLabels{g};
% 
%             % --- 1) Verificar que el campo exista en TODAS las capas ---
%             fieldMissing = false;
%             for L = 1:numel(LAYER_NAMES)
%                 if ~isfield(gData.(LAYER_NAMES{L}), condField)
%                     fprintf('    [AVISO] Campo "%s" no encontrado en capa "%s" del grupo "%s". Saltando grupo.\n', ...
%                         condField, LAYER_NAMES{L}, gLabel);
%                     fieldMissing = true;
%                     break;
%                 end
%             end
%             if fieldMissing, continue; end
% 
% 
%             % --- 2) Determinar número real de animales por capa ---
%             nAnimLayer = zeros(numel(LAYER_NAMES),1);
%             for L = 1:numel(LAYER_NAMES)
%                 mat = gData.(LAYER_NAMES{L}).(condField);
%                 nAnimLayer(L) = size(mat,1);
%             end
% 
%             % nAnim = mínimo común denominador entre capas
%             nAnim = min(nAnimLayer);
% 
%             if nAnim == 0
%                 fprintf('    [AVISO] Grupo "%s" tiene 0 animales con datos para "%s". Saltando.\n', gLabel, condField);
%                 continue;
%             end
% 
%             % --- 3) Inicializar máscara válida ---
%             validMask = true(nAnim,1);
% 
%             % --- 4) Extraer matrices recortadas y limpiar NANs ---
%             layerMats = cell(numel(LAYER_NAMES), 1);
% 
%             for L = 1:numel(LAYER_NAMES)
%                 matFull = gData.(LAYER_NAMES{L}).(condField);
% 
%                 % recortar filas a nAnim
%                 matFull = matFull(1:nAnim, :);
% 
%                 % convertir a dB si es necesario
%                 if ~condIsDB
%                     matFull = pow2db(matFull);
%                 end
% 
%                 % aplicar máscara de frecuencias
%                 layerMats{L} = matFull(:, freqIdx);
% 
%                 % marcar animales inválidos si tienen NaNs
%                 validMask = validMask & ~any(isnan(layerMats{L}), 2);
%             end
% 
%             nValid = sum(validMask);
%             if nValid == 0
%                 fprintf('    [AVISO] Grupo "%s", condicion "%s": sin datos válidos tras filtrar. Saltando.\n', gLabel, condLabel);
%                 continue;
%             end
% 
%             % --- 5) Calcular métricas ---
%             peakFqBlock   = nan(nValid, numel(LAYER_NAMES));
%             powPeakBlock  = nan(nValid, numel(LAYER_NAMES));
%             meanPowBlock  = nan(nValid, numel(LAYER_NAMES));
% 
%             for L = 1:numel(LAYER_NAMES)
% 
%                 mat = layerMats{L}(validMask, :);  % <-- AHORA SIEMPRE CONSISTENTE
% 
%                 [maxVals, maxIdx] = max(mat, [], 2);
% 
%                 powPeakBlock(:, L) = maxVals;
%                 peakFqBlock(:, L)  = freqVec(maxIdx(:)');
%                 meanPowBlock(:, L) = mean(mat, 2, 'omitnan');
%             end
% 
%             % --- 6) Acumular resultados ---
%             peakFqRows  = [peakFqRows;  peakFqBlock];
%             powPeakRows = [powPeakRows; powPeakBlock];
%             meanPowRows = [meanPowRows; meanPowBlock];
%             groupCol    = [groupCol;    repmat({gLabel}, nValid, 1)];
% 
%         end
% 
%         if isempty(peakFqRows), continue; end
% 
%         % Unidad de potencia para titulos y figuras
%         yUnitPow = 'dB';
%         if ~condIsDB, yUnitPow = 'uV^2/Hz'; end
% 
%         % Montar tablas
%         grpTab     = table(groupCol, 'VariableNames', "Group");
%         tPeakFq    = [grpTab  array2table(peakFqRows,  'VariableNames', LAYER_NAMES)];
%         tPowPeak   = [grpTab  array2table(powPeakRows, 'VariableNames', LAYER_NAMES)];
%         tMeanPow   = [grpTab  array2table(meanPowRows, 'VariableNames', LAYER_NAMES)];
% 
%         % ── Guardar Excel por animal ─────────────────────────────────────────
%         writetable(tPeakFq,  fullfile(psOutFolder, sprintf('ps_peakfq_%s_%s.xlsx',   bandName, condLabel)));
%         writetable(tPowPeak, fullfile(psOutFolder, sprintf('ps_powtopeak_%s_%s.xlsx', bandName, condLabel)));
%         writetable(tMeanPow, fullfile(psOutFolder, sprintf('ps_meanpow_%s_%s.xlsx',   bandName, condLabel)));
% 
%         % ── Resumen grupal: media ± SD entre animales ─────────────────────
%         for g = 1:nGroupsPS
%             gMask = strcmp(tMeanPow.Group, psGroupLabels{g});
%             vals  = tMeanPow(gMask, :);
%             n     = sum(gMask);
%             if n == 0, continue; end
% 
%             rowData = table( ...
%                 {bandName}, {condLabel}, {yUnitPow}, {psGroupLabels{g}}, {n}, ...
%                 mean(vals.or,  'omitnan'), std(vals.or,  'omitnan'), ...
%                 mean(vals.pyr, 'omitnan'), std(vals.pyr, 'omitnan'), ...
%                 mean(vals.rad, 'omitnan'), std(vals.rad, 'omitnan'), ...
%                 mean(vals.slm, 'omitnan'), std(vals.slm, 'omitnan'), ...
%                 'VariableNames', ["Band","Condition","Unit","Group","n", ...
%                     "mean_or","sd_or","mean_pyr","sd_pyr", ...
%                     "mean_rad","sd_rad","mean_slm","sd_slm"]);
% 
%             groupSummary = [groupSummary; rowData];
%         end
% 
%         % ── Figuras ───────────────────────────────────────────────────────
%         savePSFig(tPeakFq,  sprintf('Peak Frequency — %s — %s', bandName, condLabel), ...
%                   'Hz',       LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('peakfq_%s_%s',   bandName, condLabel));
%         savePSFig(tPowPeak, sprintf('Power to Peak — %s — %s', bandName, condLabel), ...
%                   yUnitPow,   LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('powtopeak_%s_%s', bandName, condLabel));
%         savePSFig(tMeanPow, sprintf('Mean Power — %s — %s', bandName, condLabel), ...
%                   yUnitPow,   LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('meanpow_%s_%s',   bandName, condLabel));
% 
%     end   % condiciones
% 
%     fprintf('    Banda "%s" completada.\n', bandName);
% end   % bandas
% 
% % ── Guardar tabla resumen grupal (una sola vez, todas las bandas) ─────────
% if ~isempty(groupSummary)
%     writetable(groupSummary, fullfile(psOutFolder, 'ps_group_summary.xlsx'));
%     fprintf('  Resumen grupal guardado: ps_group_summary.xlsx\n');
% end
% 
% fprintf('\n[Parte II] Completada. Resultados en: %s\n', psOutFolder);
