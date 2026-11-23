%% ========================================================================
%  RIPPLE ANALYSIS PIPELINE — Deteccion, clasificacion y analisis UMAP
%  ========================================================================
%
%  DESCRIPCION:
%    Detecta sharp-wave ripples usando CNN (rippleMasterDetector), los
%    categoriza por sink (RAD/SLM/Other) usando CSD, extrae medidas
%    (amplitud, frecuencia, duracion, entropia, etc.), y realiza analisis
%    UMAP con visualizacion y LME por grupo.
%
%  REQUISITOS:
%    - rippleMasterDetector en el path (con modelo CNN)
%    - buzcode_eacheventCSD, buzcode_eventCSD en el path
%    - run_umap (UMAP implementation)
%    - fitlme (Statistics and Machine Learning Toolbox)
%    - Archivos .thetaEpochs6to12.states.mat por sesion (de Parte I)
%    - Archivos .hippocampalLayers.channelinfo.mat por sesion
%
%  INSTRUCCIONES:
%    1. Edita bloque P3.0 (paths de CNN, thresholds, bandas)
%    2. F5 o Ctrl+Enter por seccion
%    3. Dialogos: cuantos grupos + etiqueta y carpeta por grupo + carpeta output
%    4. Todo se guarda automaticamente
%
%  SALIDA (carpeta output elegida):
%    allrippledata.mat            — estructuras de ripples por grupo
%    <grupo>_rippledata_keep.xlsx — medidas de ripples filtradas
%    all_ripples_combined.xlsx    — tabla combinada todos los grupos
%    ripplewaveform.mat           — formas de onda raw + filtered
%    Figuras: histogramas, UMAP, boxplots, etc. (SVG)
% =========================================================================

clearvars; close all; clc;

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

nGroupsStr = inputdlg('Cuantos grupos experimentales tienes?', 'Grupos', 1, {'3'});
if isempty(nGroupsStr), error('Abortado.'); end
nGroups = str2double(nGroupsStr{1});
assert(~isnan(nGroups) && nGroups >= 1, 'Numero de grupos no valido.');

groupFolders = cell(nGroups, 1);
groupLabels  = cell(nGroups, 1);

for g = 1:nGroups
    lbl = inputdlg(sprintf('Etiqueta del grupo %d (ej: SS, OS, OD, WT, KO...):', g), ...
                   sprintf('Grupo %d', g), 1, {sprintf('G%d', g)});
    if isempty(lbl), error('Abortado.'); end
    groupLabels{g} = strtrim(lbl{1});

    fprintf('>>> Selecciona la carpeta del grupo "%s":\n', groupLabels{g});
    folder = uigetdir(pwd, sprintf('Carpeta del grupo: %s', groupLabels{g}));
    if folder == 0, error('No se selecciono carpeta. Abortado.'); end
    groupFolders{g} = folder;
end

fprintf('>>> Selecciona la carpeta OUTPUT para los resultados:\n');
rippleFolder = uigetdir(groupFolders{1}, 'OUTPUT — Ripple Analysis');
if rippleFolder == 0, error('Abortado.'); end

%% ── P3.2  LISTAR SUBCARPETAS POR GRUPO ───────────────────────────────────

lfpPaths = cell(nGroups, 1);

for g = 1:nGroups
    folderInfo = dir(groupFolders{g});
    subfolders = folderInfo([folderInfo.isdir]);
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
    lfpPaths{g} = {subfolders.name};
    fprintf('  Grupo "%s": %d sesion(es) encontradas.\n', groupLabels{g}, numel(lfpPaths{g}));
end

%% ── P3.3  DETECCION Y CLASIFICACION DE RIPPLES (POR GRUPO) ──────────────

allRipplesStruct = cell(nGroups, 1);
allCSD           = cell(nGroups, 1);
allCSDEvents     = cell(nGroups, 1);
allCSDSummary    = cell(nGroups, 1);
sinkSizes        = struct();

for g = 1:nGroups
    grpLabel = groupLabels{g};
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

    % ── FIX: figuras invisibles durante deteccion batch ──────────────────
    % Evita "Graphics timeout" que ocurre cuando rippleMasterDetector llama
    % a print/saveas con figuras visibles acumuladas en el EDT de Java.
    set(0, 'DefaultFigureVisible', 'off');

    for ii = 1:nSess
        thisFolderName = fullfile(groupFolders{g}, paths{ii});
        fprintf('  [%d/%d] %s\n', ii, nSess, paths{ii});
        cd(thisFolderName);

        % Cargar canales hippocampales
        layersFile = fullfile(thisFolderName, [paths{ii} '.hippocampalLayers.channelinfo.mat']);
        load(layersFile, 'hippocampalLayers');
        rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
        radCh         = hippocampalLayers.layers{1,1}.radiatum;
        SWChannel     = hippocampalLayers.layers{1,1}.radiatum;
        slmCh         = hippocampalLayers.layers{1,1}.slm;

        % Cargar intervalos theta (para excluir)
        thetaFile = fullfile(thisFolderName, [paths{ii} '.thetaEpochs6to12.states.mat']);
        load(thetaFile, 'thisthetaEpochs');
        thetaIntervals = thisthetaEpochs.thetaRun.newintervals;

        % Cargar LFP
        lfp = getLFP('all');

        try
            % Deteccion con CNN
            % ── NOTA: con detector='customcnn', exec_env y model_file ────
            % estan HARDCODEADOS dentro de rippleMasterDetector (rama customcnn
            % llama a findRipples_CNN_modpeak con rutas fijas). Los parametros
            % CNN_PYTHON_EXE y CNN_MODEL_PATH de este script son ignorados para
            % ese detector. Si quieres usar las rutas de P3.0, cambia el
            % detector a 'cnn' o edita la rama customcnn en rippleMasterDetector.
            % ── rippleStats=false: evita que computeRippleStats_new crashee ──
            % internamente en rippleMasterDetector sin mensaje de error util.
            % Lo llamamos nosotros a continuacion en un try/catch propio.
            ripples = rippleMasterDetector_Alberto('rippleChannel', rippleChannel, 'SWChannel', SWChannel, ...
                'detector', 'customcnn', 'exec_env', CNN_PYTHON_EXE, 'model_file', CNN_MODEL_PATH, ...
                'passband', RIPPLE_PASSBAND, 'force', true, 'skipStimulationPeriods', true, ...
                'thresholds', RIPPLE_THRESHOLD, 'eventSpikeThreshold', false, ...
                'excludeIntervals', thetaIntervals, ...
                'rippleStats', false);

            % ── Llamar a computeRippleStats_new con manejo de error propio ──
            %
            % BUGS conocidos en computeRippleStats_new que afectan a la llamada:
            %
            % 1. addParameter(p,'basepath',pwd,@isstruct) — el validador es
            %    @isstruct pero pwd devuelve un char. Si pasamos basepath
            %    explicitamente como char falla la validacion. NO pasar basepath:
            %    la funcion usa pwd internamente, y como ya hicimos cd(thisFolderName)
            %    antes, pwd == thisFolderName. Correcto.
            %
            % 2. saveSummary=true por defecto: llama a mkdir+saveas aunque no
            %    queramos. Pasar saveSummary=false explicitamente.
            %
            % 3. plotOpt=[]: la figura "h = figure" se crea SIEMPRE (fuera del
            %    bloque if plotOpt). Solo se cierra si ~plotOpt. Como ya tenemos
            %    DefaultFigureVisible=off, las figuras se crean en memoria pero
            %    no se renderizan. Pasar plotOpt=false para que se cierren.
            %
            % 4. saveMat=false por defecto en la funcion, pero lo pasamos
            %    explicitamente para evitar escrituras no deseadas en batch.
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
                continue;   % saltar al siguiente ii
            end

            % ── FIX: valid_events dentro de rippleMasterDetector recorta ─────
            % timestamps/peaks/peakNormedPower pero NO rippleStats.data.
            % Si hay desajuste, sincronizamos rippleStats.data al numero de peaks.
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

            % CSD por ripple
            % outputs: [csd, lfpAvg, csd_events, lfp_events]
            % csd_events: [tiempo x (nChan-2) x nRipples]  <- nChan-2 por diff(,2,2)
            % csd.timestamps: en MUESTRAS (no ms)
            % plotCSD/plotLFP=false para evitar figuras en batch
            [csd, ~, csdEvents, ~] = buzcode_eacheventCSD(lfp, ripples.peaks, ...
                'plotCSD', false, 'plotLFP', false);

            % Resumir CSD por ripple (media temporal centrada)
            nRip = size(csdEvents, 3);
            csdSum.tsRipple    = ripples.timestamps;
            csdSum.durRipple   = diff(csdSum.tsRipple, 1, 2) * 1000;  % ms
            csdSum.rangeRipple = round(csdSum.durRipple / 2);

            % ── csd_events: [tiempo x (nChan-2) x nRipples] ─────────────
            % csd.timestamps esta en MUESTRAS (ej: -125:125 a 1250 Hz),
            % NO en ms. rangeRipple esta en ms -> hay que convertir a muestras
            % antes de usarlo como indice de ventana.
            srCSD = csd.samplingRate;            % Hz
            nCSDchan = size(csdEvents, 2);       % nChan-2 por diff(,2,2)
            csdSum.csdSummary = nan(nRip, nCSDchan);

            for jj = 1:nRip
                % Convertir duracion del ripple (ms) a muestras
                halfDur_samp = max(round(csdSum.rangeRipple(jj) / 1000 * srCSD), 1);
                idxWin = csd.timestamps >= -halfDur_samp & csd.timestamps <= halfDur_samp;
                if ~any(idxWin)
                    idxWin = true(size(csd.timestamps));  % fallback: ventana completa
                end
                slice   = csdEvents(idxWin, :, jj);      % [k x nCSDchan]
                rowMean = mean(slice, 1);                  % [1 x nCSDchan]
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
        % ── FIX: cerrar todas las figuras acumuladas (incluidas las de ──
        % rippleMasterDetector) y dar tiempo al EDT para liberarlas
        close all; drawnow;
    end

    % ── FIX: restaurar visibilidad de figuras tras deteccion batch ───────
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

% Guardar estructuras de ripples
allRipplesFile = fullfile(rippleFolder, 'allrippledata.mat');
saveVars = cellfun(@(lbl) [lbl 'ripples'], groupLabels, 'UniformOutput', false);
for g = 1:nGroups
    eval(sprintf('%s = allRipplesStruct{%d};', saveVars{g}, g));
end
save(allRipplesFile, saveVars{:}, '-mat');
fprintf('\nEstructuras de ripples guardadas: %s\n', allRipplesFile);

%% ── P3.4  EXTRAER MEDIDAS Y ENSAMBLAR TABLAS ────────────────────────────

allRipplesTables = cell(nGroups, 1);
allRipplesRaw    = cell(nGroups, 1);
allRipplesFilt   = cell(nGroups, 1);

for g = 1:nGroups
    grpLabel = groupLabels{g};
    paths    = lfpPaths{g};
    ripplesCell = allRipplesStruct{g};

    measures = [];
    rawWaves = [];
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

        % Tabla
        % ── FIX: usar string() en vez de repmat(char,...) para evitar ────
        % error de vertcat cuando los labels tienen distinta longitud.
        % repmat de char crea matriz nxLen_label; distintas longitudes
        % hacen que vertcat falle con "dimensions not consistent".
        thisID  = repmat(string(paths{ii}), numel(peakNormPow), 1);
        expCond = repmat(string(grpLabel),  numel(peakNormPow), 1);

        thisTable = table(thisID, expCond, peakNormPow, peakFq, peakAmp, duration, interRipFq, ...
            specEnt, fastRipIdx, multitap, rippleIFstat, rippleIF_meanFq, rippleIF_sdFq, rippleIF_drift, sinkType, ...
            'VariableNames', ["ID","expcond","PeakNormPow","PeakFq","PeakAmp","Duration","interrippleFq", ...
            "spectralEntropy","FastRippleIdx","multitap","freqstat","meanfq_nonstat","std_nonstat","drift_nonstat","sinktype"]);

        measures = [measures; thisTable];

        % Waveforms
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

% Guardar originales
writetable(allRipples, fullfile(rippleFolder, 'all_ripples_original.xlsx'));

% Filtrar outliers por PeakNormPow > PEAKNORM_CUTOFF
keepIdx = allRipples.PeakNormPow < PEAKNORM_CUTOFF;
allRipplesKeep = allRipples(keepIdx, :);

% Guardar por grupo
for g = 1:nGroups
    grpIdx = strcmp(allRipples.expcond, groupLabels{g});
    grpKeepIdx = grpIdx & keepIdx;
    writetable(allRipples(grpIdx, :),    fullfile(rippleFolder, sprintf('%s_rippledata_original.xlsx', groupLabels{g})));
    writetable(allRipples(grpKeepIdx, :), fullfile(rippleFolder, sprintf('%s_rippledata_keep.xlsx',     groupLabels{g})));
end

% Guardar tabla combinada filtrada
writetable(allRipplesKeep, fullfile(rippleFolder, 'all_ripples_keep.xlsx'));

% Guardar waveforms
waveformFile = fullfile(rippleFolder, 'ripplewaveform.mat');
save(waveformFile, 'allRipplesRaw', 'allRipplesFilt', '-mat');

fprintf('\nTablas de ripples guardadas en: %s\n', rippleFolder);

%% ── P3.5  AÑADIR COLUMNA SUBJECT (OPCIONAL) ──────────────────────────────
% Si tu dataset tiene un mapeo de sesiones a sujetos (animales), añadelo aqui.
% Ejemplo: extraer ID de animal del nombre de sesion.
% De lo contrario, usar la sesion misma como subject para el LME.

if ~ismember('subject', allRipplesKeep.Properties.VariableNames)
    % Crear subject desde ID (asume que ID contiene identificador de animal)
    % Ej: si ID = "C089_190520_sess1", subject = "C089"
    allRipplesKeep.subject = extractBefore(allRipplesKeep.ID, '_');
    % Si no hay patron, simplemente usar ID como subject:
    % allRipplesKeep.subject = allRipplesKeep.ID;
end

%% ── P3.6  OUTLIER REMOVAL (JERARQUICO ROBUSTO) ──────────────────────────

data = allRipplesKeep;
fprintf('\n=== Removing outliers (hierarchical robust MAD) ===\n');

% ── FIX: verificar que las columnas esperadas existen antes de proceder ──
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

        % Elegir nivel jerarquico segun N
        if numel(x) >= OUTLIER_MIN_N_REC
            xproc = x; validUse = valid;
        else
            % Subir a subject
            subjID  = unique(data.subject(idxRec));
            idxSub  = strcmp(data.subject, subjID);
            validSub = idxSub & ~isnan(xAll);
            if isLog, validSub = validSub & (xAll > 0); end
            xSub = xAll(validSub);
            if numel(xSub) >= OUTLIER_MIN_N_SUB
                xproc = xSub; validUse = validSub;
            else
                % Subir a group
                grpID   = unique(data.expcond(idxRec));
                idxGrp  = strcmp(data.expcond, grpID);
                validGrp = idxGrp & ~isnan(xAll);
                if isLog, validGrp = validGrp & (xAll > 0); end
                xGrp = xAll(validGrp);
                if numel(xGrp) >= OUTLIER_MIN_N_GROUP
                    xproc = xGrp; validUse = validGrp;
                else
                    % Global
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

% ── FIX: definir colores para cualquier numero de grupos ──────────────────
if nGroups <= 3
    cols = [1.00 0.55 0.75; 0.65 0.85 1.00; 0.05 0.25 0.65];
    cols = cols(1:nGroups, :);
else
    cols = lines(nGroups);
end

alphaVal = 0.35;
nbins    = 30;

% ── FIX: diagnostico — verificar que las variables tienen datos validos ───
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

    % ── FIX: saltar si la columna no existe ──────────────────────────────
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
    for g = 1:nGroups
        grp = groupLabels{g};
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

saveas(fig, fullfile(rippleFolder, 'all_ripples_histograms.svg'));  % evita dependencia de Ghostscript
close(fig);

