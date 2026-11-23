%% ========================================================================
%  PARTE VI — CROSS-FREQUENCY COUPLING (acoplamiento fase theta – amplitud gamma)
%  ========================================================================
%
%  Para cada sesion, capa (or/pyr/rad/slm) y condicion (run / non-run):
%    · Modulation Index (Tort et al. 2010) fase-theta -> amplitud-gamma
%      (low / medium / high gamma, mismas bandas que la Parte II)
%    · Z-score y p-valor contra surrogates (desplazamiento circular de la amplitud)
%    · Fase preferida de acoplamiento (grados) y distribucion fase-amplitud
%    · Comodulograma (fase 4-14 Hz x amplitud 25-135 Hz) por capa y condicion
%
%  REQUISITOS:
%    · Haber corrido la Parte I (usa <rec>.thetaEpochs6to12.states.mat)
%    · <rec>.hippocampalLayers.channelinfo.mat en cada sesion
%    · buzcode (getLFP), Signal Processing Toolbox, Statistics Toolbox (boxplot)
%
%  SALIDA (carpeta elegida al inicio):
%    cfc_MI.xlsx          — tabla larga: ID, ExpCond, Layer, Condition, GammaBand,
%                           PhaseRef, Dur_s, MI, MI_z, MI_p, PrefPhase_deg
%    cfc_results.mat      — tabla + comodulogramas + distribuciones fase-amplitud
%    figures/cfc/         — boxplots MI, comodulogramas medios, fase-amplitud
%    + Excels para rehacer figuras en otras apps (ver bloque V.3):
%      cfc_phaseamp_*.xlsx (ondas fase-amplitud), cfc_comodulogram_*.xlsx (heatmaps),
%      cfc_MI_groupsummary.xlsx
%
%  PARALELO: parfor por sesion (Parallel Computing Toolbox). C.useParallel / C.nWorkers.
%
%  CONVENCION DE FASE: 0 grados = pico del theta (angulo de Hilbert), rango [-180,180].
%  Se puede ejecutar de forma INDEPENDIENTE (no hace clearvars).
% =========================================================================

%% ── V.0  PARAMETROS ──────────────────────────────────────────────────────

C.phaseBand    = [6 12];                         % banda de fase (igual que detectThetaEpochs)
C.gammaNames   = {'lowgamma','mediumgamma','highgamma'};
C.gammaBands   = {[30 60], [60 90], [90 120]};   % igual que Parte II
C.layerNames   = {'or','pyr','rad','slm'};
% nombres alternativos del campo en hippocampalLayers.layers{1,1} (se usa el primero que exista)
C.layerFields  = {{'oriens','or'}, {'pyramidal','pyr'}, {'radiatum','rad'}, {'slm'}};
C.condNames    = {'run','nonrun'};
C.phaseRefLayer = 'self';   % 'self' = fase y amplitud del mismo canal | 'slm' / 'rad' / 'pyr' / 'or' = fase de referencia fija
C.nBins        = 18;        % bins de fase (Tort)
C.nSurr        = 200;       % surrogates por medida
C.minShiftS    = 1;         % desplazamiento minimo del surrogate (s)
C.minRunDur    = 1;         % s — epocas de run mas cortas se descartan (igual que Parte I)
C.minTotalDur  = 10;        % s — duracion total minima por condicion para calcular MI
C.doComod      = true;      % false = salta comodulogramas (mucho mas rapido)
C.comodPhase   = 4:0.5:14;    C.comodPhaseBW = 2;    % centros y ancho de banda de fase (Hz)
C.comodAmp     = 25:2.5:135;  C.comodAmpBW   = 20;   % centros y ancho de banda de amplitud (Hz)
C.seed         = 1;
C.useParallel  = true;      % parfor por sesion
C.nWorkers     = 8;         % [] = todos los cores. Cada worker carga ~4 LFP completos: baja si te quedas sin RAM

rng(C.seed);

%% ── V.1  CONFIGURACION: GRUPOS Y CARPETAS ────────────────────────────────

cfcReuse = false;
if exist('groupFolders','var') && exist('groupLabels','var') && exist('nGroups','var')
    ans_ = questdlg(sprintf('Reutilizar los %d grupos ya definidos?', nGroups), 'CFC — Grupos', 'Si','No','Si');
    cfcReuse = strcmp(ans_, 'Si');
end

if cfcReuse
    cfcFolders = groupFolders;  cfcLabels = groupLabels;  nGroupsC = nGroups;
else
    s_ = inputdlg('Cuantos grupos para el analisis CFC?', 'CFC — Grupos', 1, {'2'});
    if isempty(s_), error('Abortado.'); end
    nGroupsC = str2double(s_{1});
    assert(~isnan(nGroupsC) && nGroupsC >= 1, 'Numero de grupos no valido.');
    cfcFolders = cell(nGroupsC,1);  cfcLabels = cell(nGroupsC,1);
    for g = 1:nGroupsC
        l_ = inputdlg(sprintf('Etiqueta del grupo %d:', g), sprintf('CFC Grupo %d', g), 1, {sprintf('G%d', g)});
        if isempty(l_), error('Abortado.'); end
        cfcLabels{g} = strtrim(l_{1});
        fprintf('>>> Selecciona la carpeta del grupo "%s":\n', cfcLabels{g});
        f_ = uigetdir(pwd, sprintf('Carpeta del grupo: %s', cfcLabels{g}));
        if f_ == 0, error('No se selecciono carpeta. Abortado.'); end
        cfcFolders{g} = f_;
    end
end

fprintf('>>> Selecciona la carpeta OUTPUT para los resultados de CFC:\n');
cfcOut = uigetdir(cfcFolders{1}, 'OUTPUT — Parte V CFC');
if cfcOut == 0, error('Abortado.'); end
cfcFigFolder = fullfile(cfcOut, 'figures', 'cfc');
if ~exist(cfcFigFolder, 'dir'), mkdir(cfcFigFolder); end

%% ── V.2  LISTA DE SESIONES + POOL PARALELO ───────────────────────────────

jobG = [];  jobName = {};  jobPath = {};
for g = 1:nGroupsC
    ent = dir(cfcFolders{g});
    ent = ent([ent.isdir]);
    ent = ent(~ismember({ent.name}, {'.','..'}));
    for ii = 1:numel(ent)
        jobG(end+1)    = g;                                      %#ok<SAGROW>
        jobName{end+1} = ent(ii).name;                           %#ok<SAGROW>
        jobPath{end+1} = fullfile(cfcFolders{g}, ent(ii).name);  %#ok<SAGROW>
    end
    fprintf('Grupo "%s": %d sesion(es)\n', cfcLabels{g}, sum(jobG == g));
end
nJobs    = numel(jobName);
jobLabel = cfcLabels(jobG);
assert(nJobs > 0, 'No se encontraron sesiones.');

if C.useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(C.nWorkers), pool = parpool('local'); else, pool = parpool('local', C.nWorkers); end
    end
    fprintf('Pool paralelo: %d workers\n', pool.NumWorkers);
end

dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('%s\n', msg));

results = cell(nJobs, 1);
tAll = tic;
parfor j = 1:nJobs
    tJ = tic;
    Sj = cfcSessionSafe(jobPath{j}, jobName{j}, jobLabel{j}, C, j);
    results{j} = Sj;
    if Sj.ok, st = 'OK'; else, st = ['FALLO: ' Sj.msg]; end
    send(dq, sprintf('  [%d/%d] %s — %s (%.0f s)', j, nJobs, jobName{j}, st, toc(tJ)));
end
fprintf('\nCFC completado en %.1f min\n', toc(tAll)/60);

cfcRows  = cell(0, 11);
sessList = {};
for j = 1:nJobs
    if ~results{j}.ok
        warning('CFC fallo en %s: %s', jobName{j}, results{j}.msg);
        continue;
    end
    cfcRows = [cfcRows; results{j}.rows];  %#ok<AGROW>
    sessList{end+1} = results{j};          %#ok<SAGROW>
end
clear results;
assert(~isempty(cfcRows), 'No se obtuvo ningun resultado CFC.');

%% ── V.3  GUARDAR EXCELS (todo lo necesario para rehacer figuras) ─────────
%
%  cfc_MI.xlsx                       MI/z/p/fase preferida por sesion x capa x condicion x banda
%  cfc_MI_groupsummary.xlsx          media y SD por grupo
%  cfc_phaseamp_persession.xlsx      curvas fase-amplitud por sesion (formato largo)
%  cfc_phaseamp_groupmean.xlsx       curvas medias +- SEM por grupo; 1 hoja por condicion_banda,
%                                    2 ciclos (-180..540) y columnas <grupo>_<capa>_mean/_sem/_n
%  cfc_comodulogram_persession.xlsx  comodulogramas por sesion (formato largo)
%  cfc_comodulogram_groupmean.xlsx   matriz media por grupo (1 hoja por condicion_grupo_capa;
%                                    fila 1 = freq. fase, columna A = freq. amplitud)

cfcTable = cell2table(cfcRows, 'VariableNames', ...
    {'ID','ExpCond','Layer','Condition','GammaBand','PhaseRef','Dur_s','MI','MI_z','MI_p','PrefPhase_deg'});

xl = @(name) fullfile(cfcOut, name);
for nm = {'cfc_MI.xlsx','cfc_MI_groupsummary.xlsx','cfc_phaseamp_persession.xlsx', ...
          'cfc_phaseamp_groupmean.xlsx','cfc_comodulogram_persession.xlsx','cfc_comodulogram_groupmean.xlsx'}
    if exist(xl(nm{1}), 'file'), delete(xl(nm{1})); end
end

nL = numel(C.layerNames);  nB = numel(C.gammaNames);  nC = numel(C.condNames);
ctr  = -180 + ((1:C.nBins) - 0.5) * 360 / C.nBins;
grpS = cellfun(@(s) s.expCond, sessList, 'UniformOutput', false);

% 1) MI
writetable(cfcTable, xl('cfc_MI.xlsx'));
writetable(groupsummary(cfcTable, {'ExpCond','Layer','Condition','GammaBand'}, {'mean','std'}, {'MI','MI_z'}), ...
    xl('cfc_MI_groupsummary.xlsx'));

% 2) Fase-amplitud por sesion (largo)
paRows = cell(0, 7);
for i = 1:numel(sessList)
    S = sessList{i};
    for k = 1:nL
        for b = 1:nB
            for c = 1:nC
                v = squeeze(S.ampDist(k,b,c,:));
                if any(isnan(v)), continue; end
                paRows = [paRows; [repmat({S.id, S.expCond, C.layerNames{k}, C.condNames{c}, C.gammaNames{b}}, C.nBins, 1), ...
                    num2cell(ctr(:)), num2cell(v)]];  %#ok<AGROW>
            end
        end
    end
end
writetable(cell2table(paRows, 'VariableNames', {'ID','ExpCond','Layer','Condition','GammaBand','Phase_deg','AmpNorm'}), ...
    xl('cfc_phaseamp_persession.xlsx'));

% 3) Fase-amplitud media +- SEM por grupo (ancho, 2 ciclos)
xx = [ctr, ctr + 360]';
for c = 1:nC
    for b = 1:nB
        Tw = table(xx, 'VariableNames', {'Phase_deg'});
        for gi = 1:nGroupsC
            sel = find(strcmp(grpS, cfcLabels{gi}));
            for k = 1:nL
                M = nan(numel(sel), C.nBins);
                for q = 1:numel(sel)
                    M(q,:) = squeeze(sessList{sel(q)}.ampDist(k,b,c,:))';
                end
                M = M(~any(isnan(M),2), :);
                if isempty(M), continue; end
                m  = mean(M,1)';
                se = std(M,0,1)' / sqrt(size(M,1));
                base = sprintf('%s_%s', cfcLabels{gi}, C.layerNames{k});
                Tw.(matlab.lang.makeValidName([base '_mean'])) = [m; m];
                Tw.(matlab.lang.makeValidName([base '_sem']))  = [se; se];
                Tw.(matlab.lang.makeValidName([base '_n']))    = repmat(size(M,1), 2*C.nBins, 1);
            end
        end
        writetable(Tw, xl('cfc_phaseamp_groupmean.xlsx'), 'Sheet', sprintf('%s_%s', C.condNames{c}, C.gammaNames{b}));
    end
end

% 4) Comodulogramas
if C.doComod
    [PF, AF] = meshgrid(C.comodPhase, C.comodAmp);
    nPx = numel(PF);
    cmRows = cell(0, 7);
    for i = 1:numel(sessList)
        S = sessList{i};
        for k = 1:nL
            for c = 1:nC
                m = S.comod{k,c};
                if isempty(m), continue; end
                cmRows = [cmRows; [repmat({S.id, S.expCond, C.layerNames{k}, C.condNames{c}}, nPx, 1), ...
                    num2cell(PF(:)), num2cell(AF(:)), num2cell(m(:))]];  %#ok<AGROW>
            end
        end
    end
    writetable(cell2table(cmRows, 'VariableNames', {'ID','ExpCond','Layer','Condition','PhaseFreq_Hz','AmpFreq_Hz','MI'}), ...
        xl('cfc_comodulogram_persession.xlsx'));

    for c = 1:nC
        for gi = 1:nGroupsC
            for k = 1:nL
                M = groupMeanComod(sessList, cfcLabels{gi}, k, c);
                if isempty(M), continue; end
                grpShort = cfcLabels{gi};  grpShort = grpShort(1:min(end, 12));
                shName = regexprep(sprintf('%s_%s_%s', C.condNames{c}, grpShort, C.layerNames{k}), '[\[\]\:\*\?\/\\]', '_');
                cellOut = [{[]}, num2cell(C.comodPhase(:)'); num2cell(C.comodAmp(:)), num2cell(M)];
                writecell(cellOut, xl('cfc_comodulogram_groupmean.xlsx'), 'Sheet', shName);
            end
        end
    end
end

save(xl('cfc_results.mat'), 'cfcTable', 'sessList', 'C', 'cfcLabels', '-v7.3');
fprintf('\n[Parte V] Excels y .mat guardados en: %s\n', cfcOut);

%% ── V.4  FIGURAS ─────────────────────────────────────────────────────────

plotMIBoxes(cfcTable, C, cfcLabels, cfcFigFolder);
plotPhaseAmp(sessList, C, cfcLabels, cfcFigFolder);
if C.doComod, plotComod(sessList, C, cfcLabels, cfcFigFolder); end
fprintf('[Parte V] Completada. Figuras en: %s\n', cfcFigFolder);


%% ========================================================================
%  FUNCIONES LOCALES
%% ========================================================================

function S = cfcSession(recPath, recName, expCond, C)
    cd(recPath);
    nL = numel(C.layerNames);  nB = numel(C.gammaNames);  nC = numel(C.condNames);
    S.id = recName;  S.expCond = expCond;
    S.rows    = cell(0, 11);
    S.ampDist = nan(nL, nB, nC, C.nBins);
    S.comod   = cell(nL, nC);

    % --- Epocas (generadas en Parte I) ---
    tf = load(fullfile(recPath, [recName '.thetaEpochs6to12.states.mat']), 'thisthetaEpochs');
    te = tf.thisthetaEpochs;
    runInt = te.thetaRun.intervals;
    runInt = runInt((runInt(:,2) - runInt(:,1)) >= C.minRunDur, :);
    ints   = {runInt, te.intervals_nonTheta_nonRipples};

    % --- Canales por capa ---
    hl = load(fullfile(recPath, [recName '.hippocampalLayers.channelinfo.mat']), 'hippocampalLayers');
    hl = hl.hippocampalLayers;
    Ly = hl.layers{1,1};
    chan = nan(1, nL);
    for k = 1:nL, chan(k) = pickLayerChannel(Ly, C.layerFields{k}); end
    slmK = find(strcmp(C.layerNames, 'slm'));
    if isnan(chan(slmK))
        chan(slmK) = double(hl.channels(end));
        warning('SLM no encontrado en %s. Usando ultimo canal (%d).', recName, chan(slmK));
    end

    % --- LFP ---
    x = cell(1, nL);  fs = [];  t0 = 0;
    for k = 1:nL
        if isnan(chan(k)), continue; end
        lfp  = getLFP(chan(k), 'noPrompts', true);
        x{k} = double(lfp.data(:,1));
        fs   = lfp.samplingRate;
        if isfield(lfp, 'timestamps') && ~isempty(lfp.timestamps), t0 = lfp.timestamps(1); end
    end
    assert(~isempty(fs), 'No se pudo cargar ningun LFP.');
    maxF = max(max(cellfun(@max, C.gammaBands)), max(C.comodAmp) + C.comodAmpBW/2);
    assert(fs/2 > maxF, 'fs=%g Hz insuficiente para %g Hz.', fs, maxF);
    nS = numel(x{find(~cellfun(@isempty, x), 1)});

    % --- Mascaras por condicion ---
    masks = cell(1, nC);  valid = false(1, nC);  durs = zeros(1, nC);
    for c = 1:nC
        m = false(nS, 1);  iv = ints{c};
        for j = 1:size(iv,1)
            s1 = max(1,  floor((iv(j,1) - t0) * fs) + 1);
            e1 = min(nS, floor((iv(j,2) - t0) * fs) + 1);
            if e1 >= s1, m(s1:e1) = true; end
        end
        masks{c} = m;  durs(c) = sum(m) / fs;  valid(c) = durs(c) >= C.minTotalDur;
    end
    if ~any(valid), error('Duracion total < %g s en todas las condiciones.', C.minTotalDur); end

    % --- Senal de fase por capa ---
    sigIdx = 1:nL;
    if ~strcmp(C.phaseRefLayer, 'self')
        refK = find(strcmp(C.layerNames, C.phaseRefLayer));
        assert(~isempty(refK) && ~isempty(x{refK}), 'Capa de referencia "%s" no disponible.', C.phaseRefLayer);
        sigIdx(:) = refK;
    end

    % --- Bins de fase theta (principal) ---
    binsMain = cell(1, nL);
    for s = unique(sigIdx(~isnan(chan)))
        ph = angle(hilbert(bpf(x{s}, C.phaseBand, fs)));
        binsMain{s} = cell(1, nC);
        for c = 1:nC, binsMain{s}{c} = toBins(ph(masks{c}), C.nBins); end
    end

    % --- MI por capa x banda gamma x condicion ---
    for k = 1:nL
        if isnan(chan(k)), continue; end
        s = sigIdx(k);
        for b = 1:nB
            env = abs(hilbert(bpf(x{k}, C.gammaBands{b}, fs)));
            for c = 1:nC
                if ~valid(c), continue; end
                amp  = env(masks{c});
                bins = binsMain{s}{c};
                [MI, ad] = binMI(bins, amp, C.nBins);
                [mz, mp] = surrStats(bins, amp, C.nBins, MI, C.nSurr, round(C.minShiftS * fs));
                [~, imax] = max(ad);
                prefDeg = -180 + (imax - 0.5) * 360 / C.nBins;
                S.ampDist(k,b,c,:) = ad / sum(ad);
                S.rows(end+1,:) = {recName, expCond, C.layerNames{k}, C.condNames{c}, ...
                    C.gammaNames{b}, C.layerNames{s}, durs(c), MI, mz, mp, prefDeg}; %#ok<AGROW>
            end
        end
    end

    % --- Comodulogramas ---
    if C.doComod
        nPh = numel(C.comodPhase);  nAm = numel(C.comodAmp);
        for k = 1:nL
            if isnan(chan(k)), continue; end
            s = sigIdx(k);
            phB = cell(nC, nPh);
            for ip = 1:nPh
                f  = C.comodPhase(ip);
                ph = angle(hilbert(bpf(x{s}, [f - C.comodPhaseBW/2, f + C.comodPhaseBW/2], fs)));
                for c = 1:nC
                    if valid(c), phB{c,ip} = toBins(ph(masks{c}), C.nBins); end
                end
            end
            com = repmat({nan(nAm, nPh)}, 1, nC);
            for ia = 1:nAm
                fa  = C.comodAmp(ia);
                env = abs(hilbert(bpf(x{k}, [fa - C.comodAmpBW/2, fa + C.comodAmpBW/2], fs)));
                for c = 1:nC
                    if ~valid(c), continue; end
                    amp = env(masks{c});
                    for ip = 1:nPh
                        com{c}(ia, ip) = binMI(phB{c,ip}, amp, C.nBins);
                    end
                end
            end
            for c = 1:nC
                if valid(c), S.comod{k,c} = com{c}; end
            end
        end
    end
end

function S = cfcSessionSafe(recPath, recName, expCond, C, jobIdx)
    rng(C.seed + jobIdx);
    try
        S = cfcSession(recPath, recName, expCond, C);
        S.ok = true;  S.msg = '';
    catch ME
        S = struct('ok', false, 'msg', ME.message, 'id', recName, 'expCond', expCond);
    end
end

function M = groupMeanComod(sessList, grpLabel, k, c)
    stack = [];
    for q = 1:numel(sessList)
        if ~strcmp(sessList{q}.expCond, grpLabel), continue; end
        m = sessList{q}.comod{k,c};
        if ~isempty(m), stack = cat(3, stack, m); end
    end
    if isempty(stack), M = []; else, M = mean(stack, 3, 'omitnan'); end
end

function ch = pickLayerChannel(L, names)
    ch = NaN;
    for k = 1:numel(names)
        if isfield(L, names{k}) && ~isempty(L.(names{k})) && ~isnan(L.(names{k})(1))
            ch = double(L.(names{k})(1));  return;
        end
    end
end

function y = bpf(x, band, fs)
    [z, p, k] = butter(3, band / (fs/2), 'bandpass');
    [sos, g]  = zp2sos(z, p, k);
    y = filtfilt(sos, g, x);
end

function b = toBins(ph, nBins)
    b = uint8(min(floor((ph + pi) / (2*pi) * nBins) + 1, nBins));
end

function [MI, ampDist] = binMI(bins, amp, nBins)
    idx     = double(bins);
    ampDist = accumarray(idx, amp, [nBins 1]) ./ max(accumarray(idx, 1, [nBins 1]), 1);
    p  = ampDist / sum(ampDist);
    pp = p(p > 0);
    MI = (log(nBins) + sum(pp .* log(pp))) / log(nBins);
end

function [z, pval] = surrStats(bins, amp, nBins, MI, nSurr, minShift)
    n = numel(amp);
    if n <= 2*minShift + 1, minShift = max(floor(n/4), 1); end
    surr = nan(nSurr, 1);
    for s = 1:nSurr
        sh = randi([minShift, n - minShift]);
        surr(s) = binMI(bins, circshift(amp, sh), nBins);
    end
    z    = (MI - mean(surr)) / std(surr);
    pval = (sum(surr >= MI) + 1) / (nSurr + 1);
end

% ── Figuras ──────────────────────────────────────────────────────────────

function plotMIBoxes(T, C, groupLabels, figFolder)
    nL = numel(C.layerNames);
    for c = 1:numel(C.condNames)
        for b = 1:numel(C.gammaNames)
            f = figure('Color','w','Units','normalized','Position',[0.05 0.2 0.9 0.45],'Visible','off');
            for k = 1:nL
                ax = subplot(1, nL, k);  hold(ax, 'on');
                idx = strcmp(T.Layer, C.layerNames{k}) & strcmp(T.Condition, C.condNames{c}) & ...
                      strcmp(T.GammaBand, C.gammaNames{b}) & ~isnan(T.MI);
                title(ax, C.layerNames{k});
                if ~any(idx), continue; end
                boxplot(T.MI(idx), T.ExpCond(idx), 'GroupOrder', groupLabels, 'Parent', ax);
                for gi = 1:numel(groupLabels)
                    v = T.MI(idx & strcmp(T.ExpCond, groupLabels{gi}));
                    scatter(ax, gi + 0.15*(rand(size(v)) - 0.5), v, 18, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
                end
                if k == 1, ylabel(ax, 'MI (Tort)'); end
                set(ax, 'Box', 'off', 'FontSize', 10);
            end
            sgtitle(sprintf('Theta–%s MI — %s', C.gammaNames{b}, C.condNames{c}));
            saveas(f, fullfile(figFolder, sprintf('MI_%s_%s.png', C.gammaNames{b}, C.condNames{c})));
            close(f);
        end
    end
end

function plotPhaseAmp(sessList, C, groupLabels, figFolder)
    nL = numel(C.layerNames);  nB = numel(C.gammaNames);  nG = numel(groupLabels);
    cols = lines(nG);
    ctr  = -180 + ((1:C.nBins) - 0.5) * 360 / C.nBins;
    xx   = [ctr, ctr + 360];
    grp  = cellfun(@(s) s.expCond, sessList, 'UniformOutput', false);
    for c = 1:numel(C.condNames)
        f = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.85],'Visible','off');
        for b = 1:nB
            for k = 1:nL
                ax = subplot(nB, nL, (b-1)*nL + k);  hold(ax, 'on');
                for gi = 1:nG
                    sel = find(strcmp(grp, groupLabels{gi}));
                    M = nan(numel(sel), C.nBins);
                    for q = 1:numel(sel)
                        M(q,:) = squeeze(sessList{sel(q)}.ampDist(k,b,c,:))';
                    end
                    M = M(~any(isnan(M),2), :);
                    if isempty(M), continue; end
                    m  = mean(M,1);
                    se = std(M,0,1) / sqrt(size(M,1));
                    fill(ax, [xx fliplr(xx)], [[m+se m+se] fliplr([m-se m-se])], cols(gi,:), ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    plot(ax, xx, [m m], 'Color', cols(gi,:), 'LineWidth', 1.5, 'DisplayName', groupLabels{gi});
                end
                yline(ax, 1/C.nBins, 'k--', 'HandleVisibility', 'off');
                xlim(ax, [-180 540]);  set(ax, 'XTick', -180:180:540, 'Box', 'off', 'FontSize', 9);
                if b == 1, title(ax, C.layerNames{k}); end
                if k == 1, ylabel(ax, sprintf('%s\namp. norm.', C.gammaNames{b})); end
                if b == nB, xlabel(ax, 'Fase theta (deg)'); end
                if b == 1 && k == nL, legend(ax, 'Location', 'best'); end
            end
        end
        sgtitle(sprintf('Distribucion fase-amplitud — %s', C.condNames{c}));
        saveas(f, fullfile(figFolder, sprintf('phaseamp_%s.png', C.condNames{c})));
        close(f);
    end
end

function plotComod(sessList, C, groupLabels, figFolder)
    nL = numel(C.layerNames);  nG = numel(groupLabels);
    grp = cellfun(@(s) s.expCond, sessList, 'UniformOutput', false);
    for c = 1:numel(C.condNames)
        M = cell(nG, nL);
        for gi = 1:nG
            sel = find(strcmp(grp, groupLabels{gi}));
            for k = 1:nL
                stack = [];
                for q = sel
                    m = sessList{q}.comod{k,c};
                    if ~isempty(m), stack = cat(3, stack, m); end
                end
                if ~isempty(stack), M{gi,k} = mean(stack, 3, 'omitnan'); end
            end
        end
        nonEmpty = ~cellfun(@isempty, M);
        if ~any(nonEmpty(:)), continue; end
        cmax = max(cellfun(@(m) max(m(:)), M(nonEmpty)));

        f = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.85],'Visible','off');
        for gi = 1:nG
            for k = 1:nL
                ax = subplot(nG, nL, (gi-1)*nL + k);
                if isempty(M{gi,k}), axis(ax, 'off'); continue; end
                imagesc(ax, C.comodPhase, C.comodAmp, M{gi,k});
                axis(ax, 'xy');  clim(ax, [0 cmax]);  colormap(ax, 'jet');
                if gi == 1, title(ax, C.layerNames{k}); end
                if k == 1, ylabel(ax, sprintf('%s\nFreq. amplitud (Hz)', groupLabels{gi}), 'Interpreter', 'none'); end
                if gi == nG, xlabel(ax, 'Freq. fase (Hz)'); end
                colorbar(ax);
            end
        end
        sgtitle(sprintf('Comodulograma medio (MI) — %s', C.condNames{c}));
        saveas(f, fullfile(figFolder, sprintf('comodulogram_%s.png', C.condNames{c})));
        close(f);
    end
end