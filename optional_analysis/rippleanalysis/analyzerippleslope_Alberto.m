%% ========================================================================
%  RIPPLE SLOPE CHECK - control de la capa anatomica del canal de ripples
%  ========================================================================
%
%  OBJETIVO
%    Comprobar si el canal de ripples (pyramidal) cae en la misma posicion
%    anatomica en todos los registros, SIN repetir la deteccion, y si las
%    medidas de ripples son comparables entre grupos tras ajustar por esa
%    posicion (n = registro anidado en animal).
%
%  ENTRADA
%    all_ripples_original.xlsx   -> ID, expcond, sinktype, PeakNormPow...
%    ripplewaveform.mat          -> allRipplesRaw (misma orden de filas)
%    (Parte B opcional) <sesion>.ripples.events.mat, hippocampalLayers y LFP
%
%  ESTRUCTURA
%    R0  parametros (TODOS aqui)
%    R1-R6  Parte A: slope/peakDev por ripple, resumen, tests, figura, tablas
%    Parte B: perfil en profundidad (RUN_DEPTH_PROFILE = true la recalcula;
%             false carga depth_profile.mat guardado)
%    R7  ajuste por posicion: animales -> colinealidad -> Spearman ->
%        sensibilidad -> modelos mixtos
%    Funciones locales al final (obligatorio en scripts de MATLAB)
%
%  SALIDA (carpeta slope_check)
%    ripple_slope_*.xlsx, ripple_slope_check.png, depth_profile_*.{xlsx,png,mat},
%    animal_map.xlsx, ripple_slope_adjusted_models.xlsx
%
%  Requiere: Statistics and Machine Learning Toolbox, Signal Processing
%  Toolbox, MATLAB >= R2020a (boxchart).
% =========================================================================

%clearvars; close all; clc;

%% -- R0. PARAMETROS ------------------------------------------------------

% Parte A
SLOPE_WIN         = [-0.05 0.05];  % s alrededor del pico del ripple
PEAKDEV_WIN       = 0.01;          % s (+-) alrededor del pico para peakDev
LP_CUTOFF         = 40;            % Hz, paso bajo antes del slope ([] = sin filtro)
PEAKNORM_CUTOFF   = 50;            % mismo corte que en la Parte III
MIN_RIP_SESS      = 20;            % minimo de ripples para incluir sesion en tests
T_FALLBACK        = [];            % vector de tiempos (s) si los nombres de columna no son numericos

% Parte B
RUN_DEPTH_PROFILE = false;         % true = recalcular (lento, pide carpetas); false = cargar depth_profile.mat
N_RIP_PROFILE     = 300;           % ripples por sesion para el perfil
WIN_PROF          = 0.3;           % s, semiventana de extraccion
REL               = -8:8;          % canales relativos al elegido

% R7: ajuste por posicion
ANIMAL_COL        = '';            % columna de T con el animal ('' = usar animal_map.xlsx / dialogo)
ANIMAL_MODE       = 'ask';         % 'ask' | 'use' | 'list' | 'excel' | 'auto'
MEASURE_VARS      = {'PeakNormPow'};   % columnas de all_ripples_original.xlsx a testear
MODEL_LEVEL       = 'ripple';      % 'ripple'   -> (1|animal)+(1|animal:reg) sobre ripples
                                   % 'registro' -> mediana por registro, (1|animal)
COVARS            = {'medSlope','pctPosSlope','pctRad'};  % covariables a nivel registro (columnas de S)
RAD_MIN           = 70;            % sensibilidad sin perfil: % minimo de sinks rad
POL_RANGE         = [30 70];       % sensibilidad sin perfil: rango de % slope>0
OFFSET_MAX        = 3;             % sensibilidad con perfil: |offsetMaxRipple| maximo (canales)
ZC_MAX            = 2;             % sensibilidad con perfil: |zcRel| maximo (canales)
MIN_REG_SENS      = 3;             % minimo de registros por grupo para correr la sensibilidad

rng(1);
clear profTab R rows                % evita restos de ejecuciones anteriores

%% -- R1. CARGAR DATOS ----------------------------------------------------

rippleFolder = uigetdir(pwd, 'Carpeta con all_ripples_original.xlsx y ripplewaveform.mat');
if rippleFolder == 0, error('Abortado.'); end
outFolder = fullfile(rippleFolder, 'slope_check');
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

T = readtable(fullfile(rippleFolder, 'all_ripples_original.xlsx'));
T.ID       = string(T.ID);
T.expcond  = string(T.expcond);
T.sinktype = string(T.sinktype);

W = load(fullfile(rippleFolder, 'ripplewaveform.mat'), 'allRipplesRaw');
X = table2array(W.allRipplesRaw);

assert(size(X,1) == height(T), ...
    ['Filas de allRipplesRaw (%d) != filas de all_ripples_original.xlsx (%d). ' ...
     'No se puede alinear ripple a ripple.'], size(X,1), height(T));

% Vector de tiempos desde los nombres de columna
t = str2double(string(W.allRipplesRaw.Properties.VariableNames));
if any(isnan(t))
    if isempty(T_FALLBACK)
        error('Nombres de columna no numericos. Define T_FALLBACK (vector de tiempos en s).');
    end
    t = T_FALLBACK(:)';
    assert(numel(t) == size(X,2), 'T_FALLBACK no coincide con el numero de columnas.');
end
if max(abs(t)) > 5, t = t / 1000; end     % ms -> s
fs = 1 / median(diff(t));
fprintf('Ripples: %d | ventana: [%.3f %.3f] s | fs ~ %.0f Hz\n', size(X,1), t(1), t(end), fs);

%% -- R2. SLOPE Y PEAKDEV POR RIPPLE --------------------------------------

if ~isempty(LP_CUTOFF)
    [bL, aL] = butter(3, LP_CUTOFF/(fs/2), 'low');
    ok  = all(isfinite(X), 2);
    Xlp = nan(size(X));
    Xlp(ok,:) = filtfilt(bL, aL, X(ok,:)')';
else
    Xlp = X;
end

iw = t >= SLOPE_WIN(1) & t <= SLOPE_WIN(2);
assert(nnz(iw) >= 5, 'SLOPE_WIN deja menos de 5 muestras dentro de la ventana guardada.');
tc = t(iw)';  tc = tc - mean(tc);
Xw = Xlp(:, iw);
T.slope   = (Xw - mean(Xw, 2)) * tc / sum(tc.^2);          % unidades LFP / s
T.peakDev = mean(Xlp(:, abs(t) <= PEAKDEV_WIN), 2) - mean(Xlp, 2);

% Mismo filtro que Parte III (PeakNormPow); el filtro MAD no se guarda en disco
keep = T.PeakNormPow < PEAKNORM_CUTOFF & isfinite(T.slope);
Tk   = T(keep, :);
Tk.slopeZ = (Tk.slope - median(Tk.slope)) / (1.4826 * mad(Tk.slope, 1));
fprintf('Ripples tras filtro PeakNormPow: %d de %d\n', height(Tk), height(T));

groups = unique(Tk.expcond, 'stable');
nG     = numel(groups);
ids    = unique(Tk.ID, 'stable');
nS     = numel(ids);

if nG <= 3
    cols = [1.00 0.55 0.75; 0.65 0.85 1.00; 0.05 0.25 0.65];
    cols = cols(1:nG, :);
else
    cols = lines(nG);
end

%% -- R3. RESUMEN POR SESION ----------------------------------------------

sGrp = strings(nS,1);  sN = zeros(nS,1);
sMed = nan(nS,1);  sMAD = nan(nS,1);  sIQR = nan(nS,1);  sPos = nan(nS,1);
sDev = nan(nS,1);  sDevPos = nan(nS,1);
sRad = nan(nS,1);  sSlm = nan(nS,1);  sOth = nan(nS,1);

for i = 1:nS
    m = Tk.ID == ids(i);
    s = Tk.slope(m);
    sGrp(i)    = Tk.expcond(find(m, 1));
    sN(i)      = nnz(m);
    sMed(i)    = median(s);
    sMAD(i)    = 1.4826 * mad(s, 1);
    sIQR(i)    = iqr(s);
    sPos(i)    = 100 * mean(s > 0);
    sDev(i)    = median(Tk.peakDev(m));
    sDevPos(i) = 100 * mean(Tk.peakDev(m) > 0);
    sRad(i)    = 100 * mean(Tk.sinktype(m) == "rad");
    sSlm(i)    = 100 * mean(Tk.sinktype(m) == "slm");
    sOth(i)    = 100 * mean(Tk.sinktype(m) == "other");
end

S = table(ids, sGrp, sN, sMed, sMAD, sIQR, sPos, sDev, sDevPos, sRad, sSlm, sOth, ...
    'VariableNames', ["ID","expcond","nRip","medSlope","madSlope","iqrSlope", ...
    "pctPosSlope","medDev","pctPosDev","pctRad","pctSlm","pctOther"]);

S.polaridad = repmat("mixta", nS, 1);
S.polaridad(S.pctPosSlope >= 80) = "arriba";
S.polaridad(S.pctPosSlope <= 20) = "abajo";

fprintf('\n=== Resumen por sesion ===\n');
disp(S);

%% -- R4. RESUMEN POR GRUPO Y TESTS ---------------------------------------

gN = zeros(nG,1); gSess = zeros(nG,1); gMed = nan(nG,1); gMAD = nan(nG,1);
gPos = nan(nG,1); gUp = zeros(nG,1); gDown = zeros(nG,1); gMix = zeros(nG,1);
gICC = nan(nG,1);

for g = 1:nG
    m  = Tk.expcond == groups(g);
    ms = S.expcond == groups(g);
    gN(g)    = nnz(m);
    gSess(g) = nnz(ms);
    gMed(g)  = median(Tk.slope(m));
    gMAD(g)  = 1.4826 * mad(Tk.slope(m), 1);
    gPos(g)  = 100 * mean(Tk.slope(m) > 0);
    gUp(g)   = sum(S.polaridad(ms) == "arriba");
    gDown(g) = sum(S.polaridad(ms) == "abajo");
    gMix(g)  = sum(S.polaridad(ms) == "mixta");
    gICC(g)  = sessionICC(Tk(m, :));
end

G = table(groups, gSess, gN, gMed, gMAD, gPos, gUp, gDown, gMix, gICC, ...
    'VariableNames', ["expcond","nSesiones","nRipples","medSlope","madSlope", ...
    "pctPosSlope","sesArriba","sesAbajo","sesMixta","ICC_sesion"]);
fprintf('\n=== Resumen por grupo ===\n');
disp(G);

% Slope por grupo x tipo de sink
GS = groupsummary(Tk, {'expcond','sinktype'}, 'median', {'slope','peakDev'});
fprintf('\n=== Slope por grupo x sink ===\n');
disp(GS);

fprintf('\n=== Tests entre grupos ===\n');
Sv = S(S.nRip >= MIN_RIP_SESS, :);
gV = categorical(Sv.expcond);

try
    p = vartestn(Tk.slope, categorical(Tk.expcond), 'TestType', 'BrownForsythe', 'Display', 'off');
    fprintf('  Varianza del slope (ripple-level, Brown-Forsythe): p = %.4f  [pseudoreplicacion: orientativo]\n', p);
catch ME, fprintf('  vartestn ripple-level fallo: %s\n', ME.message); end

try
    p = vartestn(Sv.medSlope, gV, 'TestType', 'BrownForsythe', 'Display', 'off');
    fprintf('  Varianza de medianas de slope entre sesiones (Brown-Forsythe): p = %.4f\n', p);
catch ME, fprintf('  vartestn session-level fallo: %s\n', ME.message); end

testVars = {'medSlope','madSlope','pctPosSlope','pctPosDev','pctRad','pctSlm','pctOther'};
for q = 1:numel(testVars)
    try
        p = kruskalwallis(Sv.(testVars{q}), gV, 'off');
        fprintf('  %-12s Kruskal-Wallis (sesiones): p = %.4f\n', testVars{q}, p);
    catch ME, fprintf('  %s fallo: %s\n', testVars{q}, ME.message); end
end

try
    Tm = Tk; Tm.expcond = categorical(Tm.expcond); Tm.ID = categorical(Tm.ID);
    lme = fitlme(Tm, 'slopeZ ~ expcond + (1|ID)');
    fprintf('\n  Modelo mixto slopeZ ~ expcond + (1|ID):\n');
    disp(anova(lme));
catch ME, fprintf('  fitlme fallo: %s\n', ME.message); end

%% -- R5. FIGURA ----------------------------------------------------------

[~, gIdxS] = ismember(S.expcond, groups);
[~, ord]   = sortrows([gIdxS, S.medSlope]);
lo = prctile(Tk.slope, 1);  hi = prctile(Tk.slope, 99);

fig = figure('Color','w','Units','normalized','Position',[0.03 0.05 0.94 0.85]);
tiledlayout(2, 3, 'TileSpacing', 'compact');

% (a) histograma de slope por grupo
nexttile; hold on
edges = linspace(lo, hi, 61);
for g = 1:nG
    histogram(Tk.slope(Tk.expcond == groups(g)), edges, 'Normalization', 'probability', ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceColor', cols(g,:));
end
xline(0, 'k--'); legend(groups, 'Box', 'off');
xlabel('slope (unid. LFP/s)'); ylabel('Probability'); title('Slope por ripple'); box off

% (b) boxchart por sesion
nexttile; hold on
xcat = categorical(Tk.ID, ids(ord));
for g = 1:nG
    mk = Tk.expcond == groups(g);
    boxchart(xcat(mk), Tk.slope(mk), 'BoxFaceColor', cols(g,:), 'MarkerStyle', 'none');
end
yline(0, 'k--'); ylim([lo hi]); xtickangle(90);
ylabel('slope'); title('Slope por sesion (ordenado por grupo)'); box off

% (c) % slopes positivos por sesion
nexttile; hold on
b = bar(1:nS, S.pctPosSlope(ord), 'FaceColor', 'flat');
b.CData = cols(gIdxS(ord), :);
yline(50, 'k--'); xticks(1:nS); xticklabels(S.ID(ord)); xtickangle(90);
ylabel('% ripples con slope > 0'); title('Polaridad por sesion'); ylim([0 100]); box off

% (d) composicion de sinks por grupo
nexttile;
sinkMat = nan(nG, 3);
for g = 1:nG
    m = Tk.expcond == groups(g);
    sinkMat(g,:) = 100 * [mean(Tk.sinktype(m)=="rad"), mean(Tk.sinktype(m)=="slm"), mean(Tk.sinktype(m)=="other")];
end
bar(sinkMat, 'stacked'); xticklabels(groups); legend({'rad','slm','other'}, 'Box', 'off');
ylabel('% ripples'); title('Sink CSD por grupo'); box off

% (e) mediana de slope por sesion y grupo
nexttile; hold on
for g = 1:nG
    ms = S.expcond == groups(g);
    xj = g + 0.2*(rand(nnz(ms),1) - 0.5);
    scatter(xj, S.medSlope(ms), 50, cols(g,:), 'filled', 'MarkerEdgeColor', 'k');
    plot(g + [-0.25 0.25], median(S.medSlope(ms))*[1 1], 'k-', 'LineWidth', 2);
end
yline(0, 'k--'); xlim([0.5 nG+0.5]); xticks(1:nG); xticklabels(groups);
ylabel('mediana slope por sesion'); title('Sesiones'); box off

% (f) histograma peakDev
nexttile; hold on
lo2 = prctile(Tk.peakDev, 1); hi2 = prctile(Tk.peakDev, 99);
edges2 = linspace(lo2, hi2, 61);
for g = 1:nG
    histogram(Tk.peakDev(Tk.expcond == groups(g)), edges2, 'Normalization', 'probability', ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceColor', cols(g,:));
end
xline(0, 'k--'); xlabel('peakDev (unid. LFP)'); title('Valor en el pico - media ventana'); box off

saveas(fig, fullfile(outFolder, 'ripple_slope_check.png'));

%% -- R6. GUARDAR TABLAS --------------------------------------------------

writetable(Tk(:, {'ID','expcond','sinktype','slope','peakDev'}), fullfile(outFolder, 'ripple_slope_per_ripple.xlsx'));
writetable(S,  fullfile(outFolder, 'ripple_slope_sessions.xlsx'));
writetable(G,  fullfile(outFolder, 'ripple_slope_groups.xlsx'));
writetable(GS, fullfile(outFolder, 'ripple_slope_group_by_sink.xlsx'));
fprintf('\n[Parte A] Completada. Resultados en: %s\n', outFolder);

%% ========================================================================
%  PARTE B - PERFIL EN PROFUNDIDAD
%% ========================================================================

profFile = fullfile(outFolder, 'depth_profile.mat');

if RUN_DEPTH_PROFILE

    startDir = pwd;
    groupFolders = strings(nG, 1);
    for g = 1:nG
        fprintf('>>> Selecciona la carpeta del grupo "%s":\n', groups(g));
        f = uigetdir(rippleFolder, sprintf('Carpeta del grupo: %s', groups(g)));
        if f == 0, error('Abortado.'); end
        groupFolders(g) = string(f);
    end

    lpc = LP_CUTOFF; if isempty(lpc), lpc = 40; end
    profRip = nan(0, numel(REL));
    profSW  = nan(0, numel(REL));
    profGrp = zeros(0, 1);
    rows    = {};

    for g = 1:nG
        sess = S.ID(S.expcond == groups(g));
        for k = 1:numel(sess)
            recPath = fullfile(groupFolders(g), sess(k));
            fprintf('  [%s] %s\n', groups(g), sess(k));
            try
                cd(recPath);

                L  = load(fullfile(recPath, sess(k) + ".hippocampalLayers.channelinfo.mat"), 'hippocampalLayers');
                hl = L.hippocampalLayers;
                chans = hl.channels(:)';
                pyrCh = hl.layers{1,1}.pyramidal;

                E  = load(fullfile(recPath, sess(k) + ".ripples.events.mat"));
                fn = fieldnames(E);
                iR = find(ismember(fn, {'ripples','ripplesEvents'}), 1);
                if isempty(iR)
                    iR = find(cellfun(@(f) isstruct(E.(f)) && (isfield(E.(f),'peaks') || isfield(E.(f),'timestamps')), fn), 1);
                end
                if isempty(iR), error('Sin struct de ripples en el .mat (variables: %s).', strjoin(fn, ', ')); end
                rEv = E.(fn{iR});
                if isfield(rEv, 'peaks') && ~isempty(rEv.peaks)
                    pk = rEv.peaks(:);
                else
                    pk = mean(rEv.timestamps, 2);
                end

                lfp  = getLFP('all', 'noPrompts', true);
                fsL  = lfp.samplingRate;
                if isfield(lfp, 'channels') && ~isempty(lfp.channels)
                    lfpCh = lfp.channels(:)';
                else
                    lfpCh = 0:size(lfp.data, 2)-1;
                end
                [okc, cidx] = ismember(chans, lfpCh);
                chans = chans(okc);  cidx = cidx(okc);
                ip = find(chans == pyrCh, 1);
                if isempty(ip), error('Canal pyramidal %d no esta en hippocampalLayers.channels.', pyrCh); end

                ts   = lfp.timestamps(:);
                sIdx = round(interp1(ts, 1:numel(ts), pk, 'nearest'));
                hw   = round(WIN_PROF * fsL);
                sIdx = sIdx(isfinite(sIdx) & sIdx > hw & sIdx <= numel(ts) - hw);
                if numel(sIdx) > N_RIP_PROFILE
                    sIdx = sIdx(randperm(numel(sIdx), N_RIP_PROFILE));
                end
                idx = sIdx(:) + (-hw:hw);

                [bR, aR] = butter(3, [80 240]/(fsL/2), 'bandpass');
                [bP, aP] = butter(3, lpc/(fsL/2), 'low');
                tt  = (-hw:hw) / fsL;
                cen = abs(tt) <= 0.02;
                cpk = abs(tt) <= PEAKDEV_WIN;

                nC = numel(chans);
                rr = nan(1, nC);  sw = nan(1, nC);
                for c = 1:nC
                    col = double(lfp.data(:, cidx(c)));
                    seg = col(idx);
                    rb  = filtfilt(bR, aR, seg')';
                    lp  = filtfilt(bP, aP, seg')';
                    rr(c) = mean(sqrt(mean(rb(:, cen).^2, 2)));
                    sw(c) = mean(mean(lp(:, cpk), 2) - mean(lp, 2));
                end

                [~, imax] = max(rr);
                zc = find(sign(sw(1:end-1)) .* sign(sw(2:end)) < 0);
                if isempty(zc)
                    zcRel = NaN;
                else
                    [~, mm] = min(abs((zc + 0.5) - ip));
                    zcRel = zc(mm) + 0.5 - ip;
                end

                relRip = nan(1, numel(REL));  relSW = nan(1, numel(REL));
                for q = 1:numel(REL)
                    j = ip + REL(q);
                    if j >= 1 && j <= nC
                        relRip(q) = rr(j) / max(rr);
                        relSW(q)  = sw(j) / max(abs(sw));
                    end
                end
                profRip(end+1, :) = relRip;      %#ok<SAGROW>
                profSW(end+1, :)  = relSW;       %#ok<SAGROW>
                profGrp(end+1, 1) = g;           %#ok<SAGROW>

                rows{end+1} = table(sess(k), groups(g), nC, ip, imax - ip, sign(sw(ip)), zcRel, ...
                    'VariableNames', ["ID","expcond","nCanales","idxElegido","offsetMaxRipple","signoSW_elegido","zcRel"]); %#ok<SAGROW>
            catch ME
                warning('Sesion %s fallo en Parte B: %s', sess(k), ME.message);
            end
            close all; drawnow;
        end
    end
    cd(startDir);

    if isempty(rows), error('Parte B: ninguna sesion procesada; revisa los warnings.'); end
    profTab = vertcat(rows{:});
    fprintf('\n=== Perfil en profundidad por sesion ===\n');
    disp(profTab);
    fprintf('offsetMaxRipple = (canal con max. potencia ripple) - (canal elegido), en indices de canal.\n');
    fprintf('zcRel = posicion de la inversion de polaridad de la sharp wave relativa al canal elegido.\n');

    try
        gP = categorical(profTab.expcond);
        p1 = kruskalwallis(abs(profTab.offsetMaxRipple), gP, 'off');
        p2 = kruskalwallis(profTab.zcRel, gP, 'off');
        fprintf('  |offsetMaxRipple| entre grupos: p = %.4f\n', p1);
        fprintf('  zcRel entre grupos:            p = %.4f\n', p2);
    catch ME, fprintf('  tests perfil fallo: %s\n', ME.message); end

    fig2 = figure('Color','w','Units','normalized','Position',[0.1 0.2 0.8 0.5]);
    for pn = 1:2
        subplot(1, 2, pn); hold on
        if pn == 1, Mp = profRip; yl = 'ripple RMS / max'; else, Mp = profSW; yl = 'polaridad SW / max|SW|'; end
        for g = 1:nG
            m  = profGrp == g;
            mu = mean(Mp(m,:), 1, 'omitnan');
            se = std(Mp(m,:), 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Mp(m,:)), 1));
            errorbar(REL, mu, se, '-o', 'Color', cols(g,:), 'LineWidth', 1.5, 'MarkerFaceColor', cols(g,:));
        end
        xline(0, 'k:'); if pn == 2, yline(0, 'k--'); end
        xlabel('canal relativo al elegido (pyramidal)'); ylabel(yl); box off
        if pn == 1, legend(groups, 'Box', 'off'); end
    end
    saveas(fig2, fullfile(outFolder, 'depth_profile_check.png'));

    writetable(profTab, fullfile(outFolder, 'depth_profile_sessions.xlsx'));
    save(profFile, 'profRip', 'profSW', 'profGrp', 'profTab', 'REL', '-mat');
    fprintf('\n[Parte B] Completada. Resultados en: %s\n', outFolder);

elseif exist(profFile, 'file') == 2
    load(profFile, 'profTab');
    fprintf('\n[Parte B] Perfil cargado de %s (%d sesiones)\n', profFile, height(profTab));
else
    fprintf('\n[Parte B] Sin perfil (RUN_DEPTH_PROFILE = false y no existe %s). R7 usara la sensibilidad sin perfil.\n', profFile);
end

%% ========================================================================
%  R7 - AJUSTE POR POSICION DEL CANAL (n = registro anidado en animal)
%% ========================================================================

fprintf('\n=== R7. Ajuste por posicion del canal ===\n');

% --- 7.1 Animal por registro ---------------------------------------------
if ~isempty(ANIMAL_COL) && ismember(ANIMAL_COL, Tk.Properties.VariableNames)
    Tk.animal = string(Tk.(ANIMAL_COL));
    [~, ia]   = ismember(S.ID, Tk.ID);
    S.animal  = Tk.animal(ia);
else
    S.animal  = assignAnimals(S.ID, S.expcond, fullfile(outFolder, 'animal_map.xlsx'), ANIMAL_MODE);
    [~, ix]   = ismember(Tk.ID, S.ID);
    Tk.animal = S.animal(ix);
end
fprintf('Animales: %d | registros: %d\n', numel(unique(S.animal)), nS);
for g = 1:nG
    fprintf('  %-14s %2d animales, %2d registros\n', groups(g), ...
        numel(unique(S.animal(S.expcond == groups(g)))), nnz(S.expcond == groups(g)));
end
disp(sortrows(table(S.animal, S.expcond, S.ID, 'VariableNames', {'animal','expcond','ID'})));

% --- 7.2 Covariables de perfil (si hay profTab) ----------------------------
if exist('profTab', 'var') && ~isempty(profTab)
    [tf, ix] = ismember(S.ID, profTab.ID);
    S.absZc     = nan(nS,1);  S.absZc(tf)     = abs(profTab.zcRel(ix(tf)));
    S.absOffset = nan(nS,1);  S.absOffset(tf) = abs(profTab.offsetMaxRipple(ix(tf)));
    COVARS = unique([COVARS, {'absZc','absOffset'}], 'stable');
    fprintf('Covariables de perfil anadidas: absZc, absOffset (%d registros con perfil)\n', nnz(tf));
end

% --- 7.3 Medidas y z-score de covariables ----------------------------------
MEASURE_VARS = MEASURE_VARS(ismember(MEASURE_VARS, Tk.Properties.VariableNames));
if isempty(MEASURE_VARS)
    error('Ninguna MEASURE_VARS existe en la tabla. Revisa los nombres de columna.');
end
for q = 1:numel(MEASURE_VARS)
    mv  = MEASURE_VARS{q};
    tmp = nan(nS, 1);
    for i = 1:nS
        tmp(i) = median(Tk.(mv)(Tk.ID == ids(i)), 'omitnan');
    end
    S.(['med_' mv]) = tmp;
end
for c = 1:numel(COVARS)
    x = S.(COVARS{c});
    S.(['z_' COVARS{c}]) = (x - mean(x, 'omitnan')) / std(x, 'omitnan');
end

inclBase = S.nRip >= MIN_RIP_SESS;
Sv       = S(inclBase, :);
Sv.grp   = categorical(Sv.expcond, groups);

% --- 7.4 Colinealidad grupo-covariable --------------------------------------
fprintf('\n--- Colinealidad grupo-covariable (R2 de covariable ~ grupo) ---\n');
for c = 1:numel(COVARS)
    tb = table(Sv.grp, Sv.(COVARS{c}), 'VariableNames', {'grp','x'});
    lm = fitlm(tb, 'x ~ grp');
    fprintf('  %-12s R2 = %.2f  (VIF ~ %.1f)\n', COVARS{c}, lm.Rsquared.Ordinary, 1/(1 - lm.Rsquared.Ordinary));
    for g = 1:nG
        xg = Sv.(COVARS{c})(Sv.expcond == groups(g));
        fprintf('      %-14s n=%2d  rango [%.3g, %.3g]\n', groups(g), numel(xg), min(xg), max(xg));
    end
end
fprintf('  R2 alto o rangos que no se solapan -> la covariable no separa grupo de posicion: declarar como limitacion.\n');

% --- 7.5 Spearman medida vs covariables -------------------------------------
fprintf('\n--- Spearman medida (mediana/registro) vs covariables ---\n');
for q = 1:numel(MEASURE_VARS)
    y = Sv.(['med_' MEASURE_VARS{q}]);
    for c = 1:numel(COVARS)
        x = Sv.(COVARS{c});
        [rho, pr] = corr(y, x, 'Type', 'Spearman', 'Rows', 'complete');
        fprintf('  %-12s vs %-12s global: rho = %+.2f (p = %.4f)\n', MEASURE_VARS{q}, COVARS{c}, rho, pr);
        for g = 1:nG
            mg = Sv.expcond == groups(g);
            if nnz(mg & isfinite(x) & isfinite(y)) >= 5
                [rho, pr] = corr(y(mg), x(mg), 'Type', 'Spearman', 'Rows', 'complete');
                fprintf('      %-14s rho = %+.2f (p = %.4f, n=%d)\n', groups(g), rho, pr, nnz(mg));
            end
        end
    end
end

% --- 7.6 Sensibilidad: registros con posicion "limpia" -----------------------
usaPerfil = ismember('absOffset', S.Properties.VariableNames);
if usaPerfil
    inclSens = inclBase & S.absOffset <= OFFSET_MAX & S.absZc <= ZC_MAX;   % NaN -> excluido
    critS = sprintf('|offsetMaxRipple| <= %g, |zcRel| <= %g', OFFSET_MAX, ZC_MAX);
else
    inclSens = inclBase & S.pctRad >= RAD_MIN & S.pctPosSlope >= POL_RANGE(1) & S.pctPosSlope <= POL_RANGE(2);
    critS = sprintf('pctRad >= %g, pctPosSlope en [%g %g]', RAD_MIN, POL_RANGE);
end
nSens = arrayfun(@(g) nnz(inclSens & S.expcond == g), groups);
fprintf('\n--- Sensibilidad: %s ---\n', critS);
for g = 1:nG, fprintf('  %-14s %d de %d registros\n', groups(g), nSens(g), gSess(g)); end
runSens = all(nSens >= MIN_REG_SENS);
if ~runSens
    warning('Algun grupo queda con < %d registros: sensibilidad no ejecutada (limitacion a declarar).', MIN_REG_SENS);
end

% --- 7.7 Modelos: grp / grp+cov / grp*cov ------------------------------------
if strcmpi(MODEL_LEVEL, 'ripple')
    RE = '(1|animal) + (1|animal:reg)';
else
    RE = '(1|animal)';
end
fprintf('\n--- Modelos (%s-level): y ~ grp [+ cov | * cov] + %s ---\n', MODEL_LEVEL, RE);

rows = {};
for q = 1:numel(MEASURE_VARS)
    mv = MEASURE_VARS{q};
    zc = cellfun(@(x) ['z_' x], COVARS, 'UniformOutput', false);

    D0 = buildD(Tk, S, groups, mv, zc, MODEL_LEVEL, inclBase);
    r0 = fitAdj(D0, '', RE);
    fprintf('\n  [%s] coeficientes sin ajuste (referencia = %s):\n', mv, groups(1));
    disp(r0.coef0);
    fprintf('  mediana por registro de %s por grupo (todos | sensibilidad):\n', mv);
    yy = S.(['med_' mv]);
    for g = 1:nG
        mg = S.expcond == groups(g);
        fprintf('    %-14s %.4g | %.4g\n', groups(g), median(yy(mg & inclBase), 'omitnan'), median(yy(mg & inclSens), 'omitnan'));
    end

    pSens = NaN;
    if runSens
        Ds    = buildD(Tk, S, groups, mv, zc, MODEL_LEVEL, inclSens);
        rs    = fitAdj(Ds, '', RE);
        pSens = rs.pGrp0;
        fprintf('\n  [%s] coeficientes en sensibilidad:\n', mv);
        disp(rs.coef0);
    end

    for c = 1:numel(COVARS)
        rc = fitAdj(D0, zc{c}, RE);
        rows{end+1} = table(string(mv), string(COVARS{c}), r0.n, r0.pGrp0, rc.pGrp1, ...
            rc.pCov, rc.pInt, pSens, ...
            'VariableNames', ["medida","covariable","n","p_grupo_sin_ajuste", ...
            "p_grupo_ajustado","p_covariable","p_interaccion","p_grupo_sensibilidad"]); %#ok<SAGROW>
    end
end

R = vertcat(rows{:});
disp(R);
fprintf(['Lectura: si p_grupo_ajustado y p_grupo_sensibilidad mantienen la conclusion de p_grupo_sin_ajuste,\n' ...
         'el efecto no depende de la posicion del canal. p_interaccion < 0.05 -> el efecto de la posicion difiere por grupo.\n']);
writetable(R, fullfile(outFolder, 'ripple_slope_adjusted_models.xlsx'));
fprintf('\n[R7] Completada. Resultados en: %s\n', outFolder);

%% ========================================================================
%  FUNCIONES LOCALES (siempre al final del script)
%% ========================================================================

function icc = sessionICC(tbl)
% Fraccion de la varianza del slope (z) explicada por la sesion.
    icc = NaN;
    try
        tbl.ID = categorical(tbl.ID);
        lme = fitlme(tbl, 'slopeZ ~ 1 + (1|ID)');
        [psi, mse] = covarianceParameters(lme);
        icc = psi{1} / (psi{1} + mse);
    catch
    end
end

function animal = assignAnimals(ids, grp, mapFile, mode)
% Asigna sesion -> animal. Parte de una propuesta automatica y del mapa guardado.
% mode: 'ask' (pregunta), 'use' (mapa guardado), 'auto' (propuesta), 'list', 'excel'.
    ids = string(ids(:));  grp = string(grp(:));
    animal  = animalFromID(ids);
    haveMap = false;
    if exist(mapFile, 'file') == 2
        try
            Mm = readtable(mapFile);
            [tf, ix] = ismember(ids, string(Mm.ID));
            man = string(Mm.animal);
            animal(tf) = man(ix(tf));
            haveMap = all(tf);
        catch ME
            warning('No pude leer %s (%s).', mapFile, ME.message);
        end
    end

    if strcmpi(mode, 'ask')
        if haveMap
            c = questdlg(sprintf('Hay un mapa de animales guardado (%d animales).', numel(unique(animal))), ...
                'Animales', 'Usar', 'Editar en lista', 'Editar en Excel', 'Usar');
            if isempty(c), c = 'Usar'; end
        else
            c = questdlg('Asignar sesiones a animales (la propuesta automatica es el punto de partida).', ...
                'Animales', 'Editar en lista', 'Editar en Excel', 'Usar propuesta', 'Editar en lista');
            if isempty(c), c = 'Usar propuesta'; end
        end
        switch c
            case 'Usar',            mode = 'use';
            case 'Usar propuesta',  mode = 'auto';
            case 'Editar en lista', mode = 'list';
            otherwise,              mode = 'excel';
        end
    end

    switch lower(mode)
        case 'list'
            try
                animal = editInList(ids, grp, animal);
                writetable(table(ids, grp, animal, 'VariableNames', {'ID','expcond','animal'}), mapFile);
            catch ME
                warning('La lista fallo (%s). Paso a Excel.', ME.message);
                animal = editInExcel(ids, grp, animal, mapFile);
            end
        case 'excel'
            animal = editInExcel(ids, grp, animal, mapFile);
        otherwise
            if ~haveMap && strcmpi(mode, 'use')
                warning('No hay mapa guardado completo: uso la propuesta automatica.');
            end
    end
end

function animal = editInList(ids, grp, animal)
% Bucle: selecciona (Ctrl/Shift) sesiones del MISMO animal -> "Unir" -> nombre -> OK.
%   Nombre vacio = cada sesion seleccionada pasa a ser su propio animal.
%   "Cancel" en la lista = terminar.
    [~, ord] = sort(grp + "|" + ids);
    iv = 1;
    while true
        items = arrayfun(@(j) sprintf('%-28s | %-12s | %s', animal(j), grp(j), ids(j)), ord, 'UniformOutput', false);
        [k, ok] = listdlg('ListString', items, 'SelectionMode', 'multiple', ...
            'Name', sprintf('%d sesiones -> %d animales', numel(ids), numel(unique(animal))), ...
            'PromptString', {'Selecciona (Ctrl/Shift) sesiones del MISMO animal y pulsa Unir.', ...
                             'Cancel = terminar y guardar.'}, ...
            'ListSize', [780 520], 'OKString', 'Unir', 'InitialValue', iv);
        if ~ok || isempty(k), break; end
        sel = ord(k);
        a = inputdlg({'Nombre del animal (vacio = una sesion por animal):'}, 'Animal', 1, ...
                     {suggestName(ids(sel), animal(sel))});
        if isempty(a), continue; end
        a = strtrim(string(a{1}));
        if strlength(a) == 0, animal(sel) = ids(sel); else, animal(sel) = a; end
        iv = k(1);
    end
end

function animal = editInExcel(ids, grp, animal, mapFile)
% Escribe el mapa, lo abre en Excel, espera a que lo edites y lo relee.
    writetable(table(ids, grp, animal, 'VariableNames', {'ID','expcond','animal'}), mapFile);
    try, winopen(mapFile); catch, fprintf('Abre manualmente: %s\n', mapFile); end
    waitfor(msgbox({'Edita la columna "animal" (mismo nombre = mismo animal).', ...
                    'Guarda, CIERRA el archivo y pulsa OK.'}, 'Animales en Excel'));
    Mm = readtable(mapFile);
    [tf, ix] = ismember(ids, string(Mm.ID));
    if ~all(tf), warning('%d sesiones no estan en el Excel; mantengo su propuesta.', nnz(~tf)); end
    new = string(Mm.animal);
    new = new(ix(tf));
    bad = ismissing(new) | strlength(new) == 0;
    tmp = animal(tf);
    tmp(~bad) = new(~bad);
    animal(tf) = tmp;
end

function nm = suggestName(idsSel, animalSel)
% Prefijo comun cortado en limite de token '_'; con 1 sesion, su animal actual.
    idsSel = string(idsSel);
    if numel(idsSel) == 1, nm = char(animalSel); return; end
    s = cellstr(idsSel);
    p = s{1};
    for j = 2:numel(s)
        m = min(numel(p), numel(s{j}));
        d = find(p(1:m) ~= s{j}(1:m), 1);
        if isempty(d), p = p(1:m); else, p = p(1:d-1); end
    end
    p  = regexprep(p, '_+$', '');
    ok = ~isempty(p) && all(cellfun(@(t) numel(t) == numel(p) || t(numel(p)+1) == '_', s));
    if ~ok, p = regexprep(p, '_[^_]*$', ''); end
    if isempty(p), p = char(animalSel(1)); end
    nm = p;
end

function a = animalFromID(id)
% Propuesta automatica de animal desde el nombre de sesion (se corrige en la lista/Excel).
    id = string(id); a = id;
    for k = 1:numel(id)
        s = id(k);
        t = regexp(s, '(E\d{3}|F\d{3}|B\d{3})', 'match', 'once');
        if ismissing(t), t = regexp(s, '^Chrnaviejo', 'match', 'once'); end
        if ismissing(t), t = regexp(s, '^.*?_\d+(?=_|$)', 'match', 'once'); end
        if ~ismissing(t), a(k) = t; end
    end
end

function D = buildD(Tk, S, groups, mv, zcols, level, incl)
% Tabla para modelos. incl = logico (nS x 1) de registros incluidos.
    if strcmpi(level, 'ripple')
        [~, ix] = ismember(Tk.ID, S.ID);
        keepR   = ismember(Tk.ID, S.ID(incl));
        D = table(Tk.(mv), categorical(Tk.expcond, groups), categorical(Tk.animal), categorical(Tk.ID), ...
            'VariableNames', {'y','grp','animal','reg'});
        for c = 1:numel(zcols), D.(zcols{c}) = S.(zcols{c})(ix); end
        D = D(keepR, :);
    else
        Sx = S(incl, :);
        D = table(Sx.(['med_' mv]), categorical(Sx.expcond, groups), categorical(Sx.animal), categorical(Sx.ID), ...
            'VariableNames', {'y','grp','animal','reg'});
        for c = 1:numel(zcols), D.(zcols{c}) = Sx.(zcols{c}); end
    end
    D = D(isfinite(D.y), :);
    D.grp    = removecats(D.grp);
    D.animal = removecats(D.animal);
    D.reg    = removecats(D.reg);
end

function r = fitAdj(D, cv, RE)
% cv = '' -> solo y ~ grp ; cv = 'z_x' -> ajustado (aditivo) e interaccion.
    r = struct('n', height(D), 'pGrp0', NaN, 'pGrp1', NaN, 'pCov', NaN, 'pInt', NaN, 'coef0', []);
    try
        l0 = fitlme(D, ['y ~ grp + ' RE]);  a = anovaS(l0);  r.coef0 = l0.Coefficients;
        r.pGrp0 = pTerm(a, 'grp');
    catch ME, warning('Modelo base fallo: %s', ME.message); end
    if isempty(cv), return; end
    D = D(isfinite(D.(cv)), :);
    try
        a = anovaS(fitlme(D, ['y ~ grp + ' cv ' + ' RE]));
        r.pGrp1 = pTerm(a, 'grp');
        r.pCov  = pTerm(a, cv);
    catch ME, warning('Modelo ajustado fallo (%s): %s', cv, ME.message); end
    try
        a = anovaS(fitlme(D, ['y ~ grp * ' cv ' + ' RE]));
        r.pInt = pTerm(a, ['grp:' cv]);
    catch ME, warning('Modelo con interaccion fallo (%s): %s', cv, ME.message); end
end

function a = anovaS(lme)
    try, a = anova(lme, 'DFMethod', 'satterthwaite');
    catch, a = anova(lme); end
end

function p = pTerm(a, term)
    i = strcmp(a.Term, term);
    if any(i), p = a.pValue(find(i, 1)); else, p = NaN; end
end