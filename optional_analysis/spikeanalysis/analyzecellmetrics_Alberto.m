%% ========================================================================
%  spike_analysis_pipeline_v2.m
%  Proyecto hormonas sexuales — Pipeline completo de analisis de spike metrics
%  ========================================================================
%
%  ── FLUJO DE TRABAJO ────────────────────────────────────────────────────
%    S0. PARAMETROS     — editar solo este bloque
%    S1. RECOMPUTO      — recalcula event_metrics/state_metrics por sesion
%    S2. CARGA          — carga cell_metrics, spikes, hippocampalLayers
%    S3. LIMPIEZA       — elimina ruido/indeterminados, detecta outliers
%    S4. CONTEOS        — analisis de celulas y ratios INH/EXC
%    S5. LME EN BATCH   — todos los subgrupos x todas las medidas
%
%  DEPENDENCIAS:
%    ProcessCellMetrics          (CellExplorer)
%    flagOutliersPerMeasure.m
%    analyzeCellCounts.m
%    analyzeCellTypeRatio.m
%    analyzeSubtypeRatio.m
%    analyze_ripple_measure_lme_figures_generalized.m
%
%  INSTRUCCIONES:
%    1. Ajusta el bloque S0 si cambias medidas, thresholds o sujetos
%    2. F5 = corre todo | Ctrl+Enter seccion a seccion
%    3. Los dialogos de carpetas aparecen al inicio (S1.1)
%
%  SALIDA (carpeta output elegida al inicio):
%    all_cellmetrics_raw.xlsx          — tabla combinada antes de QC
%    cellsdata_gt_clean.mat            — tabla limpia clasificacion GT
%    cellsdata_put_clean.mat           — tabla limpia clasificacion putativa
%    gt/  y  put/                      — figuras y reportes por sistema
%
% =========================================================================

clearvars; close all; clc;

%% ========================================================================
%  S0. PARAMETROS — editar solo este bloque
%% ========================================================================

%% ── S0.1  Identificadores de sujetos por grupo ───────────────────────────
%  Se extraen automaticamente del nombre de sesion: todo lo que precede
%  al primer bloque de digitos de fecha (6 digitos consecutivos) o a
%  sufijos numericos finales (_1, _2, _one, _two, ...).
%  No es necesario editar este bloque.

SUBJECTS = struct();   % vacio — la asignacion es dinamica en S3
%% ── S0.2  Etiquetas y colores de grupos experimentales ───────────────────

EXP_CONDS    = ["malewt", "maleaut", "femalewt", "femaleaut"];
GROUP_LABELS = {'Male WT', 'Male AUT', 'Female WT', 'Female AUT'};
GROUP_COLORS = [
    132,   197,  247;   % malewt    — azul oscuro
    0, 0, 128;  % maleaut   — azul claro
    247,  180, 210;  % femalewt  — rosa oscuro
    180, 60, 100   % femaleaut — rosa claro
] / 255;
%% ── S0.3  Lista de medidas ───────────────────────────────────────────────
%  Lista fija: 14 medidas base + 3 de inmovilidad + significancia de ripple.
%  NO hay autodeteccion — solo se extraen y analizan estas medidas.
%
%  Agrupadas por tematica:
%    Intrinsecas generales:  ab_ratio, burstIndex_Mizuseki2012, burstIndex_Royer2012,
%                            cv2, firingRate, firingRateCV, thetaModulationIndex,
%                            troughToPeak, population_modIndex
%    Theta / locomocion:     burstIndex_WAKEtheta, cv2_WAKEtheta, firingRate_WAKEtheta
%    Ripple-relacionadas:    ripples_modulationIndex, ripples_modulationRatio,
%                            ripples_modulationSignificanceLevel
%    Inmovilidad (sin SWR):  burstIndex_QWake_noRipples_ThDt,
%                            cv2_QWake_noRipples_ThDt,
%                            firingRate_QWake_noRipples_ThDt

CORE_MEASURES = [ ...
    "ab_ratio", ...
    "burstIndex_Mizuseki2012", ...
    "burstIndex_Royer2012", ...
    "cv2", ...
    "firingRate", ...
    "firingRateCV", ...
    "thetaModulationIndex", ...
    "troughToPeak", ...
    "population_modIndex", ...
    "burstIndex_WAKEtheta", ...
    "cv2_WAKEtheta", ...
    "firingRate_WAKEtheta", ...
    "ripples_modulationIndex", ...
    "ripples_modulationRatio", ...
    "ripples_modulationSignificanceLevel", ...
    "burstIndex_QWake_noRipples_ThDt", ...
    "cv2_QWake_noRipples_ThDt", ...
    "firingRate_QWake_noRipples_ThDt" ...
];

%% ── S0.4  Outlier removal ────────────────────────────────────────────────

% Variables con log-transform (estrictamente positivas, sesgadas a la derecha)
LOG_VARS       = {'firingRate','firingRate_WAKEtheta','firingRate_QWake_noRipples_ThDt'};

% Umbral MAD (numero de MAD escalados para marcar como outlier)
OUTLIER_THRESH = 4.0;

%% ── S0.5  Subtipos inhibitorios para analisis de ratio ───────────────────

INH_SUBTYPES = {'PV+','SST+','ID2+','all_inh'};

%% ── S0.6  Patrones de exclusion (estados de sueno — estudio wake-only) ───

EXCLUDE_PATTERNS = {'NREM','REM','Sleep','sleep','nrem','rem'};

%% ── S0.7  Recomputar ripple metrics ─────────────────────────────────────
%  true  → primer run tras cambiar deteccion de ripples
%  false → cuando cell_metrics ya estan actualizados

RECOMPUTE_RIPPLES = true;

%% ── S0.8  Capas hipocampales ─────────────────────────────────────────────

LAYER_NAMES = ["oriens","pyr","rad","slm"];

%% ========================================================================
%  S1. CONFIGURACION INTERACTIVA: GRUPOS, CARPETAS Y OUTPUT
%% ========================================================================

fprintf('\n=== PASO 1: Configuracion de grupos y carpetas ===\n');

%% ── S1.1  Numero de grupos ───────────────────────────────────────────────

nGroupsStr = inputdlg('Cuantos grupos experimentales tienes?','Grupos',1,{'3'});
if isempty(nGroupsStr), error('Abortado por el usuario.'); end
nGroups = str2double(nGroupsStr{1});
assert(~isnan(nGroups) && nGroups >= 1,'Numero de grupos no valido.');

groupFolders = cell(nGroups,1);
groupNames   = cell(nGroups,1);

for g = 1:nGroups
    lbl = inputdlg(sprintf('Etiqueta del grupo %d (ej: SS, OS, OD...):', g), ...
        sprintf('Grupo %d',g), 1, {char(EXP_CONDS(min(g,numel(EXP_CONDS))))});
    if isempty(lbl), error('Abortado.'); end
    groupNames{g} = strtrim(lbl{1});

    fprintf('>>> Selecciona la carpeta del grupo "%s":\n', groupNames{g});
    folder = uigetdir(pwd, sprintf('Carpeta del grupo: %s', groupNames{g}));
    if folder == 0, error('No se selecciono carpeta. Abortado.'); end
    groupFolders{g} = folder;
end

%% ── S1.2  Carpeta de salida ──────────────────────────────────────────────

fprintf('>>> Selecciona la carpeta OUTPUT para los resultados:\n');
OUTPUT_DIR = uigetdir(groupFolders{1},'OUTPUT — Spike Analysis Pipeline');
if OUTPUT_DIR == 0, error('Abortado.'); end
if ~exist(OUTPUT_DIR,'dir'), mkdir(OUTPUT_DIR); end

fprintf('  Grupos: %s\n', strjoin(groupNames,', '));
fprintf('  Output: %s\n', OUTPUT_DIR);

%% ── S1.3  Listar sesiones por grupo ─────────────────────────────────────

sessionLists = cell(nGroups,1);
for g = 1:nGroups
    d = dir(groupFolders{g});
    d = d([d.isdir] & ~ismember({d.name},{'.','..'}));
    sessionLists{g} = {d.name};
    fprintf('  Grupo "%s": %d sesion(es) encontradas.\n', groupNames{g}, numel(sessionLists{g}));
end

%% ========================================================================
%  S2. RECOMPUTO DE RIPPLE METRICS (opcional)
%  Ejecuta ProcessCellMetrics con event_metrics + state_metrics.
%  Usa solo el archivo de estados theta 6to12; ignora el resto dinamicamente.
%% ========================================================================

if RECOMPUTE_RIPPLES
    fprintf('\n=== PASO 2: Recomputo de ripple metrics ===\n');
    for g = 1:nGroups
        sessions = sessionLists{g};
        fprintf('  Grupo "%s": %d sesiones\n', groupNames{g}, numel(sessions));
        for s = 1:numel(sessions)
            sesPath    = fullfile(groupFolders{g}, sessions{s});
            rippleFile = dir(fullfile(sesPath,'*.ripples.events.mat'));
            if isempty(rippleFile)
                fprintf('    OMITIDO %s — sin ripples.events.mat\n', sessions{s});
                continue
            end
            fprintf('    Procesando %s ... ', sessions{s});
            try
                prevDir   = pwd;
                cd(sesPath);

                % Detectar dinamicamente que archivos de estados existen.
                % Solo se usa el que contiene '6to12' (theta); el resto se ignora.
                allStates   = dir(fullfile(sesPath,'*.states.mat'));
                allStates   = {allStates.name};
                keepPattern = '6to12';
                ignoreList  = allStates(~contains(allStates, keepPattern));
                ignoreTypes = cellfun(@(f) strsplit(f,'.'), ignoreList,'UniformOutput',false);
                ignoreTypes = cellfun(@(p) p{end-2}, ignoreTypes,'UniformOutput',false);
                ignoreTypes = unique(ignoreTypes);

                ProcessCellMetrics( ...
                    'metrics',              {'event_metrics','state_metrics'}, ...
                    'forceReload',          true, ...
                    'ignoreStateTypes',     ignoreTypes, ...
                    'saveMat',              true, ...
                    'showFigures',          false, ...
                    'showGUI',              false, ...
                    'sessionSummaryFigure', false);
                cd(prevDir);
                fprintf('OK\n');
            catch ME
                if exist('prevDir','var'), cd(prevDir); end
                fprintf('FALLIDO: %s\n', ME.message);
            end
        end
    end
else
    fprintf('\n=== PASO 2: Recomputo omitido (RECOMPUTE_RIPPLES = false) ===\n');
end

%% ========================================================================
%  S3. CARGA Y CONSTRUCCION DE TABLA UNIFICADA
%  Por cada sesion: cell_metrics, spikes, hippocampalLayers.
%  Asigna capa hipocampal, tipo celular (GT y putativo), y columna subject.
%% ========================================================================

fprintf('\n=== PASO 3: Carga de cell metrics ===\n');

allData     = [];
varlist     = [];   % se fija tras el primer load exitoso

for g = 1:nGroups
    grpName  = groupNames{g};
    sessions = sessionLists{g};
    nSess    = numel(sessions);

    % Recuperar lista de sujetos del grupo (si esta definida en S0.1)
    if isfield(SUBJECTS, grpName)
        subjIDs = SUBJECTS.(grpName);
    else
        subjIDs = {};
    end

    fprintf('\n==== Procesando grupo "%s": %d sesion(es) ====\n', grpName, nSess);

    for s = 1:numel(sessions)
        sesName = sessions{s};
        sesPath = fullfile(groupFolders{g}, sesName);

        cmFile  = fullfile(sesPath, [sesName '.cell_metrics.cellinfo.mat']);
        spkFile = fullfile(sesPath, [sesName '.spikes.cellinfo.mat']);
        layFile = fullfile(sesPath, [sesName '.hippocampalLayers.channelinfo.mat']);

        fprintf('  [%d/%d] %s\n', s, nSess, sesName);

        if ~exist(cmFile,'file')
            fprintf('    cell_metrics no encontrado. Omitiendo.\n');
            continue
        end

        try
            load(cmFile,  'cell_metrics');
            load(spkFile, 'spikes');
            load(layFile, 'hippocampalLayers');
        catch ME
            fprintf('    ERROR cargando archivos: %s. Omitiendo.\n', ME.message);
            continue
        end

        % ── Fijar lista de variables en el primer load exitoso ────────────
        if isempty(varlist)
            varlist = buildVarList(cell_metrics, CORE_MEASURES, EXCLUDE_PATTERNS);
            fprintf('    Lista de medidas fijada: %d / %d encontradas\n', ...
                numel(varlist), numel(CORE_MEASURES));
        end

        nCells = cell_metrics.general.cellCount;

        % ── Extraer medidas numericas ─────────────────────────────────────
        thismeasures = nan(nCells, numel(varlist));
        for jj = 1:numel(varlist)
            vn = varlist(jj);
            if isfield(cell_metrics, vn)
                vals = cell_metrics.(vn);
                if isnumeric(vals) && numel(vals) == nCells
                    thismeasures(:,jj) = vals(:);
                end
            end
        end

        % ── Tipos celulares: GT y putativo (sistemas separados, nunca mezclados) ──
        %   celltype_gt      : ground_truth_classification (CAMK2_DEEP, CAMK2_SUP, PV+, SST+, ID2+)
        %   celltype_putative: putativeCellType (Pyramidal Cell, Wide Interneuron, Narrow Interneuron)
        thisfileclass_gt  = extractCellTypes(cell_metrics, nCells, 'gt');
        thisfileclass_put = extractCellTypes(cell_metrics, nCells, 'putative');

        % ── Asignar capa hipocampal a cada celula ─────────────────────────
        spikech          = spikes.phy_maxWaveformCh1(:);
        layers           = hippocampalLayers.layers{1,1}.all;
        spikelayers      = zeros(nCells,1);
        for kk = 1:nCells
            [~, spikelayers(kk)] = min(abs(layers - spikech(kk)));
        end
        spikelayersnames = strings(nCells,1);
        spikelayersnames(spikelayers==1) = LAYER_NAMES(1);
        spikelayersnames(spikelayers==2) = LAYER_NAMES(2);
        spikelayersnames(spikelayers==3) = LAYER_NAMES(3);
        spikelayersnames(spikelayers==4) = LAYER_NAMES(4);

        % ── Asignar ID de sujeto ──────────────────────────────────────────
        thisSubject = repmat(string(extractSubjectID(sesName)), nCells, 1);
        for ss = 1:numel(subjIDs)
            if contains(sesName, subjIDs{ss})
                thisSubject(:) = subjIDs{ss};
                break
            end
        end
        % Fallback: si no coincide ningun ID, extraer prefijo del nombre de sesion
        if all(thisSubject == "unknown")
            try
                prefix = extractBefore(string(sesName),'_');
                if ~ismissing(prefix), thisSubject(:) = prefix; end
            catch
                % mantener "unknown"
            end
        end

        % ── Construir tabla de esta sesion ────────────────────────────────
        thisID      = repmat(string(sesName),  nCells, 1);
        thisExpCond = repmat(string(grpName),  nCells, 1);

        thisMeasTab = array2table(thismeasures, 'VariableNames', cellstr(varlist));
        thisMetaTab = table(thisID, thisExpCond, ...
            thisfileclass_gt, thisfileclass_put, ...
            spikech, spikelayersnames, thisSubject, ...
            'VariableNames', ["ID","expcond","celltype_gt", ...
                              "celltype_putative","spikech","layer","subject"]);

        thisRec = [thisMetaTab, thisMeasTab];
        allData = [allData; thisRec]; %#ok<AGROW>

        fprintf('    OK: %d celulas\n', nCells);

    end % sesiones
end % grupos

fprintf('\n  Total unidades cargadas: %d\n', height(allData));

% Guardar tabla cruda
rawFile = fullfile(OUTPUT_DIR,'all_cellmetrics_raw.xlsx');
writetable(allData, rawFile);
fprintf('  Tabla cruda guardada: %s\n', rawFile);

%% ========================================================================
%  S4. LIMPIEZA DE DATOS
%  Dos tablas limpias paralelas, una por sistema de clasificacion.
%    GT:      elimina Noisy_unit, Undetermined
%    Putativo: elimina Unknown, etiquetas vacias
%  La deteccion de outliers se aplica a ambas con la misma lista de medidas.
%% ========================================================================

fprintf('\n=== PASO 4: Limpieza de datos ===\n');

% Lista de medidas validas (excluye estados de sueno y columnas inexistentes)
validMeasures = cellstr(CORE_MEASURES);
validMeasures = validMeasures( ...
    ismember(validMeasures, allData.Properties.VariableNames) & ...
    ~cellfun(@(n) any(contains(n, EXCLUDE_PATTERNS)), validMeasures));

fprintf('  Medidas para analisis (%d):\n', numel(validMeasures));
for mi = 1:numel(validMeasures)
    fprintf('    %2d. %s\n', mi, validMeasures{mi});
end

% ── Tabla clasificacion GT ────────────────────────────────────────────────
% Usa celltype_gt; elimina Noisy_unit y Undetermined.
% Renombra celltype_gt → celltype para compatibilidad con funciones downstream.
removeTypes_gt = {'Noisy_unit','Undetermined'};
keepGT         = ~ismember(string(allData.celltype_gt), removeTypes_gt);
data_gt        = allData(keepGT,:);
data_gt.celltype = data_gt.celltype_gt;
data_gt = removevars(data_gt, {'celltype_gt','celltype_putative'});
fprintf('  Tabla GT: %d unidades (eliminadas %d Noisy/Undetermined)\n', ...
    height(data_gt), sum(~keepGT));

[data_gt_clean, rep_gt] = flagOutliersPerMeasure(data_gt, ...
    'measureCols', validMeasures, 'logVars', LOG_VARS, ...
    'thresh', OUTLIER_THRESH, 'recCol','ID','subCol','subject', ...
    'groupCol','expcond','verbose',false);

% ── Tabla clasificacion putativa ─────────────────────────────────────────
% Usa celltype_putative; elimina Unknown y etiquetas vacias.
removeTypes_put = {'Unknown',''};
keepPut         = ~ismember(string(allData.celltype_putative), removeTypes_put);
data_put        = allData(keepPut,:);
data_put.celltype = data_put.celltype_putative;
data_put = removevars(data_put, {'celltype_gt','celltype_putative'});
fprintf('  Tabla putativa: %d unidades (eliminadas %d Unknown)\n', ...
    height(data_put), sum(~keepPut));

[data_put_clean, rep_put] = flagOutliersPerMeasure(data_put, ...
    'measureCols', validMeasures, 'logVars', LOG_VARS, ...
    'thresh', OUTLIER_THRESH, 'recCol','ID','subCol','subject', ...
    'groupCol','expcond','verbose',false);

% Guardar ambas tablas limpias
save(fullfile(OUTPUT_DIR,'cellsdata_gt_clean.mat'),  'data_gt_clean',  'rep_gt',  '-v7.3');
save(fullfile(OUTPUT_DIR,'cellsdata_put_clean.mat'), 'data_put_clean', 'rep_put', '-v7.3');
fprintf('  Tablas limpias guardadas en: %s\n', OUTPUT_DIR);

% Exportar a Excel tambien
writetable(data_gt_clean,  fullfile(OUTPUT_DIR,'cellmetrics_gt_clean.xlsx'));
writetable(data_put_clean, fullfile(OUTPUT_DIR,'cellmetrics_put_clean.xlsx'));

%% ========================================================================
%  S5 y S6. ANALISIS DE CONTEOS/RATIOS + LME EN BATCH
%  Se ejecuta una vez por cada sistema de clasificacion:
%    GT       → output/gt/    subgrupos: exc, inh, deep, sup, pv, sst, id2
%    Putativo → output/put/   subgrupos: pyr, inh, wide, narrow
%% ========================================================================

expCondsCell = cellstr(EXP_CONDS);

% Confirmar lista de medidas para batch (numericas, no-sueno, presentes en tabla)
coreCell  = cellstr(CORE_MEASURES);
allVN     = data_gt_clean.Properties.VariableNames;
inTab     = ismember(coreCell, allVN);
isNum     = false(1,numel(coreCell));
for mi = 1:numel(coreCell)
    if inTab(mi), isNum(mi) = isnumeric(data_gt_clean.(coreCell{mi})); end
end
isWake        = ~cellfun(@(n) any(contains(n, EXCLUDE_PATTERNS)), coreCell);
measureNames  = coreCell(inTab & isNum & isWake);
nMeasures     = numel(measureNames);

fprintf('\n  Medidas confirmadas para batch (%d):\n', nMeasures);
for mi = 1:nMeasures, fprintf('    %2d. %s\n', mi, measureNames{mi}); end

% ── Autodetectar etiquetas putativas presentes en los datos ───────────────
% pyr  → cualquier etiqueta con 'Pyramidal'
% inh  → cualquier etiqueta con 'Interneuron'
% wide → subconjunto de inh con 'Wide'
% narr → subconjunto de inh con 'Narrow'
putativeTypes = unique(string(data_put_clean.celltype));
putativeTypes = putativeTypes(putativeTypes ~= "" & putativeTypes ~= "Unknown");

pyrTypes  = putativeTypes(contains(putativeTypes,'Pyramidal','IgnoreCase',true));
intTypes  = putativeTypes(contains(putativeTypes,'Interneuron','IgnoreCase',true));
wideTypes = intTypes(contains(intTypes,'Wide','IgnoreCase',true));
narrTypes = intTypes(contains(intTypes,'Narrow','IgnoreCase',true));

fprintf('\n  Etiquetas putativas detectadas en datos:\n');
fprintf('    Excitatorias (pyr):   %s\n', strjoin(pyrTypes,', '));
fprintf('    Inhibitorias (inh):   %s\n', strjoin(intTypes,', '));
fprintf('    Wide interneuron:     %s\n', strjoin(wideTypes,', '));
fprintf('    Narrow interneuron:   %s\n', strjoin(narrTypes,', '));

% Construir subgrupos putativos dinamicamente
putativeSubgroups = struct('name',{},'types',{});
if ~isempty(pyrTypes)
    putativeSubgroups(end+1).name  = 'pyr';
    putativeSubgroups(end).types   = cellstr(pyrTypes);
end
if ~isempty(intTypes)
    putativeSubgroups(end+1).name  = 'inh';
    putativeSubgroups(end).types   = cellstr(intTypes);
end
if ~isempty(wideTypes)
    putativeSubgroups(end+1).name  = 'wide';
    putativeSubgroups(end).types   = cellstr(wideTypes);
end
if ~isempty(narrTypes)
    putativeSubgroups(end+1).name  = 'narrow';
    putativeSubgroups(end).types   = cellstr(narrTypes);
end

% ── Definir los dos sistemas de clasificacion ────────────────────────────
classSystems = struct( ...
    'name',    {'gt',                            'putative'}, ...
    'label',   {'Clasificacion ground-truth',    'Tipo celular putativo'}, ...
    'data',    {data_gt_clean,                   data_put_clean}, ...
    'subgroups', { ...
        struct( ...
            'name',  {'exc',               'inh',              'deep',        'sup',         'pv',   'sst',   'id2'}, ...
            'types', {{'CAMK2_DEEP','CAMK2_SUP'}, {'PV+','SST+','ID2+'}, {'CAMK2_DEEP'}, {'CAMK2_SUP'}, {'PV+'}, {'SST+'}, {'ID2+'}}), ...
        putativeSubgroups ...
    } ...
);

failLog = {};

for cs = 1:numel(classSystems)
    sys      = classSystems(cs);
    sysDir   = fullfile(OUTPUT_DIR, sys.name);
    batchDir = fullfile(sysDir, 'lme_batch');
    figDir   = fullfile(sysDir, 'figures');
    T        = sys.data;
    sgs      = sys.subgroups;

    fprintf('\n=== Clasificacion: %s ===\n', sys.label);
    for d = {sysDir, batchDir, figDir}
        if ~exist(d{1},'dir'), mkdir(d{1}); end
    end

    % ── PASO 5: Analisis de conteos y ratios ─────────────────────────────
    fprintf('  Paso 5: Conteos y ratios\n');

    analyzeCellCounts(T, GROUP_COLORS, expCondsCell, GROUP_LABELS, ...
        fullfile(sysDir,'counts_total'));

    % analyzeCellTypeRatio y analyzeSubtypeRatio usan etiquetas CAMK2/PV+/SST+/ID2+
    % y solo son validos para la clasificacion GT.
    if strcmp(sys.name,'gt')
        analyzeCellTypeRatio(T, GROUP_COLORS, expCondsCell, GROUP_LABELS, ...
            fullfile(sysDir,'ratio_inh_exc'));
        for st = 1:numel(INH_SUBTYPES)
            subtype = INH_SUBTYPES{st};
            analyzeSubtypeRatio(T, GROUP_COLORS, expCondsCell, GROUP_LABELS, ...
                subtype, fullfile(sysDir, sprintf('ratio_%s', strrep(subtype,'+','plus'))));
        end
    else
        fprintf('    Ratios omitidos (etiquetas putativas no compatibles)\n');
    end

    % Guardar figuras de conteos
    saveFigs(figDir, sprintf('%s_counts', sys.name));

    % ── PASO 6: LME en batch ─────────────────────────────────────────────
    fprintf('  Paso 6: LME en batch (%d subgrupos x %d medidas)\n', numel(sgs), nMeasures);

    results = struct('subgroup',{},'measure',{},'stats',{},'status',{});
    ri      = 0;

    for s = 1:numel(sgs)
        sg    = sgs(s);
        sgDir = fullfile(batchDir, sg.name);
        if ~exist(sgDir,'dir'), mkdir(sgDir); end

        Tsub   = T(ismember(string(T.celltype), sg.types), :);
        nUnits = height(Tsub);
        fprintf('    Subgrupo: %-10s | %d unidades\n', sg.name, nUnits);

        if nUnits == 0
            fprintf('      AVISO: sin unidades — omitiendo\n');
            failLog{end+1} = sprintf('[%s] %s: sin unidades', sys.name, sg.name); %#ok<SAGROW>
            continue
        end

        for m = 1:nMeasures
            mn   = measureNames{m};
            vals = Tsub.(mn);
            ri   = ri + 1;
            results(ri).subgroup = sg.name;
            results(ri).measure  = mn;

            if all(isnan(vals))
                results(ri).status = 'omitido_todoNaN';
                results(ri).stats  = [];
                failLog{end+1} = sprintf('[%s] %s/%s: todo NaN', sys.name, sg.name, mn); %#ok<SAGROW>
                continue
            end

            mDir    = fullfile(sgDir, mn);
            origDir = pwd;
            if ~exist(mDir,'dir'), mkdir(mDir); end
            try
                cd(mDir);
                lmeOut = analyze_ripple_measure_lme_figures_generalized( ...
                    Tsub, mn, 'expcond', GROUP_COLORS, EXP_CONDS, GROUP_LABELS);
                cd(origDir);
                results(ri).stats  = lmeOut;
                results(ri).status = 'ok';

                % Guardar figuras LME
                saveFigs(fullfile(figDir, sg.name), sprintf('%s_%s', sg.name, mn));

            catch ME
                cd(origDir);
                fprintf('      FALLIDO %s/%s: %s\n', sg.name, mn, ME.message);
                results(ri).status = sprintf('fallido: %s', ME.message);
                results(ri).stats  = [];
                failLog{end+1} = sprintf('[%s] %s/%s: %s', sys.name, sg.name, mn, ME.message); %#ok<SAGROW>
            end
        end % medidas
    end % subgrupos

    % ── Reporte resumen ───────────────────────────────────────────────────
    reportFile = fullfile(sysDir, sprintf('reporte_%s.txt', sys.name));
    writeFinalReport(results, sys.label, measureNames, GROUP_LABELS, reportFile);
    fprintf('  Reporte guardado: %s\n', reportFile);

    % Log de incidencias del batch
    fid = fopen(fullfile(batchDir,'incidencias_batch.txt'),'w');
    fprintf(fid,'Clasificacion: %s\nGenerado: %s\n\n', sys.label, datestr(now));
    sysLog = failLog(contains(failLog, sprintf('[%s]', sys.name)));
    if isempty(sysLog)
        fprintf(fid,'Sin incidencias.\n');
    else
        fprintf(fid,'Incidencias (%d):\n', numel(sysLog));
        for i = 1:numel(sysLog), fprintf(fid,'  - %s\n', sysLog{i}); end
    end
    fclose(fid);

end % sistemas de clasificacion

%% ========================================================================
%  RESUMEN FINAL
%% ========================================================================

fprintf('\n========================================================\n');
fprintf('[Pipeline Spike Analysis] Completado.\n');
fprintf('  Grupos procesados   : %s\n', strjoin(groupNames,', '));
fprintf('  Clasificacion GT    : %d unidades tras QC\n', height(data_gt_clean));
fprintf('  Clasificacion put.  : %d unidades tras QC\n', height(data_put_clean));
fprintf('  Medidas analizadas  : %d\n', nMeasures);
fprintf('  Resultados en       : %s\n', OUTPUT_DIR);
fprintf('========================================================\n\n');


%% ========================================================================
%  FUNCIONES AUXILIARES LOCALES
%% ========================================================================

function saveFigs(figFolder, prefix)
%SAVEFIGS  Guarda todas las figuras abiertas como SVG y cierra.
%   Mismo patron que en cell_metrics_processing (saveCellMetricsFigs).

    if ~exist(figFolder,'dir'), mkdir(figFolder); end
    figs = findobj('Type','figure');
    for i = 1:numel(figs)
        f    = figs(i);
        name = sprintf('%s_Fig%d', prefix, f.Number);
        saveas(f, fullfile(figFolder, [name '.svg']));
    end
    close all;

end % saveFigs


function writeFinalReport(results, sysLabel, measureNames, groupLabels, reportFile)
%WRITEFINALREPORT  Escribe un resumen de hallazgos significativos.
%   Solo incluye medidas con omnibus p < 0.10 (tendencia o significativo).
%   p omnibus de stats.contextA.p_value; pairwise de stats.contextA.posthoc.

    SIG_ALPHA   = 0.05;
    TREND_ALPHA = 0.10;

    fid = fopen(reportFile,'w');

    fprintf(fid,'================================================================\n');
    fprintf(fid,'  REPORTE RESUMEN — SPIKE ANALYSIS\n');
    fprintf(fid,'  Clasificacion : %s\n', sysLabel);
    fprintf(fid,'  Generado      : %s\n', datestr(now));
    fprintf(fid,'================================================================\n\n');
    fprintf(fid,'  Significancia: * p<0.05  ** p<0.01  *** p<0.001  (t) tendencia p<0.10\n');
    fprintf(fid,'  Grupos: %s\n\n', strjoin(groupLabels,' | '));
    fprintf(fid,'  Valores p pairwise: correccion Holm.\n\n');

    % Subgrupos unicos en orden de aparicion
    allSubgroups = {};
    for r = 1:numel(results)
        if ~ismember(results(r).subgroup, allSubgroups)
            allSubgroups{end+1} = results(r).subgroup; %#ok<AGROW>
        end
    end

    fprintf(fid,'================================================================\n');
    fprintf(fid,'  HALLAZGOS SIGNIFICATIVOS  (omnibus p < %.2f)\n', TREND_ALPHA);
    fprintf(fid,'================================================================\n\n');

    hasSig = false;
    for sg = 1:numel(allSubgroups)
        sgName    = allSubgroups{sg};
        sgResults = results(strcmp({results.subgroup}, sgName));
        block     = '';

        for r = 1:numel(sgResults)
            res = sgResults(r);
            if ~strcmp(res.status,'ok') || isempty(res.stats), continue; end

            [omniP, pairLines] = extractStats(res.stats);
            if isnan(omniP) || omniP >= TREND_ALPHA, continue; end

            hasSig = true;
            omniSig = getSigStars(omniP, SIG_ALPHA, TREND_ALPHA);
            block   = [block, sprintf('  %-42s omnibus p = %.4f%s\n', ...
                res.measure, omniP, omniSig)]; %#ok<AGROW>
            for pp = 1:numel(pairLines)
                block = [block, sprintf('    %s\n', pairLines{pp})]; %#ok<AGROW>
            end
            block = [block, '\n']; %#ok<AGROW>
        end

        if ~isempty(block)
            fprintf(fid,'  Subgrupo: %s\n', upper(sgName));
            fprintf(fid,'  %s\n', repmat('-',1,60));
            fprintf(fid,'%s', block);
        end
    end

    if ~hasSig
        fprintf(fid,'  Ninguna medida alcanzo p < %.2f en ningun subgrupo.\n', TREND_ALPHA);
    end

    fprintf(fid,'\n================================================================\n');
    fprintf(fid,'  FIN DEL REPORTE\n');
    fprintf(fid,'================================================================\n');
    fclose(fid);

end % writeFinalReport


function [omniP, pairLines] = extractStats(stats)
%EXTRACTSTATS  Extrae p omnibus y lineas pairwise de la estructura LME conocida.
%   stats.contextA.p_value          — p omnibus (escalar)
%   stats.contextA.posthoc          — solo presente si p_value < 0.05
%   stats.contextA.posthoc.contrasts — etiquetas de comparaciones
%   stats.contextA.posthoc.p_holm   — p corregidos Holm
%   stats.contextA.posthoc.p_raw    — p sin corregir (fallback)

    omniP     = NaN;
    pairLines = {};

    if isempty(stats) || ~isstruct(stats), return; end

    if isfield(stats,'contextA') && isfield(stats.contextA,'p_value')
        omniP = double(stats.contextA.p_value);
    end

    if ~isfield(stats,'contextA') || ~isfield(stats.contextA,'posthoc'), return; end
    ph = stats.contextA.posthoc;

    if     isfield(ph,'p_holm'), pvals = ph.p_holm(:);
    elseif isfield(ph,'p_raw'),  pvals = ph.p_raw(:);
    else,  return
    end

    if isfield(ph,'contrasts')
        labels = ph.contrasts;
    else
        labels = arrayfun(@(i) sprintf('Comparacion %d',i),(1:numel(pvals))','UniformOutput',false);
    end

    SIG_ALPHA   = 0.05;
    TREND_ALPHA = 0.10;

    for i = 1:numel(pvals)
        pv  = pvals(i);
        lbl = '';
        if iscell(labels) && i <= numel(labels)
            lbl = labels{i};
            if ~ischar(lbl), lbl = char(string(lbl)); end
        end
        pvSig = getSigStars(pv, SIG_ALPHA, TREND_ALPHA);
        pairLines{end+1} = sprintf('%-36s p_holm=%.4f%s', lbl, pv, pvSig); %#ok<AGROW>
    end

end % extractStats


function stars = getSigStars(p, sigAlpha, trendAlpha)
%GETSIGSTARS  Devuelve la cadena de asteriscos segun nivel de p.

    if     p < 0.001,          stars = ' ***';
    elseif p < 0.01,           stars = ' **';
    elseif p < sigAlpha,       stars = ' *';
    elseif p < trendAlpha,     stars = ' (t)';
    else,                      stars = '';
    end

end % getSigStars


function cellTypes = extractCellTypes(cell_metrics, nCells, system)
%EXTRACTCELLTYPES  Extrae etiquetas de tipo celular segun el sistema indicado.
%
%   system = 'gt'       → ground_truth_classification
%                          (CAMK2_DEEP, CAMK2_SUP, PV+, SST+, ID2+, Noisy_unit, Undetermined)
%   system = 'putative' → putativeCellType
%                          (Pyramidal Cell, Wide Interneuron, Narrow Interneuron)
%
%   Devuelve cell array columna de strings, longitud nCells.

    if nargin < 3, system = 'gt'; end
    cellTypes = [];

    if strcmp(system,'gt')
        if isfield(cell_metrics,'ground_truth_classification')
            gt = cell_metrics.ground_truth_classification;
            if     isstruct(gt) && isfield(gt,'cell_subtypes'), cellTypes = gt.cell_subtypes;
            elseif isstruct(gt) && isfield(gt,'cell_types'),    cellTypes = gt.cell_types;
            elseif isstruct(gt) && isfield(gt,'labels'),        cellTypes = gt.labels;
            end
        end
        if isempty(cellTypes)
            warning('extractCellTypes: ground_truth_classification no encontrado — etiquetando como Unknown.');
            cellTypes = repmat({'Unknown'}, nCells, 1);
        end

    elseif strcmp(system,'putative')
        if isfield(cell_metrics,'putativeCellType')
            cellTypes = cell_metrics.putativeCellType;
        end
        if isempty(cellTypes)
            warning('extractCellTypes: putativeCellType no encontrado — etiquetando como Unknown.');
            cellTypes = repmat({'Unknown'}, nCells, 1);
        end
    end

    % Estandarizar a cell array columna de strings
    if ischar(cellTypes)
        cellTypes = cellstr(cellTypes);
    elseif isnumeric(cellTypes)
        cellTypes = arrayfun(@(x) num2str(x), cellTypes, 'UniformOutput', false);
    end
    cellTypes = cellTypes(:);

    % Guardia ante discordancia de longitud
    if numel(cellTypes) ~= nCells
        warning('extractCellTypes [%s]: discordancia de longitud (%d vs %d) — rellenando con Unknown.', ...
            system, numel(cellTypes), nCells);
        out     = repmat({'Unknown'}, nCells, 1);
        n       = min(numel(cellTypes), nCells);
        out(1:n) = cellTypes(1:n);
        cellTypes = out;
    end

end % extractCellTypes


function varlist = buildVarList(cell_metrics, coreMeasures, excludePatterns)
%BUILDVARLIST  Devuelve el subconjunto de coreMeasures que existe en cell_metrics
%   como campo numerico escalar. NO hay autodeteccion — lista explicita y fija.
%   Medidas NREM/REM se excluyen incondicionalmente (estudio wake-only).

    nCells = cell_metrics.general.cellCount;
    keep   = false(1, numel(coreMeasures));

    for v = 1:numel(coreMeasures)
        vn = coreMeasures(v);

        % Exclusion hard de estados de sueno
        if any(contains(vn, excludePatterns))
            fprintf('      EXCLUIDA (estado sueno): %s\n', vn);
            continue
        end

        if isfield(cell_metrics, vn)
            val = cell_metrics.(vn);
            if isnumeric(val) && numel(val) == nCells
                keep(v) = true;
            else
                fprintf('      NOTA: %s existe pero tipo/tamano incorrecto — omitiendo\n', vn);
            end
        else
            fprintf('      NOTA: %s no encontrado en cell_metrics\n', vn);
        end
    end

    varlist = coreMeasures(keep);

end % buildVarList

function subjectID = extractSubjectID(sesName)
%EXTRACTSUBJECTID  Extrae el ID de sujeto del nombre de sesion.
%
%  Patrones soportados:
%    Shank3male_71wt_250925_1_16ch      → Shank3male_71wt
%    Shank3male_4311wt_250617_one       → Shank3male_4311wt
%    Shank3female_37wt_1_251126_16ch    → Shank3female_37wt
%    Shank3female_31he_2_251125_16ch    → Shank3female_31he

    % Patron 1: corta antes de _N_FECHA  (numero sesion ANTES de fecha de 6 dig)
    %   ej: Shank3female_37wt _1_ 251126_16ch
    tokens = regexp(sesName, '^(.+?)_\d+_\d{6}', 'tokens', 'once');
    if ~isempty(tokens)
        subjectID = tokens{1};
        return
    end

    % Patron 2: corta antes de _FECHA  (numero sesion DESPUES de fecha)
    %   ej: Shank3male_71wt _250925_ 1_16ch
    tokens = regexp(sesName, '^(.+?)_\d{6}', 'tokens', 'once');
    if ~isempty(tokens)
        subjectID = tokens{1};
        return
    end

    % Patron 3: sufijos textuales sin fecha  (_one, _two, _16ch, _malo...)
    sessionSuffixes = {'_one','_two','_three','_four','_five', ...
                       '_16ch','_32ch','_64ch','_malo'};
    s = sesName;
    for i = 1:numel(sessionSuffixes)
        if endsWith_ci(s, sessionSuffixes{i})
            s = s(1 : end-numel(sessionSuffixes{i}));
        end
    end
    s = regexprep(s, '_\d+$', '');
    subjectID = s;

end % extractSubjectID


function tf = endsWith_ci(str, suffix)
%ENDSWITH_CI  Comprueba si str termina con suffix (case-insensitive).
    tf = strcmpi(str(max(1,end-numel(suffix)+1):end), suffix);
end