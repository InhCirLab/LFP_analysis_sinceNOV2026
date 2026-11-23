%% ===========================================================================
%  spike_analysis_pipeline.m
%  Sex hormones project — full spike metrics analysis pipeline
%
%  WORKFLOW:
%    0. CONFIG  — edit this block only
%    1. RIPPLE RECOMPUTE  — rerun ProcessCellMetrics event_metrics per session
%    2. LOAD & BUILD TABLE  — load cell_metrics across all sessions/groups
%    3. CLEAN DATA  — remove noise/undetermined, flag outliers per measure
%    4. COUNT & RATIO ANALYSES  — cell counts, INH/EXC ratio, subtype ratios
%    5. BATCH LME ANALYSIS  — all subgroups x all measures
%
%  DEPENDENCIES:
%    ProcessCellMetrics (CellExplorer)
%    flagOutliersPerMeasure.m
%    analyzeCellCounts.m
%    analyzeCellTypeRatio.m
%    analyzeSubtypeRatio.m
%    analyze_ripple_measure_lme_figures_generalized.m
% ===========================================================================

clear; clc;

%% ===========================================================================
%  0. CONFIG — edit this block
%% ===========================================================================

% --- Group folder paths (one per experimental condition) ---
% Each folder contains session subfolders.
% Leave empty to use uigetdir prompt.
GROUP_FOLDERS = struct( ...
    'SS', '', ...   % e.g. 'C:\data\SS'
    'OS', '', ...   % e.g. 'C:\data\OS'
    'OD', ''  ...   % e.g. 'C:\data\OD'
);

% --- Output root ---
OUTPUT_DIR = '';    % leave empty to prompt

% --- Subject ID lists (used to assign subject column from session folder names) ---
SUBJECTS = struct( ...
    'SS', {{"C089","C090","C091","C199","E075"}}, ...
    'OS', {{"C093","C094","C095","D094"}}, ...
    'OD', {{"D096","D097","E078","E079"}} ...
);

% --- Group display settings ---
EXP_CONDS   = ["SS","OS","OD"];
GROUP_LABELS = {'Ctrl','ORX','ORX+DHT'};
GROUP_COLORS = [0,0,128; 132,197,247; 7,126,149] / 255;

% --- Recompute ripple metrics? ---
% Set true on first run after changing ripple detection.
% Set false once cell_metrics files are up to date.
RECOMPUTE_RIPPLES = true;

% --- Core measure list ---
% Exact set: existing 14 measures + 3 immobility-specific additions + ripple significance.
% NO auto-detection — only these measures are extracted and analyzed.
%
% Grouped by theme:
%   General intrinsic:     ab_ratio, burstIndex_Mizuseki2012, burstIndex_Royer2012,
%                          cv2, firingRate, firingRateCV, thetaModulationIndex, troughToPeak,
%                          population_modIndex
%   Theta / locomotion:    burstIndex_WAKEtheta, cv2_WAKEtheta, firingRate_WAKEtheta
%   Ripple-related:        ripples_modulationIndex, ripples_modulationRatio,
%                          ripples_modulationSignificanceLevel
%   Immobility (no SWR):   burstIndex_QWake_noRipples_ThDt, cv2_QWake_noRipples_ThDt,
%                          firingRate_QWake_noRipples_ThDt
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

% --- Log-transform candidates for outlier detection ---
% Only strictly positive, right-skewed measures.
LOG_VARS = {'firingRate','firingRate_WAKEtheta','firingRate_QWake_noRipples_ThDt'};

% --- Outlier threshold ---
OUTLIER_THRESH = 4.0;

% --- Subtype analysis: inhibitory subtypes to test ---
INH_SUBTYPES = {'PV+','SST+','ID2+','all_inh'};

%% ===========================================================================
%  Resolve paths via prompt if not set
%% ===========================================================================

groups = fieldnames(GROUP_FOLDERS);
for g = 1:numel(groups)
    if isempty(GROUP_FOLDERS.(groups{g}))
        fprintf('Select folder for group %s:\n', groups{g});
        GROUP_FOLDERS.(groups{g}) = uigetdir(pwd, ['Select folder: ' groups{g}]);
    end
end

if isempty(OUTPUT_DIR)
    OUTPUT_DIR = uigetdir(pwd, 'Select output directory');
end
if ~exist(OUTPUT_DIR,'dir'); mkdir(OUTPUT_DIR); end

%% ===========================================================================
%  1. RECOMPUTE RIPPLE METRICS
%  Runs ProcessCellMetrics event_metrics only, using existing ripples.events.mat
%% ===========================================================================

if RECOMPUTE_RIPPLES
    fprintf('\n=== STEP 1: Recomputing ripple metrics ===\n');
    for g = 1:numel(groups)
        grpFolder = GROUP_FOLDERS.(groups{g});
        sessions  = getSessionFolders(grpFolder);
        fprintf('  Group %s: %d sessions\n', groups{g}, numel(sessions));
        for s = 1:numel(sessions)
            sesPath = fullfile(grpFolder, sessions{s});
            rippleFile = dir(fullfile(sesPath, '*.ripples.events.mat'));
            if isempty(rippleFile)
                fprintf('    SKIP %s — no ripples.events.mat found\n', sessions{s});
                continue
            end
            fprintf('    Processing %s ... ', sessions{s});
            try
                prevDir = pwd; cd(sesPath);

                % Whitelist: only use the 6to12 theta states file.
                % Dynamically build ignoreStateTypes from whatever else exists.
                allStates   = dir(fullfile(sesPath,'*.states.mat'));
                allStates   = {allStates.name};
                keepPattern = '6to12';
                ignoreList  = allStates(~contains(allStates, keepPattern));
                % Strip to state type name only (between first and last dot)
                ignoreTypes = cellfun(@(f) strsplit(f,'.'), ignoreList, ...
                    'UniformOutput',false);
                ignoreTypes = cellfun(@(p) p{end-2}, ignoreTypes, ...
                    'UniformOutput',false);
                ignoreTypes = unique(ignoreTypes);

                ProcessCellMetrics( ...
                    'metrics',         {'event_metrics','state_metrics'}, ...
                    'forceReload',     true, ...
                    'ignoreStateTypes', ignoreTypes, ...
                    'saveMat',         true, ...
                    'showFigures',     false, ...
                    'showGUI',         false, ...
                    'sessionSummaryFigure', false);
                cd(prevDir);
                fprintf('OK\n');
            catch ME
                if exist('prevDir','var'); cd(prevDir); end
                fprintf('FAILED: %s\n', ME.message);
            end
        end
    end
else
    fprintf('\n=== STEP 1: Ripple recompute skipped (RECOMPUTE_RIPPLES=false) ===\n');
end

%% ===========================================================================
%  2. LOAD AND BUILD TABLE
%  Iterates all group folders/sessions, loads cell_metrics, spikes,
%  hippocampalLayers, builds unified table with subject column.
%  Auto-detects available state-specific metrics on first load.
%% ===========================================================================

fprintf('\n=== STEP 2: Loading cell metrics ===\n');

allData       = [];
varlist       = [];   % finalised after first successful load
layersnames   = ["oriens","pyr","rad","slm"];

for g = 1:numel(groups)
    grpName   = groups{g};
    grpFolder = GROUP_FOLDERS.(grpName);
    sessions  = getSessionFolders(grpFolder);
    subjIDs   = SUBJECTS.(grpName);

    fprintf('  Group %s: %d sessions\n', grpName, numel(sessions));

    for s = 1:numel(sessions)
        sesName  = sessions{s};
        sesPath  = fullfile(grpFolder, sesName);

        cmFile  = fullfile(sesPath, [sesName '.cell_metrics.cellinfo.mat']);
        spkFile = fullfile(sesPath, [sesName '.spikes.cellinfo.mat']);
        layFile = fullfile(sesPath, [sesName '.hippocampalLayers.channelinfo.mat']);

        if ~exist(cmFile,'file')
            fprintf('    SKIP %s — no cell_metrics file\n', sesName);
            continue
        end

        try
            load(cmFile,  'cell_metrics');
            load(spkFile, 'spikes');
            load(layFile, 'hippocampalLayers');
        catch ME
            fprintf('    SKIP %s — load error: %s\n', sesName, ME.message);
            continue
        end

        % --- Finalise varlist on first successful load ---
        if isempty(varlist)
            varlist = buildVarList(cell_metrics, CORE_MEASURES);
            fprintf('    Measure list finalised: %d / %d measures found\n', ...
                numel(varlist), numel(CORE_MEASURES));
            logVarsExtended = LOG_VARS;   % fixed — no auto-detection
        end

        % --- Extract measures ---
        nCells = cell_metrics.general.cellCount;
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

        % --- Cell type: extract BOTH classification systems ---
        % celltype_gt      : ground_truth_classification (CAMK2_DEEP, CAMK2_SUP, PV+, SST+, ID2+)
        % celltype_putative: putativeCellType (Pyramidal Cell, Wide Interneuron, Narrow Interneuron)
        % These are kept separate and never mixed.
        thisfileclass_gt  = extractCellTypes(cell_metrics, nCells, 'gt');
        thisfileclass_put = extractCellTypes(cell_metrics, nCells, 'putative');

        % --- Spike channel and layer assignment ---
        spikech     = spikes.phy_maxWaveformCh1(:);
        layers      = hippocampalLayers.layers{1,1}.all;
        spikelayers = zeros(nCells,1);
        for kk = 1:nCells
            [~, spikelayers(kk)] = min(abs(layers - spikech(kk)));
        end
        spikelayersnames = strings(nCells,1);
        spikelayersnames(spikelayers==1) = layersnames(1);
        spikelayersnames(spikelayers==2) = layersnames(2);
        spikelayersnames(spikelayers==3) = layersnames(3);
        spikelayersnames(spikelayers==4) = layersnames(4);

        % --- Subject assignment ---
        thisSubject = repmat("unknown", nCells, 1);
        for ss = 1:numel(subjIDs)
            if contains(sesName, subjIDs{ss})
                thisSubject(:) = subjIDs{ss};
                break
            end
        end

        % --- Build per-session table ---
        thisID       = repmat(string(sesName),  nCells, 1);
        thisExpcond  = repmat(string(grpName),  nCells, 1);
        thisMeasTab  = array2table(thismeasures, 'VariableNames', cellstr(varlist));
        thisMetaTab  = table(thisID, thisExpcond, ...
                             thisfileclass_gt, thisfileclass_put, ...
                             spikech, spikelayersnames, thisSubject, ...
                             'VariableNames', ["ID","expcond","celltype_gt", ...
                                              "celltype_putative","spikech","layer","subject"]);
        thisRec  = [thisMetaTab, thisMeasTab];
        allData  = [allData; thisRec]; %#ok<AGROW>

        fprintf('    Loaded %s: %d cells\n', sesName, nCells);
    end
end

fprintf('  Total units loaded: %d\n', height(allData));

% Save raw table
rawFile = fullfile(OUTPUT_DIR, 'all_cellmetrics_raw.xlsx');
writetable(allData, rawFile);
fprintf('  Raw table saved: %s\n', rawFile);

%% ===========================================================================
%  3. CLEAN DATA
%  Two parallel clean tables — one per classification system.
%  GT:      removes Noisy_unit, Undetermined
%  Putative: removes Unknown, empty labels
%  Outlier flagging applied to both using the same measure list.
%% ===========================================================================

fprintf('\n=== STEP 3: Cleaning data ===\n');

EXCLUDE_PATTERNS = {'NREM','REM','Sleep','sleep','nrem','rem'};

% --- Measure list (shared by both classification systems) ---
validMeasures = cellstr(CORE_MEASURES);
validMeasures = validMeasures( ...
    ismember(validMeasures, allData.Properties.VariableNames) & ...
    ~cellfun(@(n) any(contains(n, EXCLUDE_PATTERNS)), validMeasures));

fprintf('  Measures entering analysis (%d):\n', numel(validMeasures));
for mi = 1:numel(validMeasures)
    fprintf('    %2d. %s\n', mi, validMeasures{mi});
end

% --- GT classification table ---
% Uses celltype_gt; removes Noisy_unit and Undetermined.
% Rename celltype_gt → celltype so downstream functions work unchanged.
removeTypes_gt = {'Noisy_unit','Undetermined'};
keepGT  = ~ismember(string(allData.celltype_gt), removeTypes_gt);
data_gt = allData(keepGT, :);
data_gt.celltype = data_gt.celltype_gt;   % active celltype column
data_gt = removevars(data_gt, {'celltype_gt','celltype_putative'});
fprintf('  GT table: %d units (removed %d Noisy/Undetermined)\n', ...
    height(data_gt), sum(~keepGT));

[data_gt_clean, rep_gt] = flagOutliersPerMeasure(data_gt, ...
    'measureCols', validMeasures, 'logVars', LOG_VARS, ...
    'thresh', OUTLIER_THRESH, 'recCol','ID','subCol','subject', ...
    'groupCol','expcond','verbose',false);

% --- Putative classification table ---
% Uses celltype_putative; removes Unknown/empty.
removeTypes_put = {'Unknown',''};
keepPut  = ~ismember(string(allData.celltype_putative), removeTypes_put);
data_put = allData(keepPut, :);
data_put.celltype = data_put.celltype_putative;
data_put = removevars(data_put, {'celltype_gt','celltype_putative'});
fprintf('  Putative table: %d units (removed %d Unknown)\n', ...
    height(data_put), sum(~keepPut));

[data_put_clean, rep_put] = flagOutliersPerMeasure(data_put, ...
    'measureCols', validMeasures, 'logVars', LOG_VARS, ...
    'thresh', OUTLIER_THRESH, 'recCol','ID','subCol','subject', ...
    'groupCol','expcond','verbose',false);

% Save both
save(fullfile(OUTPUT_DIR,'cellsdata_gt_clean.mat'),  'data_gt_clean',  'rep_gt',  '-v7.3');
save(fullfile(OUTPUT_DIR,'cellsdata_put_clean.mat'), 'data_put_clean', 'rep_put', '-v7.3');
fprintf('  Clean tables saved to output dir.\n');

%% ===========================================================================
%  4 & 5. COUNT/RATIO ANALYSES + BATCH LME
%  Run once for each classification system.
%  GT      → output/gt/    subgroups: exc, inh, deep, sup, pv, sst, id2
%  Putative → output/put/   subgroups: pyr, inh, wide, narrow
%% ===========================================================================

expCondsCell = cellstr(EXP_CONDS);

% Confirmed measure list for batch (same for both systems)
coreCell   = cellstr(CORE_MEASURES);
allVN      = data_gt_clean.Properties.VariableNames;  % same columns in both
inTab      = ismember(coreCell, allVN);
isNum      = false(1,numel(coreCell));
for mi = 1:numel(coreCell)
    if inTab(mi); isNum(mi) = isnumeric(data_gt_clean.(coreCell{mi})); end
end
isWake     = ~cellfun(@(n) any(contains(n, EXCLUDE_PATTERNS)), coreCell);
measureNames = coreCell(inTab & isNum & isWake);
nMeasures    = numel(measureNames);

fprintf('\n  Confirmed measure list for batch (%d):\n', nMeasures);
for mi = 1:nMeasures; fprintf('    %2d. %s\n', mi, measureNames{mi}); end

% Auto-detect putative cell type labels from the actual data.
% Map to exc/inh semantically: anything containing 'Pyramidal' → excitatory;
% anything containing 'Interneuron' or 'interneuron' → inhibitory subtype.
putativeTypes = unique(string(data_put_clean.celltype));
putativeTypes = putativeTypes(putativeTypes ~= "" & putativeTypes ~= "Unknown");

pyrTypes  = putativeTypes(contains(putativeTypes, 'Pyramidal', 'IgnoreCase', true));
intTypes  = putativeTypes(contains(putativeTypes, 'Interneuron', 'IgnoreCase', true));
wideTypes = intTypes(contains(intTypes, 'Wide', 'IgnoreCase', true));
narrTypes = intTypes(contains(intTypes, 'Narrow', 'IgnoreCase', true));

fprintf('\n  Putative cell type labels detected in data:\n');
fprintf('    Excitatory (pyr): %s\n', strjoin(pyrTypes, ', '));
fprintf('    Inhibitory (inh): %s\n', strjoin(intTypes, ', '));
fprintf('    Wide interneuron: %s\n', strjoin(wideTypes, ', '));
fprintf('    Narrow interneuron: %s\n', strjoin(narrTypes, ', '));

% Build putative subgroups from detected labels
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

% Probe LME output field names on first available subgroup/measure
% so extractStatsFromResult knows which fields to look for.
% Result stored in LME_FIELD_PROBE and passed to writeFinalReport.
LME_FIELD_PROBE = [];
probeT = data_gt_clean(ismember(string(data_gt_clean.celltype), ...
    {'CAMK2_DEEP','CAMK2_SUP'}), :);
if height(probeT) > 0 && ~isempty(measureNames)
    for pm = 1:numel(measureNames)
        if ~all(isnan(probeT.(measureNames{pm})))
            try
                origDir = pwd;
                tmpDir  = fullfile(OUTPUT_DIR, 'tmp_probe');
                if ~exist(tmpDir,'dir'); mkdir(tmpDir); end
                cd(tmpDir);
                probeOut = analyze_ripple_measure_lme_figures_generalized( ...
                    probeT, measureNames{pm}, 'expcond', ...
                    GROUP_COLORS, EXP_CONDS, GROUP_LABELS);
                cd(origDir);
                if ~isempty(probeOut)
                    LME_FIELD_PROBE = probeOut;
                    fprintf('\n  LME output fields detected: %s\n', ...
                        strjoin(fieldnames(probeOut), ', '));
                end
                rmdir(tmpDir,'s');
            catch
                cd(origDir);
            end
            break
        end
    end
end

% Define both classification systems using auto-detected putative subgroups
classSystems = struct( ...
    'name',     {'gt',                                  'putative'}, ...
    'label',    {'Ground-truth classification',         'Putative cell type'}, ...
    'data',     {data_gt_clean,                         data_put_clean}, ...
    'subgroups', { ...
        struct('name',  {'exc',                 'inh',              'deep',         'sup',          'pv',   'sst',   'id2'}, ...
               'types', {{'CAMK2_DEEP','CAMK2_SUP'}, {'PV+','SST+','ID2+'}, {'CAMK2_DEEP'}, {'CAMK2_SUP'}, {'PV+'}, {'SST+'}, {'ID2+'}}), ...
        putativeSubgroups ...
    } ...
);

failLog = {};

for cs = 1:numel(classSystems)
    sys      = classSystems(cs);
    sysDir   = fullfile(OUTPUT_DIR, sys.name);
    batchDir = fullfile(sysDir, 'lme_batch');
    T        = sys.data;
    sgs      = sys.subgroups;

    fprintf('\n=== Classification: %s ===\n', sys.label);
    if ~exist(sysDir,'dir'); mkdir(sysDir); end

    % -- Step 4: count and ratio analyses --
    fprintf('  Step 4: Count and ratio analyses\n');

    analyzeCellCounts(T, GROUP_COLORS, expCondsCell, GROUP_LABELS, ...
        fullfile(sysDir, 'counts_total'));

    % analyzeCellTypeRatio and analyzeSubtypeRatio use hardcoded
    % CAMK2/PV+/SST+/ID2+ labels — only valid for GT classification.
    if strcmp(sys.name, 'gt')
        analyzeCellTypeRatio(T, GROUP_COLORS, expCondsCell, GROUP_LABELS, ...
            fullfile(sysDir, 'ratio_inh_exc'));
        for st = 1:numel(INH_SUBTYPES)
            subtype = INH_SUBTYPES{st};
            analyzeSubtypeRatio(T, GROUP_COLORS, expCondsCell, GROUP_LABELS, ...
                subtype, fullfile(sysDir, sprintf('ratio_%s', strrep(subtype,'+','plus'))));
        end
    else
        fprintf('    Skipping ratio analyses (putative labels not compatible)\n');
    end

    % -- Step 5: batch LME --
    fprintf('  Step 5: Batch LME (%d subgroups x %d measures)\n', numel(sgs), nMeasures);
    if ~exist(batchDir,'dir'); mkdir(batchDir); end

    % Pre-allocate results collector
    % results(s,m).stats = whatever analyze_ripple_measure_lme_figures_generalized returns
    results = struct('subgroup',{},'measure',{},'stats',{},'status',{});
    ri = 0;

    for s = 1:numel(sgs)
        sg    = sgs(s);
        sgDir = fullfile(batchDir, sg.name);
        if ~exist(sgDir,'dir'); mkdir(sgDir); end

        Tsub   = T(ismember(string(T.celltype), sg.types), :);
        nUnits = height(Tsub);
        fprintf('    Subgroup: %-8s | %d units\n', sg.name, nUnits);

        if nUnits == 0
            fprintf('      WARNING: no units — skipping\n');
            failLog{end+1} = sprintf('[%s] %s: no units', sys.name, sg.name); %#ok<SAGROW>
            continue
        end

        for m = 1:nMeasures
            mn   = measureNames{m};
            vals = Tsub.(mn);
            ri   = ri + 1;
            results(ri).subgroup = sg.name;
            results(ri).measure  = mn;

            if all(isnan(vals))
                results(ri).status = 'skipped_allNaN';
                results(ri).stats  = [];
                failLog{end+1} = sprintf('[%s] %s/%s: all NaN', sys.name, sg.name, mn); %#ok<SAGROW>
                continue
            end

            mDir = fullfile(sgDir, mn);
            if ~exist(mDir,'dir'); mkdir(mDir); end
            origDir = pwd;
            try
                cd(mDir);
                lmeOut = analyze_ripple_measure_lme_figures_generalized( ...
                    Tsub, mn, 'expcond', GROUP_COLORS, EXP_CONDS, GROUP_LABELS);
                cd(origDir);
                results(ri).stats  = lmeOut;
                results(ri).status = 'ok';
            catch ME
                cd(origDir);
                fprintf('      FAILED %s/%s: %s\n', sg.name, mn, ME.message);
                results(ri).status = sprintf('failed: %s', ME.message);
                results(ri).stats  = [];
                failLog{end+1} = sprintf('[%s] %s/%s: %s', sys.name, sg.name, mn, ME.message); %#ok<SAGROW>
            end
        end
    end

    % -- Final summary report --
    reportFile = fullfile(sysDir, sprintf('summary_report_%s.txt', sys.name));
    writeFinalReport(results, sys.label, measureNames, GROUP_LABELS, reportFile, LME_FIELD_PROBE);
    fprintf('  Summary report: %s\n', reportFile);

    % Batch issues log
    fid = fopen(fullfile(batchDir,'batch_issues.txt'),'w');
    fprintf(fid,'Classification: %s\nGenerated: %s\n\n', sys.label, datestr(now));
    sysLog = failLog(contains(failLog, sprintf('[%s]', sys.name)));
    if isempty(sysLog)
        fprintf(fid,'No issues.\n');
    else
        fprintf(fid,'Issues (%d):\n', numel(sysLog));
        for i = 1:numel(sysLog); fprintf(fid,'  - %s\n', sysLog{i}); end
    end
    fclose(fid);

end

fprintf('\nPipeline complete.\n');


%% ===========================================================================
%  LOCAL FUNCTIONS
%% ===========================================================================

function writeFinalReport(results, sysLabel, measureNames, groupLabels, reportFile, ~)
% Write Section 1 only — significant findings overview.
% Omnibus p from stats.p_value, pairwise from stats.posthoc.

SIG_ALPHA   = 0.05;
TREND_ALPHA = 0.10;

fid = fopen(reportFile, 'w');

fprintf(fid,'================================================================\n');
fprintf(fid,'  SPIKE ANALYSIS SUMMARY REPORT\n');
fprintf(fid,'  Classification : %s\n', sysLabel);
fprintf(fid,'  Generated      : %s\n', datestr(now));
fprintf(fid,'================================================================\n\n');
fprintf(fid,'  Significance: * p<0.05  ** p<0.01  *** p<0.001  (t) trend p<0.10\n');
fprintf(fid,'  Groups: %s\n\n', strjoin(groupLabels, ' | '));
fprintf(fid,'  Pairwise p-values: Holm-corrected.\n\n');

% Unique subgroups in order
allSubgroups = {};
for r = 1:numel(results)
    if ~ismember(results(r).subgroup, allSubgroups)
        allSubgroups{end+1} = results(r).subgroup; %#ok<AGROW>
    end
end

fprintf(fid,'================================================================\n');
fprintf(fid,'  SIGNIFICANT FINDINGS  (omnibus p < %.2f)\n', TREND_ALPHA);
fprintf(fid,'================================================================\n\n');

hasSig = false;
for sg = 1:numel(allSubgroups)
    sgName    = allSubgroups{sg};
    sgResults = results(strcmp({results.subgroup}, sgName));
    block     = '';

    for r = 1:numel(sgResults)
        res = sgResults(r);
        if ~strcmp(res.status,'ok') || isempty(res.stats); continue; end

        [omniP, pairLines] = extractStats(res.stats);
        if isnan(omniP) || omniP >= TREND_ALPHA; continue; end

        hasSig = true;
        % inline sig stars to avoid cross-function scope issue
        if     omniP < 0.001; omniSig = ' ***';
        elseif omniP < 0.01;  omniSig = ' **';
        elseif omniP < 0.05;  omniSig = ' *';
        elseif omniP < 0.10;  omniSig = ' (t)';
        else;                 omniSig = '';
        end
        block  = [block, sprintf('  %-42s omnibus p = %.4f%s\n', ...
            res.measure, omniP, omniSig)]; %#ok<AGROW>
        for pp = 1:numel(pairLines)
            block = [block, sprintf('    %s\n', pairLines{pp})]; %#ok<AGROW>
        end
        block = [block, '\n']; %#ok<AGROW>
    end

    if ~isempty(block)
        fprintf(fid,'  Subgroup: %s\n', upper(sgName));
        fprintf(fid,'  %s\n', repmat('-',1,60));
        fprintf(fid,'%s', block);
    end
end

if ~hasSig
    fprintf(fid,'  No measures reached p < %.2f in any subgroup.\n', TREND_ALPHA);
end

fprintf(fid,'\n================================================================\n');
fprintf(fid,'  END OF REPORT\n');
fprintf(fid,'================================================================\n');
fclose(fid);
end % writeFinalReport


function [omniP, pairLines] = extractStats(stats)
% Extract omnibus p and pairwise lines from the known output structure:
%   stats.contextA.p_value          — omnibus p (scalar)
%   stats.contextA.posthoc          — only present when p_value < 0.05
%   stats.contextA.posthoc.contrasts — cell array of labels {nPairs x 1}
%   stats.contextA.posthoc.p_holm   — Holm-corrected p-values [nPairs x 1]
%   stats.contextA.posthoc.p_raw    — raw p-values [nPairs x 1]

omniP     = NaN;
pairLines = {};

if isempty(stats) || ~isstruct(stats); return; end

% Omnibus p — always in contextA
if isfield(stats,'contextA') && isfield(stats.contextA,'p_value')
    omniP = double(stats.contextA.p_value);
end

% Pairwise — only present when omnibus p < 0.05
if ~isfield(stats,'contextA') || ~isfield(stats.contextA,'posthoc')
    return
end
ph = stats.contextA.posthoc;

% p-values: prefer Holm-corrected
if isfield(ph,'p_holm')
    pvals = ph.p_holm(:);
elseif isfield(ph,'p_raw')
    pvals = ph.p_raw(:);
else
    return
end

% labels
if isfield(ph,'contrasts')
    labels = ph.contrasts;
else
    labels = arrayfun(@(i) sprintf('Comparison %d',i), ...
        (1:numel(pvals))', 'UniformOutput', false);
end

for i = 1:numel(pvals)
    pv  = pvals(i);
    lbl = '';
    if iscell(labels) && i <= numel(labels)
        lbl = labels{i};
        if ~ischar(lbl); lbl = char(string(lbl)); end
    end
    if     pv < 0.001; pvSig = ' ***';
    elseif pv < 0.01;  pvSig = ' **';
    elseif pv < 0.05;  pvSig = ' *';
    elseif pv < 0.10;  pvSig = ' (t)';
    else;              pvSig = '';
    end
    pairLines{end+1} = sprintf('%-36s p_holm=%.4f%s', lbl, pv, pvSig); %#ok<AGROW>
end
end % extractStats


function [omniField, pairField, pairPField, pairLField] = detectLMEFields(~)
% Field names are now hardcoded from the known output structure.
% This function is kept for interface compatibility but returns empty —
% extractStats uses direct field access instead.
omniField  = '';
pairField  = '';
pairPField = '';
pairLField = '';
end % detectLMEFields

function cellTypes = extractCellTypes(cell_metrics, nCells, system)
% Extract cell type labels from cell_metrics for a specific classification system.
%
% system = 'gt'       → ground_truth_classification.cell_subtypes
%                        (CAMK2_DEEP, CAMK2_SUP, PV+, SST+, ID2+, Noisy_unit, Undetermined)
% system = 'putative' → putativeCellType
%                        (Pyramidal Cell, Wide Interneuron, Narrow Interneuron)
%
% Returns a cell column array of strings, length nCells.

    if nargin < 3; system = 'gt'; end
    cellTypes = [];

    if strcmp(system, 'gt')
        % Ground truth classification (HippoCookBook cellTypeClassifier output)
        if isfield(cell_metrics,'ground_truth_classification')
            gt = cell_metrics.ground_truth_classification;
            if isstruct(gt) && isfield(gt,'cell_subtypes')
                cellTypes = gt.cell_subtypes;
            elseif isstruct(gt) && isfield(gt,'cell_types')
                cellTypes = gt.cell_types;
            elseif isstruct(gt) && isfield(gt,'labels')
                cellTypes = gt.labels;
            end
        end
        if isempty(cellTypes)
            warning('extractCellTypes: ground_truth_classification not found — labelling all as Unknown.');
            cellTypes = repmat({'Unknown'}, nCells, 1);
        end

    elseif strcmp(system, 'putative')
        % Standard CellExplorer putative cell type
        if isfield(cell_metrics,'putativeCellType')
            cellTypes = cell_metrics.putativeCellType;
        end
        if isempty(cellTypes)
            warning('extractCellTypes: putativeCellType not found — labelling all as Unknown.');
            cellTypes = repmat({'Unknown'}, nCells, 1);
        end
    end

    % Standardise to cell column array of strings
    if ischar(cellTypes)
        cellTypes = cellstr(cellTypes);
    elseif isnumeric(cellTypes)
        cellTypes = arrayfun(@(x) num2str(x), cellTypes, 'UniformOutput', false);
    end
    cellTypes = cellTypes(:);

    % Guard length mismatch
    if numel(cellTypes) ~= nCells
        warning('extractCellTypes [%s]: length mismatch (%d vs %d) — padding with Unknown.', ...
            system, numel(cellTypes), nCells);
        out = repmat({'Unknown'}, nCells, 1);
        n   = min(numel(cellTypes), nCells);
        out(1:n) = cellTypes(1:n);
        cellTypes = out;
    end
end


function sessions = getSessionFolders(rootPath)
% Returns names of immediate subfolders (= recording sessions).
    d = dir(rootPath);
    d = d([d.isdir] & ~ismember({d.name},{'.','..'}));
    sessions = {d.name};
end


function varlist = buildVarList(cell_metrics, coreMeasures)
% Returns the subset of coreMeasures that exist as numeric scalar fields
% in cell_metrics. NO auto-detection — list is explicit and fixed.
% NREM/REM measures are excluded unconditionally (wake-only study).
% Missing measures are reported as notes but do not abort.

    EXCLUDE_PATTERNS = {'NREM','REM','Sleep','sleep','nrem','rem'};

    nCells  = cell_metrics.general.cellCount;
    keep    = false(1, numel(coreMeasures));

    for v = 1:numel(coreMeasures)
        vn = coreMeasures(v);

        % Hard exclude NREM/REM
        if any(contains(vn, EXCLUDE_PATTERNS))
            fprintf('      EXCLUDED (sleep state): %s\n', vn);
            continue
        end

        if isfield(cell_metrics, vn)
            val = cell_metrics.(vn);
            if isnumeric(val) && numel(val) == nCells
                keep(v) = true;
            else
                fprintf('      NOTE: %s exists but wrong type/size — skipping\n', vn);
            end
        else
            fprintf('      NOTE: %s not found in cell_metrics\n', vn);
        end
    end

    varlist = coreMeasures(keep);
end
