%% ========================================================================
%  LFP ANALYSIS PIPELINE — Post-HippoCookbook preprocessing (REORDERED)
%  ========================================================================
%
%  ORIGINAL PROBLEM
%    intervals_nonTheta_nonRipples (used by Part I nonrun and by
%    ps1to120 -> allrec / loc / nonlocnonrip) is computed in detectThetaEpochs
%    by READING whatever <rec>.ripples.events.mat exists at that moment. In
%    the old order (I -> III -> II) the new ripples from Part III were
%    generated AFTER saving .thetaEpochs6to12.states.mat, so Part I (nonrun)
%    and Part II used the OLD ripples.
%
%  NEW EXECUTION ORDER
%    1. P1  PASS 1    : theta + run (generates thetaRun.newintervals, which
%                       the CNN needs to exclude locomotion). No Excel.
%    2. P3  RIPPLES   : CNN detection -> saves NEW <rec>.ripples.events.mat
%    3. P1R REFRESH   : re-runs detectThetaEpochs (now it reads the new
%                       ripples), recomputes nonrun, rewrites
%                       .thetaEpochs6to12.states.mat and saves Excel/figures.
%                       Prints per session whether intervals_nonTheta_nonRipples
%                       changed relative to the old version.
%    4. P2  POWER SPECTRUM : ps1to120_Paria reads the already-refreshed .states.mat
%
%  ps1to120_Paria.m needs NO changes.
%
%  The helper functions (extractEpochPower, storeRunPower, buildRow,
%  saveFigSession, saveGroupFig, savePSFig, robust_zscore_matrix,
%  getAnimalID ...) are the usual ones, on the path.
%
%  --- LOCAL CHANGE (Paria) -------------------------------------------------
%  P3.5: the animal ID is extracted with getAnimalID (the last chunk after
%  the final '-'), e.g. 2026-06-09_16-47-11-F016 -> F016. Previously it used
%  extractBefore(ID,'_'), which returned "2026" for these names and broke the
%  hierarchical outlier removal (all sessions treated as one animal).
% =========================================================================

%clearvars; close all; clc;

%% ── SERIAL EXECUTION ─────────────────────────────────────────────────────
%  This version runs EVERYTHING serially (one recording / one session at a
%  time). No parallel pool is opened. Slower, but easier to follow and debug:
%  you see each session as it runs, in order, and get the full error for the
%  exact session that fails.
%
%  To switch back to parallel later:
%    1) reopen a pool here, e.g.
%         N_WORKERS = 8;
%         if license('test','Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
%             parpool('local', min(N_WORKERS, maxNumCompThreads));
%         end
%    2) change the two loops marked "SERIAL (was parfor)" back to parfor:
%         - in runThetaStage below
%         - in ps1to120_Paria.m

%% ========================================================================
%  PART I — THETA & GAMMA DURING LOCOMOTION
%% ========================================================================

%% ── P1.0  PARAMETERS ─────────────────────────────────────────────────────

P.THETA_BAND  = [4  12];
P.LGAMMA_BAND = [20 60];
P.HGAMMA_BAND = [60 100];
P.N_CYCLES    = 8;
P.MIN_RUN_DUR = 1;   % seconds — running epochs shorter than this are discarded

P.DETECT_OPTS = {'saveMat', false, 'force', true, ...
                 'bandpass', [6 12], 'delta_bandpass', [1 3.5], ...
                 'powerThreshold', 1};

P.powerVarNames = ["events (n)", "totaldur (s)", "meandur (s)", ...
    "medianTheta (uv)",  "medianLgamma (uv)",  "medianHgamma (uv)", ...
    "medianTheta (db)",  "medianLgamma (db)",  "medianHgamma (db)", ...
    "zTheta (au)",       "zLgamma (au)",       "zHgamma (au)",      ...
    "zmodTheta (au)",    "zmodLgamma (au)",    "zmodHgamma (au)"];

%% ── P1.1  CONFIGURATION: GROUPS AND FOLDERS ──────────────────────────────

nGroupsStr = inputdlg('How many experimental groups do you have?', 'Groups', 1, {'1'});
if isempty(nGroupsStr), error('Aborted.'); end
nGroups = str2double(nGroupsStr{1});
assert(~isnan(nGroups) && nGroups >= 1, 'Invalid number of groups.');

groupFolders = cell(nGroups, 1);
groupLabels  = cell(nGroups, 1);

for g = 1:nGroups
    lbl = inputdlg(sprintf('Label for group %d (e.g. PV, SST, WT...):', g), ...
                   sprintf('Group %d', g), 1, {sprintf('G%d', g)});
    if isempty(lbl), error('Aborted.'); end
    groupLabels{g} = strtrim(lbl{1});

    fprintf('>>> Select the folder for group "%s":\n', groupLabels{g});
    folder = uigetdir(pwd, sprintf('Group folder: %s', groupLabels{g}));
    if folder == 0, error('No folder selected. Aborted.'); end
    groupFolders{g} = folder;
end

fprintf('>>> Select the OUTPUT folder for Part I results:\n');
outFolder = uigetdir(groupFolders{1}, 'OUTPUT — Part I');
if outFolder == 0, error('Aborted.'); end

figFolder = fullfile(outFolder, 'figures');
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

%% ── P1.2  PASS 1: THETA + RUN (only to generate newintervals) ───────────
%  Needed because the CNN (Part III) excludes thetaRun.newintervals.
%  Does not save Excel or figures: that happens in P1R with the new ripples.

fprintf('\n##### PASS 1 (theta + run) #####\n');
runThetaStage(groupFolders, groupLabels, P, figFolder, true);

close all
%% ========================================================================
%  PART III — RIPPLE ANALYSIS
%% ========================================================================

%% ── P3.0  PARAMETERS ─────────────────────────────────────────────────────

% CNN model and Python executable paths
CNN_MODEL_PATH = 'C:\Users\usuario\cnn-matlab-master\cnn';
CNN_PYTHON_EXE = 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe';

% Detection parameters
RIPPLE_THRESHOLD = [2 5];        % [min max] in SD
RIPPLE_PASSBAND  = [80 240];     % Hz
PEAKNORM_CUTOFF  = 50;           % outliers with PeakNormPow > this value are discarded

% Robust outlier analysis — thresholds
OUTLIER_THRESH_MAD  = 3.5;
OUTLIER_MIN_N_REC   = 20;
OUTLIER_MIN_N_SUB   = 20;
OUTLIER_MIN_N_GROUP = 20;
OUTLIER_VARS        = {'PeakAmp','PeakNormPow','Duration','PeakFq'};
OUTLIER_LOG_VARS    = {'PeakAmp','PeakNormPow','Duration'};

%% ── P3.1  CONFIGURATION: GROUPS AND FOLDERS ──────────────────────────────

nGroupsStr3 = inputdlg('How many groups for the ripple analysis?', 'Ripples — Groups', 1, {num2str(nGroups)});
if isempty(nGroupsStr3), error('Aborted.'); end
nGroupsR = str2double(nGroupsStr3{1});
assert(~isnan(nGroupsR) && nGroupsR >= 1, 'Invalid number of groups.');

ripGroupFolders = cell(nGroupsR, 1);
ripGroupLabels  = cell(nGroupsR, 1);

for g = 1:nGroupsR
    defaultLbl = '';
    if g <= nGroups, defaultLbl = groupLabels{g}; else, defaultLbl = sprintf('G%d', g); end
    lbl3 = inputdlg(sprintf('Label for group %d:', g), ...
                    sprintf('Ripples group %d', g), 1, {defaultLbl});
    if isempty(lbl3), error('Aborted.'); end
    ripGroupLabels{g} = strtrim(lbl3{1});

    fprintf('>>> Select the folder for group "%s":\n', ripGroupLabels{g});
    folder3 = uigetdir(pwd, sprintf('Group folder: %s', ripGroupLabels{g}));
    if folder3 == 0, error('No folder selected. Aborted.'); end
    ripGroupFolders{g} = folder3;
end

fprintf('>>> Select the OUTPUT folder for ripple results:\n');
rippleFolder = uigetdir(ripGroupFolders{1}, 'OUTPUT — Part III Ripples');
if rippleFolder == 0, error('Aborted.'); end

%% ── P3.2  LIST SESSION SUBFOLDERS PER GROUP ──────────────────────────────

lfpPaths = cell(nGroupsR, 1);

for g = 1:nGroupsR
    folderInfo = dir(ripGroupFolders{g});
    subfolders = folderInfo([folderInfo.isdir]);
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
    lfpPaths{g} = {subfolders.name};
    fprintf('  Group "%s": %d session(s) found.\n', ripGroupLabels{g}, numel(lfpPaths{g}));
end

%% ── P3.3  RIPPLE DETECTION AND CLASSIFICATION (PER GROUP) ──────────────


allRipplesStruct = cell(nGroupsR, 1);
allCSD           = cell(nGroupsR, 1);
allCSDEvents     = cell(nGroupsR, 1);
allCSDSummary    = cell(nGroupsR, 1);
sinkSizes        = struct();

for g = 1:nGroupsR
    grpLabel = ripGroupLabels{g};
    paths    = lfpPaths{g};
    nSess    = numel(paths);

    fprintf('\n==== Processing group "%s": %d session(s) ====\n', grpLabel, nSess);

    ripplesCell = cell(nSess, 1);
    csdCell     = cell(nSess, 1);
    csdEvCell   = cell(nSess, 1);
    csdSumCell  = cell(nSess, 1);

    slmSize   = nan(nSess, 1);
    radSize   = nan(nSess, 1);
    otherSize = nan(nSess, 1);

    % Avoid "Graphics timeout" with invisible figures during batch detection
    set(0, 'DefaultFigureVisible', 'off');

    for ii = 1:nSess
        thisFolderName = fullfile(ripGroupFolders{g}, paths{ii});
        fprintf('  [%d/%d] %s\n', ii, nSess, paths{ii});
        cd(thisFolderName);

        % Load hippocampal channels
        layersFile = fullfile(thisFolderName, [paths{ii} '.hippocampalLayers.channelinfo.mat']);
        load(layersFile, 'hippocampalLayers');
        rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
        radCh         = hippocampalLayers.layers{1,1}.radiatum;
        SWChannel     = hippocampalLayers.layers{1,1}.radiatum;
        slmCh         = hippocampalLayers.layers{1,1}.slm;

        % Load theta intervals (generated in PASS 1 of Part I)
        thetaFile = fullfile(thisFolderName, [paths{ii} '.thetaEpochs6to12.states.mat']);
        load(thetaFile, 'thisthetaEpochs');
        thetaIntervals = thisthetaEpochs.thetaRun.newintervals;

        % Load LFP
        lfp = getLFP('all');

        try
            % CNN detection
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
                warning('computeRippleStats_new failed in session %s:\n  %s\nSession discarded.', ...
                    paths{ii}, ME_stats.message);
                fprintf('    Stack trace:\n');
                for sti = 1:numel(ME_stats.stack)
                    fprintf('      %s (line %d)\n', ME_stats.stack(sti).name, ME_stats.stack(sti).line);
                end
                ripplesCell{ii} = []; csdCell{ii} = []; csdEvCell{ii} = []; csdSumCell{ii} = [];
                close all; drawnow;
                continue;
            end

            % Fix valid_events: sync rippleStats.data to the number of peaks
            nRipClean = numel(ripples.peaks);
            nRipStats = size(ripples.rippleStats.data.peakAmplitude, 1);
            if nRipStats ~= nRipClean
                fprintf('    FIX valid_events: rippleStats.data had %d rows vs %d peaks. Trimming...\n', nRipStats, nRipClean);
                statFields = fieldnames(ripples.rippleStats.data);
                for sf = 1:numel(statFields)
                    fdata = ripples.rippleStats.data.(statFields{sf});
                    if isnumeric(fdata) && size(fdata,1) == nRipStats
                        ripples.rippleStats.data.(statFields{sf}) = fdata(1:nRipClean, :);
                    end
                end
            end

            % ── SAVE <rec>.ripples.events.mat ────────────────────────────
            %    Read by detectThetaEpochs (in P1R) and nonlocnonrip.
            %    Saved variable is named 'ripples' (buzcode convention).
            ripplesEvents.timestamps  = ripples.timestamps;
            ripplesEvents.peaks       = ripples.peaks;
            ripplesEvents.amplitude   = ripples.peakNormedPower;
            ripplesEvents.duration    = diff(ripples.timestamps, 1, 2);
            ripplesEvents.detectorinfo.detectorname = 'rippleMasterDetector';
            ripplesEvents.detectorinfo.detectionparms.passband   = RIPPLE_PASSBAND;
            ripplesEvents.detectorinfo.detectionparms.thresholds = RIPPLE_THRESHOLD;
            ripplesEventsFile = fullfile(thisFolderName, [paths{ii} '.ripples.events.mat']);
            S_rip.ripples = ripplesEvents;                 % helper struct -> variable 'ripples'
            save(ripplesEventsFile, '-struct', 'S_rip', '-mat');
            clear S_rip ripplesEvents;
            fprintf('    Saved: %s.ripples.events.mat\n', paths{ii});

            % Per-ripple CSD
            [csd, ~, csdEvents, ~] = buzcode_eacheventCSD(lfp, ripples.peaks, ...
                'plotCSD', false, 'plotLFP', false);

            % Summarize CSD per ripple (centered time average)
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
                    idxWin = true(size(csd.timestamps));  % fallback: full window
                end
                slice   = csdEvents(idxWin, :, jj);
                rowMean = mean(slice, 1);
                csdSum.csdSummary(jj, :) = rowMean;
            end

            % Categorize by sink (RAD / SLM / Other)
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
            warning('Error in session %s: %s', paths{ii}, ME.message);
            ripplesCell{ii} = [];
            csdCell{ii}     = [];
            csdEvCell{ii}   = [];
            csdSumCell{ii}  = [];
        end
        close all; drawnow;
    end

    % Restore figure visibility after batch detection
    set(0, 'DefaultFigureVisible', 'on');

    allRipplesStruct{g} = ripplesCell;
    allCSD{g}           = csdCell;
    allCSDEvents{g}     = csdEvCell;
    allCSDSummary{g}    = csdSumCell;

    grpField = matlab.lang.makeValidName(grpLabel);   % field names can't contain '-' etc.
    sinkSizes.(grpField).slm   = slmSize;
    sinkSizes.(grpField).rad   = radSize;
    sinkSizes.(grpField).other = otherSize;

    fprintf('  Group "%s" done.\n', grpLabel);
end

% Save full ripple structures
allRipplesFile = fullfile(rippleFolder, 'allrippledata.mat');
saveVars = cellfun(@(lbl) [matlab.lang.makeValidName(lbl) 'ripples'], ripGroupLabels, 'UniformOutput', false);  % variable names can't contain '-' etc.
for g = 1:nGroupsR
    eval(sprintf('%s = allRipplesStruct{%d};', saveVars{g}, g));
end
save(allRipplesFile, saveVars{:}, '-mat');
fprintf('\nRipple structures saved: %s\n', allRipplesFile);

  

%% ── P3.4  EXTRACT MEASURES AND ASSEMBLE TABLES ──────────────────────────

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

        % Extract measures
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

% Combine all groups
allRipples     = vertcat(allRipplesTables{:});
allRipplesRaw  = vertcat(allRipplesRaw{:});
allRipplesFilt = vertcat(allRipplesFilt{:});

% Save original table
writetable(allRipples, fullfile(rippleFolder, 'all_ripples_original.xlsx'));

% Filter outliers by PeakNormPow > PEAKNORM_CUTOFF
keepIdx = allRipples.PeakNormPow < PEAKNORM_CUTOFF;
allRipplesKeep = allRipples(keepIdx, :);

% Save per group
for g = 1:nGroupsR
    grpIdx = strcmp(allRipples.expcond, ripGroupLabels{g});
    grpKeepIdx = grpIdx & keepIdx;
    writetable(allRipples(grpIdx, :),     fullfile(rippleFolder, sprintf('%s_rippledata_original.xlsx', ripGroupLabels{g})));
    writetable(allRipples(grpKeepIdx, :), fullfile(rippleFolder, sprintf('%s_rippledata_keep.xlsx',     ripGroupLabels{g})));
end

% Save combined filtered table
writetable(allRipplesKeep, fullfile(rippleFolder, 'all_ripples_keep.xlsx'));

% Save waveforms
waveformFile = fullfile(rippleFolder, 'ripplewaveform.mat');
save(waveformFile, 'allRipplesRaw', 'allRipplesFilt', '-mat');

fprintf('\nRipple tables saved in: %s\n', rippleFolder);

%% ── P3.5  ADD SUBJECT COLUMN ─────────────────────────────────────────────
%  Animal ID = last chunk after the final '-' (getAnimalID).
%  e.g. 2026-06-09_16-47-11-F016 -> F016. - Esto es como lo tiene Paria

% if ~ismember('subject', allRipplesKeep.Properties.VariableNames)
%     allRipplesKeep.subject = arrayfun(@(s) string(getAnimalID(s)), allRipplesKeep.ID);
% end

%Así como lo guarda Alberto
if ~ismember('subject', allRipplesKeep.Properties.VariableNames)
    % Extrae ID de animal del nombre de sesion: "C089_190520_sess1" -> "C089"
    allRipplesKeep.subject = extractBefore(allRipplesKeep.ID, '_');
end

%% ── P3.6  OUTLIER REMOVAL (HIERARCHICAL ROBUST) ─────────────────────────

data = allRipplesKeep;
fprintf('\n=== Removing outliers (hierarchical robust MAD) ===\n');

missingVars = setdiff(OUTLIER_VARS, data.Properties.VariableNames);
if ~isempty(missingVars)
    error('Columns not found in table: %s\nAvailable columns: %s', ...
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

%% ── P3.7  HISTOGRAMS PER GROUP ───────────────────────────────────────────

vars = {'PeakFq','PeakAmp','spectralEntropy','Duration'};

if nGroupsR <= 3
    cols = [1.00 0.55 0.75; 0.65 0.85 1.00; 0.05 0.25 0.65];
    cols = cols(1:nGroupsR, :);
else
    cols = lines(nGroupsR);
end

alphaVal = 0.35;
nbins    = 30;

fprintf('\n=== Histogram variable diagnostics ===\n');
for i = 1:numel(vars)
    varname = vars{i};
    if ismember(varname, data.Properties.VariableNames)
        nvalid = sum(~isnan(data.(varname)) & isfinite(data.(varname)));
        fprintf('  %s: %d valid values\n', varname, nvalid);
    else
        fprintf('  %s: *** COLUMN NOT FOUND ***\n', varname);
        fprintf('    Available columns: %s\n', strjoin(data.Properties.VariableNames, ', '));
    end
end
fprintf('======================================\n\n');

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

fprintf('\n[Part III] Done. .ripples.events.mat files saved in each session folder.\n');
fprintf('           Analysis results in: %s\n\n', rippleFolder);


%% ========================================================================
%  PART I (REFRESH) — RECOMPUTE WITH NEW RIPPLES
%  Now detectThetaEpochs reads the freshly generated .ripples.events.mat.
%  Recomputes nonrun, rewrites .thetaEpochs6to12.states.mat (which ps1to120
%  reads), and saves Part I Excel + figures.
%% ========================================================================

%% ── P1R.1  REFRESH ───────────────────────────────────────────────────────

fprintf('\n##### REFRESH Part I (with new ripples) #####\n');
Tout = runThetaStage(groupFolders, groupLabels, P, figFolder, false);

%% ── P1R.2  SAVE EXCEL ─────────────────────────────────────────────────────

writetable(Tout.allrec,  fullfile(outFolder, 'theta_allrec.xlsx'));
writetable(Tout.run,     fullfile(outFolder, 'theta_run.xlsx'));
writetable(Tout.nonrun,  fullfile(outFolder, 'theta_nonrun.xlsx'));
writetable(Tout.loc,     fullfile(outFolder, 'locevent_table.xlsx'));
fprintf('\n[Part I] Excel files saved in: %s\n', outFolder);

%% ── P1R.3  GROUP FIGURES ──────────────────────────────────────────────────

saveGroupFig(Tout.run, Tout.nonrun, groupLabels, figFolder);
fprintf('[Part I] Done.\n\n');


%% ========================================================================
%  PART II — POWER SPECTRUM 1-120 Hz
%  (reads the .thetaEpochs6to12.states.mat already refreshed with new ripples)
%% ========================================================================

%% ── P2.0  BAND PARAMETERS ─────────────────────────────────────────────────

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

%% ── P2.1  CONFIGURATION: GROUPS AND OUTPUT FOLDER ────────────────────────

fprintf('=== PART II: Power Spectrum 1-120 Hz ===\n');

nGrpsStr2 = inputdlg('How many groups for the Power Spectrum?', 'PS — Groups', 1, {'2'});
if isempty(nGrpsStr2), error('Aborted.'); end
nGroupsPS = str2double(nGrpsStr2{1});
assert(~isnan(nGroupsPS) && nGroupsPS >= 1, 'Invalid number of groups.');

psGroupLabels = cell(nGroupsPS, 1);
psGroupData   = cell(nGroupsPS, 1);

for g = 1:nGroupsPS
    lbl2 = inputdlg(sprintf('Label for group %d (must match the one ps1to120 uses):', g), ...
                    sprintf('PS group %d', g), 1, {sprintf('G%d', g)});
    if isempty(lbl2), error('Aborted.'); end
    psGroupLabels{g} = strtrim(lbl2{1});

    fprintf('  Loading data for group "%s" — select its folder when the dialog appears...\n', psGroupLabels{g});
    psGroupData{g} = ps1to120_Paria(psGroupLabels{g});
    fprintf('  Group "%s" loaded: %d animals.\n', psGroupLabels{g}, height(psGroupData{g}.idtable));
end

fprintf('>>> Select the OUTPUT folder for Power Spectrum results:\n');
psOutFolder = uigetdir(pwd, 'OUTPUT — Part II Power Spectrum');
if psOutFolder == 0, error('Aborted.'); end

psFigFolder = fullfile(psOutFolder, 'figures', 'ps');
if ~exist(psFigFolder, 'dir'), mkdir(psFigFolder); end

freqs = psGroupData{1}.freqs;

%% ── P2.2  MAIN LOOP: BANDS × CONDITIONS ─────────────────────────────────
%  Keeps every recording (one row per recording, aligned to idtable) and
%  writes an ID column to the Excel files. Missing data stays as NaN (an
%  empty cell in Excel) instead of dropping the recording. Figures use the
%  Group + layers format (ID column removed just for plotting).

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
    fprintf('\n  Band: %s [%.1f-%.1f Hz] (%d bins)\n', ...
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

            % --- 1) Recording IDs (reference for the number of rows) ---
            ids   = getRecIDs(gData.idtable);
            nAnim = numel(ids);
            if nAnim == 0
                fprintf('    [WARN] Group "%s" has no records in idtable. Skipping.\n', gLabel);
                continue;
            end

            % --- 2) Per-layer matrices, aligned with ids ---
            %     If the layer or field does not exist (or is empty) -> all NaN
            layerMats = cell(nL, 1);
            for L = 1:nL
                lname = LAYER_NAMES{L};
                hasData = isfield(gData, lname) && isstruct(gData.(lname)) && ...
                          isfield(gData.(lname), condField) && ...
                          ~isempty(gData.(lname).(condField));
                if hasData
                    matFull = gData.(lname).(condField);
                else
                    fprintf('    [WARN] Group "%s": layer "%s" has no data for "%s". Filled with NaN.\n', ...
                        gLabel, lname, condField);
                    matFull = nan(nAnim, numel(freqs));
                end
                nr = size(matFull, 1);
                if nr ~= nAnim
                    fprintf(['    [WARN] Group "%s", layer "%s", "%s": %d rows vs %d records in idtable. ' ...
                             'Padded/trimmed at the end -> CHECK that the order matches.\n'], ...
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

            % --- 3) Metrics (NaNs are kept as gaps) ---
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

            % --- 4) Accumulate ---
            peakFqRows  = [peakFqRows;  peakFqBlock];
            powPeakRows = [powPeakRows; powPeakBlock];
            meanPowRows = [meanPowRows; meanPowBlock];
            groupCol    = [groupCol;    repmat({gLabel}, nAnim, 1)];
            idCol       = [idCol;       ids(:)];

        end

        if isempty(peakFqRows), continue; end

        yUnitPow = 'dB';
        if ~condIsDB, yUnitPow = 'uV^2/Hz'; end

        % Tables with ID (for Excel)
        idGrpTab = table(idCol, groupCol, 'VariableNames', ["ID","Group"]);
        tPeakFq  = [idGrpTab  array2table(peakFqRows,  'VariableNames', LAYER_NAMES)];
        tPowPeak = [idGrpTab  array2table(powPeakRows, 'VariableNames', LAYER_NAMES)];
        tMeanPow = [idGrpTab  array2table(meanPowRows, 'VariableNames', LAYER_NAMES)];

        % writetable leaves an empty cell for each numeric NaN
        writetable(tPeakFq,  fullfile(psOutFolder, sprintf('ps_peakfq_%s_%s.xlsx',    bandName, condLabel)));
        writetable(tPowPeak, fullfile(psOutFolder, sprintf('ps_powtopeak_%s_%s.xlsx', bandName, condLabel)));
        writetable(tMeanPow, fullfile(psOutFolder, sprintf('ps_meanpow_%s_%s.xlsx',   bandName, condLabel)));

        % Group summary (mean/SD with omitnan; n = all records in the group)
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

        % Figures: without the ID column (Group + layers format, as before)
        savePSFig(removevars(tPeakFq,  'ID'), sprintf('Peak Frequency — %s — %s', bandName, condLabel), ...
                  'Hz',     LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('peakfq_%s_%s',    bandName, condLabel));
        savePSFig(removevars(tPowPeak, 'ID'), sprintf('Power to Peak — %s — %s', bandName, condLabel), ...
                  yUnitPow, LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('powtopeak_%s_%s', bandName, condLabel));
        savePSFig(removevars(tMeanPow, 'ID'), sprintf('Mean Power — %s — %s', bandName, condLabel), ...
                  yUnitPow, LAYER_NAMES, psGroupLabels, psFigFolder, sprintf('meanpow_%s_%s',   bandName, condLabel));

    end   % conditions

    fprintf('    Band "%s" done.\n', bandName);
end   % bands

if ~isempty(groupSummary)
    writetable(groupSummary, fullfile(psOutFolder, 'ps_group_summary.xlsx'));
    fprintf('  Group summary saved: ps_group_summary.xlsx\n');
end

fprintf('\n[Part II] Done. Results in: %s\n', psOutFolder);
%% ========================================================================
%  PART II (b) — FULL POWER SPECTRA: export curves + group-average figures
%  Exports the whole 1-120 Hz spectrum per recording/layer (tidy long CSV:
%  one row per recording x layer x frequency) and draws the classic
%  group-average spectrum with shaded SEM, one panel per layer.
%  Uses psGroupData / freqs, already in memory after Part II.
%% ========================================================================

% dB spectra to export/plot: {field in psGroupData , nice label}
PS_FULL_CONDITIONS = {
    'allrec',              'AllRec';
    'locindb',             'Loc';
    'nonlocnonrippleindb', 'NonLocNonRip';
    'nonlocindb',          'NonLoc';
};

specOutFolder = fullfile(psOutFolder, 'full_spectra');
if ~exist(specOutFolder,'dir'), mkdir(specOutFolder); end
specFigFolder = fullfile(psFigFolder, 'spectra_curves');
if ~exist(specFigFolder,'dir'), mkdir(specFigFolder); end

fvec  = freqs(:).';                 % 1 x nFreq
nFreq = numel(fvec);
% group colours by label: XX = blues, XY = reds; Sham = dark, ORX = light
grpColors = zeros(nGroupsPS, 3);
fallbackC = lines(nGroupsPS);
for g = 1:nGroupsPS
    l = lower(psGroupLabels{g});
    isXX = contains(l,'xx'); isXY = contains(l,'xy');
    isSham = contains(l,'sham'); isORX = contains(l,'orx');
    if     isXX && isSham, grpColors(g,:) = [0.00 0.20 0.55];   % dark blue
    elseif isXX && isORX,  grpColors(g,:) = [0.45 0.70 0.95];   % light blue
    elseif isXY && isSham, grpColors(g,:) = [0.65 0.08 0.12];   % dark red
    elseif isXY && isORX,  grpColors(g,:) = [0.95 0.55 0.55];   % light red
    else,                  grpColors(g,:) = fallbackC(g,:);
    end
end
longBlocks = {};                    % accumulate tidy rows for export

for cc = 1:size(PS_FULL_CONDITIONS,1)
    cond = PS_FULL_CONDITIONS{cc,1};
    clab = PS_FULL_CONDITIONS{cc,2};

    fig = figure('Color','w','Position',[80 80 1200 820]);

    for L = 1:numel(LAYER_NAMES)
        lname = LAYER_NAMES{L};
        ax = subplot(2,2,L); hold(ax,'on');
        legH = gobjects(0); legTxt = {};

        for g = 1:nGroupsPS
            gData = psGroupData{g};
            gLab  = psGroupLabels{g};
            if ~isfield(gData,lname) || ~isstruct(gData.(lname)) || ...
               ~isfield(gData.(lname),cond) || isempty(gData.(lname).(cond))
                continue;
            end
            M   = gData.(lname).(cond);                 % [nRec x nFreq], dB
            ids = getRecIDs(gData.idtable);
            nR  = min(size(M,1), numel(ids));
            M   = M(1:nR,:);  ids = ids(1:nR);

            % --- tidy long rows for export ---
            FreqCol  = repmat(fvec(:), nR, 1);
            PowerCol = reshape(M.', [], 1);
            IDcol    = repelem(string(ids(:)), nFreq);
            blk = table(IDcol, repmat(string(gLab),nR*nFreq,1), ...
                        repmat(string(lname),nR*nFreq,1), ...
                        repmat(string(clab),nR*nFreq,1), FreqCol, PowerCol, ...
                        'VariableNames', ["ID","Group","Layer","Condition","Freq","Power_dB"]);
            longBlocks{end+1} = blk; %#ok<SAGROW>

            % --- group mean +/- SEM for the figure ---
            mu  = mean(M,1,'omitnan');
            nok = sum(~isnan(M),1);
            sem = std(M,0,1,'omitnan') ./ sqrt(max(nok,1));
            fin = isfinite(mu) & isfinite(sem);
            if ~any(fin), continue; end
            xx = [fvec(fin), fliplr(fvec(fin))];
            yy = [mu(fin)+sem(fin), fliplr(mu(fin)-sem(fin))];
            fill(ax, xx, yy, grpColors(g,:), 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
            h = plot(ax, fvec(fin), mu(fin), '-', 'LineWidth',1.5, 'Color', grpColors(g,:));
            legH(end+1) = h; legTxt{end+1} = gLab; %#ok<SAGROW>
        end

        title(ax, sprintf('%s — %s', clab, lname), 'Interpreter','none');
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power (dB)');
        xlim(ax,[fvec(1) fvec(end)]); grid(ax,'on'); box(ax,'off');
        if L==1 && ~isempty(legH), legend(ax, legH, legTxt, 'Location','best'); end
    end

    saveas(fig, fullfile(specFigFolder, sprintf('spectrum_%s.png', clab)));
    print(fig, fullfile(specFigFolder, sprintf('spectrum_%s.svg', clab)), '-dsvg');
    close(fig);
end

if ~isempty(longBlocks)
    fullSpectraLong = vertcat(longBlocks{:});
    writetable(fullSpectraLong, fullfile(specOutFolder, 'full_spectra_long.csv'));
    fprintf('\n[Part II-b] Full spectra (long CSV): %s\n', fullfile(specOutFolder,'full_spectra_long.csv'));
    fprintf('           Group-average spectrum figures in: %s\n', specFigFolder);
end
