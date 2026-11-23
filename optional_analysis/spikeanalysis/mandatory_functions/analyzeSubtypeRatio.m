function stats = analyzeSubtypeRatio(T, groupColors, expConditions, groupLabels, ...
                                      inhSubtype, outputDir, nPerm)
% analyzeSubtypeRatio  Count inspection + ratio analysis for one INH subtype vs EXC.
%
% USAGE:
%   stats = analyzeSubtypeRatio(T, groupColors, expConditions, groupLabels, ...
%                               inhSubtype, outputDir, nPerm)
%
% INPUTS:
%   T             - Table with columns: 'expcond', 'subject', 'spikech', 'celltype'.
%                   Each row = one unit.
%   groupColors   - [nGroups x 3] RGB matrix (0-1), same order as expConditions.
%   expConditions - Cell array of condition codes to include (matches T.expcond).
%   groupLabels   - Cell array of display labels, same order as expConditions.
%   inhSubtype    - String. One of: 'PV+', 'SST+', 'ID2+', or 'all_inh'.
%                   'all_inh' reproduces the full INH/EXC analysis from
%                   analyzeCellTypeRatio for reference.
%   outputDir     - Output directory (created if absent).
%   nPerm         - Number of permutations (default: 10000).
%
% OUTPUT:
%   stats         - Struct with fields:
%                     .counts   : per-recording, per-subject, per-group count tables
%                     .ratio    : ratio analysis results (if yield is sufficient)
%                     .Tsubj    : subject-level aggregated table
%                     .Tsess    : session-level aggregated table
%                     .yieldsOK : logical — whether ratio stats were run
%
% WORKFLOW:
%   1. ALWAYS runs count inspection first (figures + tables + console).
%      Reports N of target subtype cells per recording, per subject, per group.
%      Flags subjects/groups with zero or very low counts (< 3 cells).
%   2. Only proceeds to ratio statistics if median subtype count per subject
%      is >= 2 in ALL groups. Otherwise stops after inspection with a warning.
%      You can override this by setting stats.yieldsOK = true manually and
%      calling runRatioStats(stats, ...) directly — but interpret with caution.
%
% CELL TYPE MAPPING:
%   Excitatory : CAMK2_DEEP, CAMK2_SUP
%   Inhibitory subtypes: PV+, SST+, ID2+
%
% STATISTICS (if yield sufficient):
%   PRIMARY  — Permutation test on logit-transformed fractions (Agresti-Coull
%              correction). Omnibus F + pairwise with Bonferroni correction.
%              Effect size: Cohen's d on logit scale.
%   SECONDARY — Kruskal-Wallis + Dunn/Bonferroni. Epsilon² effect size.
%
% FIGURES:
%   Count inspection:
%     - Per-recording counts: strip plot per group, coloured by subject
%     - Per-subject counts: bar + scatter per group
%     - Per-group heatmap: subjects x groups
%   Ratio analysis (if run):
%     - Fraction boxplot + scatter per group
%   All saved as SVG + JPEG.

%% ---- 0. Setup ----------------------------------------------------------

if nargin < 7 || isempty(nPerm);    nPerm = 10000; end
if nargin < 6 || isempty(outputDir)
    outputDir = fullfile(pwd, sprintf('subtype_%s_output', strrep(inhSubtype,'+','plus')));
end
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

nGroups = numel(expConditions);
assert(numel(groupLabels)  == nGroups, 'groupLabels must match expConditions length.');
assert(size(groupColors,1) == nGroups, 'groupColors must have one row per group.');

validSubtypes = {'PV+','SST+','ID2+','all_inh'};
assert(ismember(inhSubtype, validSubtypes), ...
    'inhSubtype must be one of: PV+, SST+, ID2+, all_inh.');

% Filter to requested conditions
keepRows = ismember(T.expcond, expConditions);
T = T(keepRows, :);
assert(height(T) > 0, 'No rows remain after filtering for requested expConditions.');

% Standardise string columns — must be cellstr for strcmp and {} indexing
if iscategorical(T.expcond)  || ischar(T.expcond);  T.expcond  = cellstr(T.expcond);  end
if iscategorical(T.celltype) || ischar(T.celltype);  T.celltype = cellstr(T.celltype); end
if iscategorical(T.subject)  || ischar(T.subject)
    T.subject = cellstr(T.subject);
elseif isnumeric(T.subject)
    T.subject = arrayfun(@(x) num2str(x), T.subject, 'UniformOutput', false);
end

fprintf('\n========================================================\n');
fprintf('  analyzeSubtypeRatio: %s vs EXC\n', inhSubtype);
fprintf('========================================================\n');

%% ---- 1. Map cell types -------------------------------------------------

excTypes = {'CAMK2_DEEP','CAMK2_SUP'};
if strcmp(inhSubtype, 'all_inh')
    targetTypes = {'PV+','SST+','ID2+'};
else
    targetTypes = {inhSubtype};
end

isEXC    = ismember(T.celltype, excTypes);
isTarget = ismember(T.celltype, targetTypes);
isOther  = ~isEXC & ~isTarget;

if any(isOther)
    uTypes = unique(T.celltype(isOther));
    fprintf('  Note: %d units with other celltypes excluded from ratio: %s\n', ...
        sum(isOther), strjoin(uTypes, ', '));
end

T.broadType = repmat({'OTHER'}, height(T), 1);
T.broadType(isEXC)    = {'EXC'};
T.broadType(isTarget) = {'TARGET'};

%% ---- 2. Aggregate to session and subject level -------------------------

[sessGroups, sessSubj, sessChan] = findgroups(T.subject, T.spikech);
sessCond    = splitapply(@(x) x(1),                         T.expcond,   sessGroups);
sessNTarget = splitapply(@(x) sum(strcmp(x,'TARGET')),      T.broadType, sessGroups);
sessNEXC    = splitapply(@(x) sum(strcmp(x,'EXC')),         T.broadType, sessGroups);
sessNTotal  = sessNTarget + sessNEXC;

Tsess = table(sessSubj, sessChan, sessCond, sessNTarget, sessNEXC, sessNTotal, ...
    'VariableNames', {'subject','spikech','expcond','nTarget','nEXC','nTotal'});

% Subject level: sum across sessions
[subjGroups, subjID] = findgroups(Tsess.subject);
subjCond    = splitapply(@(x) x(1), Tsess.expcond,   subjGroups);
subjNTarget = splitapply(@sum,      Tsess.nTarget,    subjGroups);
subjNEXC    = splitapply(@sum,      Tsess.nEXC,       subjGroups);
subjNSess   = splitapply(@numel,    Tsess.nTotal,     subjGroups);
subjNTotal  = subjNTarget + subjNEXC;

% Guard zero-total subjects
valid = subjNTotal > 0;
subjID      = subjID(valid);      subjCond    = subjCond(valid);
subjNTarget = subjNTarget(valid); subjNEXC    = subjNEXC(valid);
subjNTotal  = subjNTotal(valid);  subjNSess   = subjNSess(valid);

subjFrac  = subjNTarget ./ subjNTotal;
subjRatio = subjNTarget ./ max(subjNEXC, 1);

Tsubj = table(subjID, subjCond, subjNTarget, subjNEXC, subjNTotal, ...
              subjNSess, subjFrac, subjRatio, ...
    'VariableNames', {'subject','expcond','nTarget','nEXC','nTotal', ...
                      'nSessions','targetFraction','targetRatio'});

stats.Tsubj = Tsubj;
stats.Tsess = Tsess;

%% ---- 3. Count inspection -----------------------------------------------

stats.counts = inspectCounts(Tsess, Tsubj, expConditions, groupLabels, ...
    groupColors, inhSubtype, outputDir);

%% ---- 4. Yield check ----------------------------------------------------
% Check median subtype count per subject per group.
% Stats are run only if median >= 2 in ALL groups.

fprintf('\n  Yield check (subtype counts per subject):\n');
fprintf('  %-18s %8s %8s %8s %8s\n','Group','Median','Min','Max','N_subj');
minMedian = Inf;
for g = 1:nGroups
    idx  = strcmp(Tsubj.expcond, expConditions{g});
    vals = Tsubj.nTarget(idx);
    med  = median(vals);
    fprintf('  %-18s %8.1f %8d %8d %8d\n', ...
        groupLabels{g}, med, min(vals), max(vals), sum(idx));
    if med < minMedian; minMedian = med; end
end

YIELD_THRESHOLD = 2;   % minimum median subtype cells per subject per group
if minMedian < YIELD_THRESHOLD
    warning(['analyzeSubtypeRatio: Median %s count per subject is %.1f in at least ' ...
        'one group (threshold = %d).\nRatio statistics NOT run — insufficient yield.\n' ...
        'Descriptive counts are saved. Interpret figures with caution.'], ...
        inhSubtype, minMedian, YIELD_THRESHOLD);
    stats.yieldsOK = false;
    stats.ratio    = [];
    fprintf('\n  Ratio statistics skipped due to insufficient yield.\n');
    fprintf('  Output saved to: %s\n', outputDir);
    return
end

stats.yieldsOK = true;
fprintf('\n  Yield OK (median >= %d in all groups). Proceeding to ratio stats.\n', ...
    YIELD_THRESHOLD);

%% ---- 5. Ratio statistics -----------------------------------------------

fprintf('\n--- %s / EXC fraction (subject-level, PRIMARY) ---\n', inhSubtype);
stats.ratio.fraction = runGroupStats(Tsubj, 'targetFraction', ...
    expConditions, groupLabels, ...
    sprintf('Subject-level %s fraction', inhSubtype), nPerm);

fprintf('\n--- %s / EXC ratio (subject-level, supplementary) ---\n', inhSubtype);
stats.ratio.ratio_raw = runGroupStats(Tsubj, 'targetRatio', ...
    expConditions, groupLabels, ...
    sprintf('Subject-level %s/EXC ratio', inhSubtype), nPerm);

%% ---- 6. Ratio figure ---------------------------------------------------

makeFractionFigure(Tsubj, expConditions, groupLabels, groupColors, ...
    inhSubtype, outputDir);

%% ---- 7. Export stats ---------------------------------------------------

exportStats(stats, expConditions, groupLabels, inhSubtype, outputDir);

fprintf('\nDone. Output saved to: %s\n', outputDir);

end % main


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function counts = inspectCounts(Tsess, Tsubj, expConditions, groupLabels, ...
                                groupColors, inhSubtype, outputDir)
% Inspect, print, plot and export subtype cell counts at three levels:
% per-recording (session), per-subject, per-group summary.

nGroups = numel(expConditions);
counts.Tsess = Tsess;
counts.Tsubj = Tsubj;

%% --- Console table: per-recording ---------------------------------------
fprintf('\n  Per-recording counts (%s):\n', inhSubtype);
fprintf('  %-20s %-20s %-12s %8s %8s %8s\n', ...
    'Subject','SpikeCh','Group','nTarget','nEXC','nTotal');
for g = 1:nGroups
    idx = strcmp(Tsess.expcond, expConditions{g});
    rows = Tsess(idx,:);
    for r = 1:height(rows)
        flag = '';
        if rows.nTarget(r) == 0; flag = '  <- ZERO'; end
        % subject may be cellstr or numeric depending on input type
        if iscell(rows.subject)
            subjStr = rows.subject{r};
        else
            subjStr = num2str(rows.subject(r));
        end
        fprintf('  %-20s %-20s %-12s %8d %8d %8d%s\n', ...
            subjStr, num2str(rows.spikech(r)), groupLabels{g}, ...
            rows.nTarget(r), rows.nEXC(r), rows.nTotal(r), flag);
    end
end

%% --- Console table: per-subject -----------------------------------------
fprintf('\n  Per-subject counts (%s):\n', inhSubtype);
fprintf('  %-20s %-12s %8s %8s %8s %8s %8s\n', ...
    'Subject','Group','nSess','nTarget','nEXC','nTotal','Fraction');
for g = 1:nGroups
    idx = strcmp(Tsubj.expcond, expConditions{g});
    rows = Tsubj(idx,:);
    for r = 1:height(rows)
        flag = '';
        if rows.nTarget(r) == 0; flag = '  <- ZERO';
        elseif rows.nTarget(r) < 3; flag = '  <- LOW'; end
        if iscell(rows.subject)
            subjStr = rows.subject{r};
        else
            subjStr = num2str(rows.subject(r));
        end
        fprintf('  %-20s %-12s %8d %8d %8d %8d %8.3f%s\n', ...
            subjStr, groupLabels{g}, rows.nSessions(r), ...
            rows.nTarget(r), rows.nEXC(r), rows.nTotal(r), ...
            rows.targetFraction(r), flag);
    end
end

%% --- Export CSV ---------------------------------------------------------
% Per-recording
recFile = fullfile(outputDir, sprintf('counts_%s_per_recording.csv', ...
    strrep(inhSubtype,'+','plus')));
writetable(Tsess, recFile);

% Per-subject (add fraction)
subjFile = fullfile(outputDir, sprintf('counts_%s_per_subject.csv', ...
    strrep(inhSubtype,'+','plus')));
writetable(Tsubj, subjFile);
fprintf('\n  Count tables exported to:\n    %s\n    %s\n', recFile, subjFile);

%% --- Figure 1: Per-recording strip plot ---------------------------------
% Each dot = one recording session, colour = subject (within group)

fig1 = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[0 0 12 10]);
ax1 = axes(fig1); hold(ax1,'on');

for g = 1:nGroups
    idx   = strcmp(Tsess.expcond, expConditions{g});
    rows  = Tsess(idx,:);
    subjs = unique(rows.subject);
    nSubj = numel(subjs);
    % Shade each subject slightly differently within group color
    for s = 1:nSubj
        sidx  = strcmp(rows.subject, subjs{s});
        vals  = rows.nTarget(sidx);
        shade = 0.4 + 0.6*(s/nSubj);   % lightness gradient
        col   = groupColors(g,:) * shade + (1-shade)*[1 1 1];
        col   = min(max(col,0),1);
        nv    = numel(vals);
        xJit  = g + (rand(nv,1)-0.5)*0.25;
        scatter(ax1, xJit, vals, 35, col, 'filled', ...
            'MarkerEdgeColor', groupColors(g,:), 'LineWidth', 0.5);
    end
    % Group median line
    med = median(rows.nTarget);
    plot(ax1, [g-0.3 g+0.3], [med med], '-', 'Color', groupColors(g,:)*0.7, ...
        'LineWidth', 2);
end

% Zero-count warning line
ymax = max(Tsess.nTarget) + 1;
plot(ax1, [0.3 nGroups+0.7], [0.5 0.5], '--', 'Color',[0.8 0.2 0.2], ...
    'LineWidth', 0.8);

ax1.XTick      = 1:nGroups;
ax1.XTickLabel = groupLabels;
ax1.XLim       = [0.4 nGroups+0.6];
ax1.YLim       = [-0.5 ymax];
ax1.TickDir    = 'out';
ax1.Box        = 'off';
ax1.FontSize   = 10;
ax1.FontName   = 'Arial';
ylabel(ax1, sprintf('N %s cells', inhSubtype), 'FontSize', 11);
title(ax1,  sprintf('%s counts per recording', inhSubtype), ...
    'FontSize', 11, 'FontWeight', 'normal');
text(ax1, nGroups+0.5, 0.5, 'zero', 'FontSize', 7, 'Color', [0.8 0.2 0.2], ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','right');

saveFig(fig1, outputDir, sprintf('counts_%s_per_recording', ...
    strrep(inhSubtype,'+','plus')));

%% --- Figure 2: Per-subject bar + scatter --------------------------------

fig2 = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[0 0 12 10]);
ax2 = axes(fig2); hold(ax2,'on');

for g = 1:nGroups
    idx  = strcmp(Tsubj.expcond, expConditions{g});
    vals = Tsubj.nTarget(idx);
    col  = groupColors(g,:);
    n    = numel(vals);

    bh = boxplot(ax2, vals, 'Positions', g, 'Widths', 0.45, ...
        'Colors', col, 'Symbol', '');
    set(bh, 'LineWidth', 1.5);

    xJit = g + (rand(n,1)-0.5)*0.2;
    scatter(ax2, xJit, vals, 45, col, 'filled', ...
        'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','w','LineWidth',0.5);

    % Subject labels next to dots
    subjIdx = find(idx);
    for s = 1:n
        if iscell(Tsubj.subject)
            sLabel = Tsubj.subject{subjIdx(s)};
        else
            sLabel = num2str(Tsubj.subject(subjIdx(s)));
        end
        text(ax2, xJit(s)+0.07, vals(s), sLabel, ...
            'FontSize', 6, 'Color', [0.4 0.4 0.4], 'Interpreter','none');
    end
end

% Low-yield threshold line
plot(ax2, [0.3 nGroups+0.7], [2 2], '--', 'Color',[0.8 0.5 0.1], ...
    'LineWidth', 0.8);
text(ax2, nGroups+0.5, 2, 'threshold (n=2)', 'FontSize', 7, ...
    'Color',[0.8 0.5 0.1], 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right');

ax2.XTick      = 1:nGroups;
ax2.XTickLabel = groupLabels;
ax2.XLim       = [0.4 nGroups+0.6];
ax2.YLim       = [-0.5 max(Tsubj.nTarget)+2];
ax2.TickDir    = 'out';
ax2.Box        = 'off';
ax2.FontSize   = 10;
ax2.FontName   = 'Arial';
ylabel(ax2, sprintf('N %s cells (sum across sessions)', inhSubtype), 'FontSize',11);
title(ax2,  sprintf('%s counts per subject', inhSubtype), ...
    'FontSize', 11, 'FontWeight','normal');

% n subjects per group
for g = 1:nGroups
    idx = strcmp(Tsubj.expcond, expConditions{g});
    text(ax2, g, -0.4, sprintf('n=%d', sum(idx)), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4]);
end

saveFig(fig2, outputDir, sprintf('counts_%s_per_subject', ...
    strrep(inhSubtype,'+','plus')));

%% --- Figure 3: Heatmap — subjects x groups ------------------------------
% Rows = subjects (sorted by group), columns = [nTarget, nEXC, fraction]

allSubjIDs   = Tsubj.subject;
% ensure cellstr for y-tick labels
if ~iscell(allSubjIDs)
    allSubjIDs = arrayfun(@(x) num2str(x), allSubjIDs, 'UniformOutput', false);
end
allSubjConds = Tsubj.expcond;
nSubjTotal   = height(Tsubj);

% Sort by group order
sortOrder = zeros(nSubjTotal,1);
for g = 1:nGroups
    sortOrder(strcmp(allSubjConds, expConditions{g})) = g;
end
[~, sortIdx] = sort(sortOrder);
sortedSubj  = allSubjIDs(sortIdx);
sortedConds = allSubjConds(sortIdx);
sortedNTarg = Tsubj.nTarget(sortIdx);
sortedNEXC  = Tsubj.nEXC(sortIdx);
sortedFrac  = Tsubj.targetFraction(sortIdx);

heatData = [sortedNTarg, sortedNEXC, sortedFrac];

fig3 = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[0 0 14 max(8, nSubjTotal*0.7 + 3)]);
ax3 = axes(fig3);

imagesc(ax3, heatData);
colormap(ax3, flipud(gray(256)));   % built-in, light=low dark=high
cb = colorbar(ax3);
cb.FontSize = 8;

ax3.XTick      = 1:3;
ax3.XTickLabel = {sprintf('N %s', inhSubtype), 'N EXC', 'Fraction'};
ax3.YTick      = 1:nSubjTotal;
ax3.YTickLabel = sortedSubj;
ax3.TickDir    = 'out';
ax3.FontSize   = 9;
ax3.FontName   = 'Arial';
ax3.XAxisLocation = 'top';
title(ax3, sprintf('%s count heatmap (subjects)', inhSubtype), ...
    'FontSize', 11, 'FontWeight','normal');

% Cell text values
for r = 1:nSubjTotal
    for c = 1:3
        val = heatData(r,c);
        if c < 3
            txt = sprintf('%d', val);
        else
            txt = sprintf('%.2f', val);
        end
        text(ax3, c, r, txt, 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontSize', 8, 'Color','k');
    end
end

% Group separator lines
cumN = 0;
for g = 1:nGroups-1
    cumN = cumN + sum(strcmp(sortedConds, expConditions{g}));
    plot(ax3, [0.5 3.5], [cumN+0.5 cumN+0.5], 'w-', 'LineWidth', 1.5);
end

% Group labels on left margin
cumN = 0;
for g = 1:nGroups
    gn   = sum(strcmp(sortedConds, expConditions{g}));
    yMid = cumN + gn/2 + 0.5;
    text(ax3, 0.3, yMid, groupLabels{g}, 'FontSize', 8, ...
        'Color', groupColors(g,:), 'FontWeight','bold', ...
        'HorizontalAlignment','right', 'Rotation', 0);
    cumN = cumN + gn;
end

saveFig(fig3, outputDir, sprintf('counts_%s_heatmap', ...
    strrep(inhSubtype,'+','plus')));

end % inspectCounts


function s = runGroupStats(T, dv, expConditions, groupLabels, label, nPerm)

nGroups = numel(expConditions);
s.label = label;
s.dv    = dv;

% Descriptives
for g = 1:nGroups
    idx  = strcmp(T.expcond, expConditions{g});
    vals = T.(dv)(idx);
    s.groups(g).condition = expConditions{g};
    s.groups(g).label     = groupLabels{g};
    s.groups(g).n         = sum(idx);
    s.groups(g).values    = vals;
    s.groups(g).median    = median(vals);
    s.groups(g).iqr       = iqr(vals);
    s.groups(g).mean      = mean(vals);
    s.groups(g).sd        = std(vals);
    s.groups(g).min       = min(vals);
    s.groups(g).max       = max(vals);
end

% Group index vector
allVals  = T.(dv);
groupIdx = zeros(height(T),1);
for g = 1:nGroups
    groupIdx(strcmp(T.expcond, expConditions{g})) = g;
end
assert(all(groupIdx > 0), 'Some rows have expcond not in expConditions.');
n = numel(allVals); k = nGroups;

% Agresti-Coull logit on fraction DV using nTarget/nTotal columns
if ismember('nTarget', T.Properties.VariableNames) && ...
   ismember('nTotal',  T.Properties.VariableNames)
    p_adj = (T.nTarget + 0.5) ./ (T.nTotal + 1);
else
    p_adj = allVals;
    p_adj(p_adj <= 0) = 0.5 / (max(allVals)*10 + 1);
    p_adj(p_adj >= 1) = 1 - 0.5 / (max(allVals)*10 + 1);
end
logitVals = log(p_adj ./ (1 - p_adj));

% PRIMARY: permutation F
obsF = anovaF(logitVals, groupIdx, k);
rng(42);
permF = zeros(nPerm,1);
for i = 1:nPerm
    permF(i) = anovaF(logitVals, groupIdx(randperm(n)), k);
end
s.perm.F_obs  = obsF;
s.perm.p      = mean(permF >= obsF);
s.perm.nPerm  = nPerm;

% Pairwise permutation + Bonferroni
pairs  = nchoosek(1:k, 2);
nPairs = size(pairs,1);
s.perm.pairLabels = cell(nPairs,1);
s.perm.p_raw      = zeros(nPairs,1);
s.perm.p_bonf     = zeros(nPairs,1);
s.perm.d          = zeros(nPairs,1);

for p = 1:nPairs
    g1 = pairs(p,1); g2 = pairs(p,2);
    s.perm.pairLabels{p} = sprintf('%s vs %s', groupLabels{g1}, groupLabels{g2});
    v1 = logitVals(groupIdx==g1); v2 = logitVals(groupIdx==g2);
    combined = [v1;v2]; n1 = numel(v1); n2 = numel(v2);
    obsDiff = abs(mean(v1)-mean(v2));
    permDiff = zeros(nPerm,1);
    for i = 1:nPerm
        sh = combined(randperm(n1+n2));
        permDiff(i) = abs(mean(sh(1:n1))-mean(sh(n1+1:end)));
    end
    s.perm.p_raw(p) = mean(permDiff >= obsDiff);
    poolSD = sqrt(((n1-1)*var(v1)+(n2-1)*var(v2))/(n1+n2-2));
    s.perm.d(p) = (mean(v1)-mean(v2)) / max(poolSD, eps);
end
s.perm.p_bonf = min(s.perm.p_raw * nPairs, 1);

% SECONDARY: KW
[s.kw.p, s.kw.tbl, s.kw.stats] = kruskalwallis(allVals, groupLabels(groupIdx)', 'off');
H = s.kw.tbl{2,5};
s.kw.H  = H; s.kw.df = k-1;
s.kw.epsilon2 = max((H-k+1)/(n-k), 0);
if     s.kw.epsilon2 < 0.01, s.kw.epsilon2_interp = 'negligible';
elseif s.kw.epsilon2 < 0.06, s.kw.epsilon2_interp = 'small';
elseif s.kw.epsilon2 < 0.14, s.kw.epsilon2_interp = 'medium';
else,                         s.kw.epsilon2_interp = 'large';
end

if k > 2
    [s.kw.mc,~,~,s.kw.gnames] = multcompare(s.kw.stats,'CType','bonferroni','Display','off');
    s.kw.p_bonf = s.kw.mc(:,6);
    nPKW = size(s.kw.mc,1);
    s.kw.pairLabels = cell(nPKW,1);
    s.kw.r          = zeros(nPKW,1);
    for p = 1:nPKW
        n1 = s.kw.gnames{s.kw.mc(p,1)}; n2 = s.kw.gnames{s.kw.mc(p,2)};
        s.kw.pairLabels{p} = sprintf('%s vs %s', n1, n2);
        i1 = find(strcmp(groupLabels,n1),1); i2 = find(strcmp(groupLabels,n2),1);
        if ~isempty(i1)&&~isempty(i2)
            s.kw.r(p) = rankBiserialR(s.groups(i1).values, s.groups(i2).values);
        else; s.kw.r(p) = NaN; end
    end
else
    v1 = s.groups(1).values; v2 = s.groups(2).values;
    [s.kw.p_pair,~,mw] = ranksum(v1,v2);
    s.kw.U = mw.ranksum;
    s.kw.r = rankBiserialR(v1,v2);
    s.kw.pairLabels = {sprintf('%s vs %s', groupLabels{1}, groupLabels{2})};
end

end % runGroupStats


function F = anovaF(vals, groupIdx, k)
grandMean  = mean(vals); n = numel(vals);
SS_between = 0; SS_within = 0;
for g = 1:k
    v = vals(groupIdx==g);
    if isempty(v); continue; end
    SS_between = SS_between + numel(v)*(mean(v)-grandMean)^2;
    SS_within  = SS_within  + sum((v-mean(v)).^2);
end
df_b = k-1; df_w = n-k;
if df_w <= 0 || SS_within == 0; F = 0; return; end
F = (SS_between/df_b)/(SS_within/df_w);
end


function r = rankBiserialR(v1, v2)
n1 = numel(v1); n2 = numel(v2);
if n1==0||n2==0; r=NaN; return; end
U = 0;
for i = 1:n1
    U = U + sum(v1(i)>v2) + 0.5*sum(v1(i)==v2);
end
r = 1-(2*U)/(n1*n2);
end


function makeFractionFigure(Tsubj, expConditions, groupLabels, groupColors, ...
                             inhSubtype, outputDir)
nGroups   = numel(expConditions);
jitterAmt = 0.12;

fig = figure('Visible','off','Color','w','Units','centimeters','Position',[0 0 10 12]);
ax  = axes(fig); hold(ax,'on');

for g = 1:nGroups
    idx  = strcmp(Tsubj.expcond, expConditions{g});
    vals = Tsubj.targetFraction(idx);
    col  = groupColors(g,:); n = numel(vals);
    bh   = boxplot(ax, vals, 'Positions', g, 'Widths', 0.4, ...
        'Colors', col, 'Symbol','');
    set(bh,'LineWidth',1.5);
    xJit = g + (rand(n,1)-0.5)*2*jitterAmt;
    scatter(ax, xJit, vals, 40, col, 'filled', ...
        'MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
end

ax.XTick      = 1:nGroups;
ax.XTickLabel = groupLabels;
ax.XLim       = [0.4 nGroups+0.6];
ax.YLim       = [0 1];
ax.YTick      = 0:0.25:1;
ax.TickDir    = 'out'; ax.Box = 'off';
ax.FontSize   = 10;    ax.FontName = 'Arial';
ylabel(ax, sprintf('%s fraction  (N_{%s} / N_{total})', inhSubtype, inhSubtype), ...
    'FontSize',11);
title(ax, sprintf('%s fraction per subject', inhSubtype), ...
    'FontSize',11,'FontWeight','normal');

for g = 1:nGroups
    idx = strcmp(Tsubj.expcond, expConditions{g});
    text(ax, g, -0.06, sprintf('n=%d', sum(idx)), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4]);
end

saveFig(fig, outputDir, sprintf('ratio_%s_fraction', strrep(inhSubtype,'+','plus')));
end


function exportStats(stats, expConditions, groupLabels, inhSubtype, outputDir)

fname = sprintf('ratio_%s_stats.txt', strrep(inhSubtype,'+','plus'));
fpath = fullfile(outputDir, fname);
fid   = fopen(fpath,'w');

fprintf(fid,'================================================================\n');
fprintf(fid,'  %s / EXC RATIO ANALYSIS — Statistical Report\n', inhSubtype);
fprintf(fid,'  Generated: %s\n', datestr(now));
fprintf(fid,'================================================================\n\n');
fprintf(fid,'  Subtype analyzed : %s\n', inhSubtype);
fprintf(fid,'  Excitatory types : CAMK2_DEEP, CAMK2_SUP\n\n');

if ~stats.yieldsOK
    fprintf(fid,'  *** RATIO STATISTICS NOT RUN ***\n');
    fprintf(fid,'  Median subtype count per subject was below threshold (n=2)\n');
    fprintf(fid,'  in at least one group. See count tables for details.\n\n');
    fprintf(fid,'  Exported files:\n');
    fprintf(fid,'    counts_%s_per_recording.csv\n', strrep(inhSubtype,'+','plus'));
    fprintf(fid,'    counts_%s_per_subject.csv\n',   strrep(inhSubtype,'+','plus'));
    fclose(fid);
    fprintf('Stats report saved to: %s\n', fpath);
    return
end

levels = {stats.ratio.fraction, stats.ratio.ratio_raw};

for L = 1:numel(levels)
    s = levels{L};
    fprintf(fid,'----------------------------------------------------------------\n');
    fprintf(fid,'  %s\n', s.label);
    fprintf(fid,'----------------------------------------------------------------\n\n');

    fprintf(fid,'  Descriptives (raw values):\n');
    fprintf(fid,'  %-18s %5s %8s %8s %8s %8s %8s %8s\n', ...
        'Group','N','Median','IQR','Mean','SD','Min','Max');
    for g = 1:numel(s.groups)
        fprintf(fid,'  %-18s %5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            s.groups(g).label, s.groups(g).n, ...
            s.groups(g).median, s.groups(g).iqr, ...
            s.groups(g).mean, s.groups(g).sd, ...
            s.groups(g).min, s.groups(g).max);
    end

    fprintf(fid,'\n  PRIMARY — Permutation test on logit-transformed fractions:\n');
    fprintf(fid,'    (Agresti-Coull correction; n = %d permutations)\n', s.perm.nPerm);
    fprintf(fid,'    Omnibus: F_obs = %.4f,  p_perm = %.4f', s.perm.F_obs, s.perm.p);
    if s.perm.p<0.001; fprintf(fid,' ***');
    elseif s.perm.p<0.01; fprintf(fid,' **');
    elseif s.perm.p<0.05; fprintf(fid,' *'); end
    fprintf(fid,'\n\n');
    fprintf(fid,'  Pairwise (Bonferroni-corrected):\n');
    fprintf(fid,'  %-28s %10s %10s %10s\n','Pair','p_raw','p_bonf','d (logit)');
    for p = 1:numel(s.perm.pairLabels)
        pv = s.perm.p_bonf(p); sig = '';
        if pv<0.001; sig='***'; elseif pv<0.01; sig='**'; elseif pv<0.05; sig='*'; end
        fprintf(fid,'  %-28s %10.4f %10.4f %10.3f  %s\n', ...
            s.perm.pairLabels{p}, s.perm.p_raw(p), pv, s.perm.d(p), sig);
    end

    fprintf(fid,'\n  SECONDARY — Kruskal-Wallis on raw values:\n');
    fprintf(fid,'    H(%d) = %.3f,  p = %.4f\n', s.kw.df, s.kw.H, s.kw.p);
    fprintf(fid,'    Effect size: epsilon² = %.3f (%s)\n', ...
        s.kw.epsilon2, s.kw.epsilon2_interp);
    if isfield(s.kw,'p_bonf')
        fprintf(fid,'  Dunn post-hoc (Bonferroni):\n');
        fprintf(fid,'  %-28s %10s %10s\n','Pair','p (adj)','r (RBC)');
        for p = 1:numel(s.kw.pairLabels)
            pv = s.kw.p_bonf(p); sig = '';
            if pv<0.001; sig='***'; elseif pv<0.01; sig='**'; elseif pv<0.05; sig='*'; end
            fprintf(fid,'  %-28s %10.4f %10.3f  %s\n', ...
                s.kw.pairLabels{p}, pv, s.kw.r(p), sig);
        end
    else
        fprintf(fid,'  Mann-Whitney: p = %.4f, r = %.3f\n', s.kw.p_pair, s.kw.r);
    end
    fprintf(fid,'\n  Significance: * p<0.05  ** p<0.01  *** p<0.001\n');
    fprintf(fid,'  d = Cohen''s d on logit scale; RBC = rank-biserial correlation\n\n');
end

fprintf(fid,'================================================================\n');
fprintf(fid,'  NOTES\n');
fprintf(fid,'================================================================\n');
fprintf(fid,'  Analysis only runs if median subtype count >= 2 per subject\n');
fprintf(fid,'  in all groups. This threshold is conservative given small N;\n');
fprintf(fid,'  treat any significant result as exploratory.\n\n');
fprintf(fid,'  Fraction = N_%s / (N_%s + N_EXC) per subject (counts summed\n', ...
    inhSubtype, inhSubtype);
fprintf(fid,'  across sessions before computing fraction).\n\n');
fprintf(fid,'  Groups analyzed: ');
for g = 1:numel(expConditions)
    fprintf(fid,'%s (%s)', expConditions{g}, groupLabels{g});
    if g < numel(expConditions); fprintf(fid,', '); end
end
fprintf(fid,'\n');
fclose(fid);
fprintf('Stats report saved to: %s\n', fpath);
end % exportStats


function saveFig(fig, outputDir, fname)
print(fig, fullfile(outputDir,[fname '.svg']),  '-dsvg',  '-r0');
print(fig, fullfile(outputDir,[fname '.jpg']),  '-djpeg', '-r300');
close(fig);
end