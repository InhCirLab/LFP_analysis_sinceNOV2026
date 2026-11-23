function stats = analyzeCellTypeRatio(T, groupColors, expConditions, groupLabels, outputDir)
% analyzeCellTypeRatio  Compare inhibitory/excitatory cell proportions across groups.
%
% USAGE:
%   stats = analyzeCellTypeRatio(T, groupColors, expConditions, groupLabels, outputDir)
%
% INPUTS:
%   T             - Table with columns: 'expcond', 'subject', 'spikech',
%                   'celltype'. Each row = one unit.
%   groupColors   - [nGroups x 3] RGB matrix (0-1), same order as expConditions.
%   expConditions - Cell array of condition codes to include (matches T.expcond).
%   groupLabels   - Cell array of display labels, same order as expConditions.
%   outputDir     - Output directory (created if absent).
%
% OUTPUT:
%   stats   - Struct with subject-level results for:
%               .inh_fraction : INH/(INH+EXC) per subject — PRIMARY DV
%               .inh_pyr_ratio: N_INH/N_EXC per subject  — supplementary
%
% CELL TYPE MAPPING:
%   Excitatory : CAMK2_DEEP, CAMK2_SUP
%   Inhibitory : PV+, SST+, ID2+
%   Units with unrecognized celltype are excluded with a warning.
%
% ANALYSIS:
%   PRIMARY unit of analysis: subject (between-subjects design).
%
%   DV = INH fraction (N_INH / N_total per subject). Bounded [0,1].
%   Raw fraction reported in descriptives.
%
%   PRIMARY TEST — Permutation test on logit-transformed fractions.
%     Logit: log(p/(1-p)) with Agresti-Coull correction for boundary values:
%       p_adj = (N_INH + 0.5) / (N_total + 1)
%     Omnibus statistic: one-way ANOVA F on logit values.
%     Pairwise permutation tests, Bonferroni-corrected across pairs.
%     Effect size: Cohen's d on logit scale per pair.
%     n = 10,000 permutations.
%
%   SECONDARY TEST — Kruskal-Wallis on raw fractions + Dunn/Bonferroni.
%     Reported for completeness. epsilon² effect size.
%
% FIGURES:
%   Panel 1 — INH fraction per subject, boxplot + scatter, per group.
%   Panel 2 — Stacked bar: mean N_EXC and N_INH per subject per group.
%   Saved as SVG and JPEG.

%% ---- 0. Setup ----------------------------------------------------------

if nargin < 5 || isempty(outputDir)
    outputDir = fullfile(pwd, 'celltyperation_output');
end
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

nGroups = numel(expConditions);
assert(numel(groupLabels)    == nGroups, 'groupLabels must match expConditions length.');
assert(size(groupColors,1)   == nGroups, 'groupColors must have one row per group.');

% Filter to requested conditions
keepRows = ismember(T.expcond, expConditions);
T = T(keepRows, :);
assert(height(T) > 0, 'No rows remain after filtering for requested expConditions.');

% Standardise expcond to cellstr
if iscategorical(T.expcond) || ischar(T.expcond)
    T.expcond = cellstr(T.expcond);
end
if iscategorical(T.celltype) || ischar(T.celltype)
    T.celltype = cellstr(T.celltype);
end

%% ---- 1. Map celltype to EXC / INH -------------------------------------

excTypes = {'CAMK2_DEEP','CAMK2_SUP'};
inhTypes = {'PV+','SST+','ID2+'};

isEXC = ismember(T.celltype, excTypes);
isINH = ismember(T.celltype, inhTypes);
isUnknown = ~isEXC & ~isINH;

if any(isUnknown)
    uTypes = unique(T.celltype(isUnknown));
    warning('analyzeCellTypeRatio: %d units with unrecognized celltype excluded: %s', ...
        sum(isUnknown), strjoin(uTypes, ', '));
    T = T(~isUnknown, :);
    isEXC = ismember(T.celltype, excTypes);
    isINH = ismember(T.celltype, inhTypes);
end

T.broadType = repmat({'EXC'}, height(T), 1);
T.broadType(isINH) = {'INH'};

%% ---- 2. Subject-level aggregation ------------------------------------
% Per subject: count INH and EXC units, compute fraction and ratio.
% A session = unique (subject, spikech).

% First collapse to session level to get per-session cell type counts,
% then aggregate to subject. This makes sessions-per-subject transparent.

[sessGroups, sessSubj, sessChan] = findgroups(T.subject, T.spikech);
sessCond  = splitapply(@(x) x(1), T.expcond,   sessGroups);
sessNINH  = splitapply(@(x) sum(strcmp(x,'INH')), T.broadType, sessGroups);
sessNEXC  = splitapply(@(x) sum(strcmp(x,'EXC')), T.broadType, sessGroups);
sessNTot  = sessNINH + sessNEXC;

Tsess = table(sessSubj, sessChan, sessCond, sessNINH, sessNEXC, sessNTot, ...
    'VariableNames', {'subject','spikech','expcond','nINH','nEXC','nTotal'});

% Subject-level: sum counts across sessions, then compute fraction
[subjGroups, subjID] = findgroups(Tsess.subject);
subjCond  = splitapply(@(x) x(1),  Tsess.expcond,  subjGroups);
subjNINH  = splitapply(@sum,        Tsess.nINH,     subjGroups);
subjNEXC  = splitapply(@sum,        Tsess.nEXC,     subjGroups);
subjNSess = splitapply(@numel,      Tsess.nTotal,   subjGroups);
subjNTot  = subjNINH + subjNEXC;

% Guard: subjects with zero total units (shouldn't happen post-filter)
validSubj = subjNTot > 0;
if any(~validSubj)
    warning('%d subjects with zero classified units removed.', sum(~validSubj));
end
subjID    = subjID(validSubj);
subjCond  = subjCond(validSubj);
subjNINH  = subjNINH(validSubj);
subjNEXC  = subjNEXC(validSubj);
subjNTot  = subjNTot(validSubj);
subjNSess = subjNSess(validSubj);

subjFrac  = subjNINH ./ subjNTot;           % PRIMARY DV: INH fraction [0,1]
subjRatio = subjNINH ./ max(subjNEXC, 1);   % supplementary: INH/EXC ratio

Tsubj = table(subjID, subjCond, subjNINH, subjNEXC, subjNTot, ...
              subjNSess, subjFrac, subjRatio, ...
    'VariableNames', {'subject','expcond','nINH','nEXC','nTotal', ...
                      'nSessions','inhFraction','inhRatio'});

%% ---- 3. Statistics ----------------------------------------------------

fprintf('\n--- INH Fraction (subject-level, PRIMARY) ---\n');
stats.inh_fraction = runGroupStats(Tsubj, 'inhFraction', ...
    expConditions, groupLabels, 'Subject-level INH fraction');

fprintf('\n--- INH/EXC Ratio (subject-level, supplementary) ---\n');
stats.inh_ratio = runGroupStats(Tsubj, 'inhRatio', ...
    expConditions, groupLabels, 'Subject-level INH/EXC ratio');

% Also store the aggregated table for inspection
stats.Tsubj = Tsubj;
stats.Tsess = Tsess;

%% ---- 4. Figures -------------------------------------------------------

makeFractionFigure(Tsubj, expConditions, groupLabels, groupColors, outputDir);
makeStackedBarFigure(Tsubj, expConditions, groupLabels, groupColors, outputDir);

%% ---- 5. Export stats --------------------------------------------------

exportRatioStats(stats, expConditions, groupLabels, outputDir);

fprintf('\nDone. Output saved to: %s\n', outputDir);

end % main


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function s = runGroupStats(T, dv, expConditions, groupLabels, label, nPerm)
% Primary: permutation test on logit-transformed fractions.
% Secondary: Kruskal-Wallis on raw values.
% Both omnibus + pairwise reported.

if nargin < 6; nPerm = 10000; end

nGroups = numel(expConditions);
s.label = label;
s.dv    = dv;

% --- Descriptives (on raw fraction values)
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

% --- Build group index vector
allVals  = T.(dv);
groupIdx = zeros(height(T), 1);
for g = 1:nGroups
    groupIdx(strcmp(T.expcond, expConditions{g})) = g;
end
assert(all(groupIdx > 0), 'Some rows have expcond not in expConditions.');
n = numel(allVals);
k = nGroups;

% =========================================================
% PRIMARY: Permutation test on logit-transformed fractions
% =========================================================
% Agresti-Coull correction handles p=0 and p=1 boundary cases:
%   p_adj = (N_INH + 0.5) / (N_total + 1)
% Applied here to the DV values directly since DV = inhFraction = N_INH/N_total.
% Reconstruct corrected fraction from the stored counts if available,
% otherwise apply correction directly to the value.
if ismember('nINH', T.Properties.VariableNames) && ...
   ismember('nTotal', T.Properties.VariableNames)
    p_adj  = (T.nINH + 0.5) ./ (T.nTotal + 1);
else
    % fallback: shift boundary values by small epsilon
    p_adj = allVals;
    p_adj(p_adj == 0) = 0.5 / (max(allVals)*10 + 1);
    p_adj(p_adj == 1) = 1 - 0.5 / (max(allVals)*10 + 1);
end
logitVals = log(p_adj ./ (1 - p_adj));

% Omnibus statistic: one-way ANOVA F on logit values
obsF = anovaF(logitVals, groupIdx, k);

% Permutation distribution of F
rng(42);   % reproducibility — remove if you want different seeds each run
permF = zeros(nPerm, 1);
for i = 1:nPerm
    permF(i) = anovaF(logitVals, groupIdx(randperm(n)), k);
end
s.perm.F_obs   = obsF;
s.perm.p       = mean(permF >= obsF);
s.perm.nPerm   = nPerm;
s.perm.logitVals = logitVals;   % store for pairwise use

% Pairwise permutation tests (all pairs), Bonferroni-corrected
pairs    = nchoosek(1:k, 2);
nPairs   = size(pairs, 1);
s.perm.pairLabels = cell(nPairs, 1);
s.perm.p_raw      = zeros(nPairs, 1);
s.perm.p_bonf     = zeros(nPairs, 1);
s.perm.d          = zeros(nPairs, 1);   % Cohen's d on logit scale

for p = 1:nPairs
    g1 = pairs(p,1); g2 = pairs(p,2);
    s.perm.pairLabels{p} = sprintf('%s vs %s', groupLabels{g1}, groupLabels{g2});

    v1 = logitVals(groupIdx == g1);
    v2 = logitVals(groupIdx == g2);
    combined = [v1; v2];
    n1 = numel(v1); n2 = numel(v2);

    % Observed mean difference
    obsDiff = abs(mean(v1) - mean(v2));

    % Permutation distribution of |mean diff|
    permDiff = zeros(nPerm, 1);
    for i = 1:nPerm
        shuf = combined(randperm(n1+n2));
        permDiff(i) = abs(mean(shuf(1:n1)) - mean(shuf(n1+1:end)));
    end
    s.perm.p_raw(p) = mean(permDiff >= obsDiff);

    % Cohen's d on logit scale (pooled SD)
    pooledSD = sqrt(((n1-1)*var(v1) + (n2-1)*var(v2)) / (n1+n2-2));
    if pooledSD > 0
        s.perm.d(p) = (mean(v1) - mean(v2)) / pooledSD;
    else
        s.perm.d(p) = NaN;
    end
end
% Bonferroni correction
s.perm.p_bonf = min(s.perm.p_raw * nPairs, 1);

% =========================================================
% SECONDARY: Kruskal-Wallis on raw fractions
% =========================================================
groupLabelsForKW = groupLabels(groupIdx)';
[s.kw.p, s.kw.tbl, s.kw.stats] = kruskalwallis(allVals, groupLabelsForKW, 'off');
H = s.kw.tbl{2,5};
s.kw.H  = H;
s.kw.df = k - 1;

s.kw.epsilon2 = max((H - k + 1) / (n - k), 0);
if     s.kw.epsilon2 < 0.01, s.kw.epsilon2_interp = 'negligible';
elseif s.kw.epsilon2 < 0.06, s.kw.epsilon2_interp = 'small';
elseif s.kw.epsilon2 < 0.14, s.kw.epsilon2_interp = 'medium';
else,                         s.kw.epsilon2_interp = 'large';
end

if k > 2
    [s.kw.mc, ~, ~, s.kw.gnames] = ...
        multcompare(s.kw.stats, 'CType', 'bonferroni', 'Display', 'off');
    s.kw.p_bonf     = s.kw.mc(:,6);
    nPairsKW        = size(s.kw.mc,1);
    s.kw.pairLabels = cell(nPairsKW,1);
    s.kw.r          = zeros(nPairsKW,1);
    for p = 1:nPairsKW
        n1 = s.kw.gnames{s.kw.mc(p,1)};
        n2 = s.kw.gnames{s.kw.mc(p,2)};
        s.kw.pairLabels{p} = sprintf('%s vs %s', n1, n2);
        i1 = find(strcmp(groupLabels, n1), 1);
        i2 = find(strcmp(groupLabels, n2), 1);
        if ~isempty(i1) && ~isempty(i2)
            s.kw.r(p) = rankBiserialR(s.groups(i1).values, s.groups(i2).values);
        else
            s.kw.r(p) = NaN;
        end
    end
else
    v1 = s.groups(1).values; v2 = s.groups(2).values;
    [s.kw.p_pair, ~, mw] = ranksum(v1, v2);
    s.kw.U          = mw.ranksum;
    s.kw.r          = rankBiserialR(v1, v2);
    s.kw.pairLabels = {sprintf('%s vs %s', groupLabels{1}, groupLabels{2})};
end

end % runGroupStats


function F = anovaF(vals, groupIdx, k)
% One-way ANOVA F statistic. vals and groupIdx are column vectors.
grandMean = mean(vals);
n = numel(vals);
SS_between = 0;
SS_within  = 0;
for g = 1:k
    v = vals(groupIdx == g);
    if isempty(v); continue; end
    SS_between = SS_between + numel(v) * (mean(v) - grandMean)^2;
    SS_within  = SS_within  + sum((v - mean(v)).^2);
end
df_between = k - 1;
df_within  = n - k;
if df_within <= 0 || SS_within == 0
    F = 0; return;
end
F = (SS_between / df_between) / (SS_within / df_within);
end


function r = rankBiserialR(v1, v2)
n1 = numel(v1); n2 = numel(v2);
if n1 == 0 || n2 == 0; r = NaN; return; end
U = 0;
for i = 1:n1
    U = U + sum(v1(i) > v2) + 0.5*sum(v1(i) == v2);
end
r = 1 - (2*U)/(n1*n2);
end


function makeFractionFigure(Tsubj, expConditions, groupLabels, groupColors, outputDir)
% Box + scatter of INH fraction per subject, one box per group.

nGroups   = numel(expConditions);
jitterAmt = 0.12;

fig = figure('Visible','off','Color','w','Units','centimeters','Position',[0 0 10 12]);
ax  = axes(fig); hold(ax,'on');

for g = 1:nGroups
    idx  = strcmp(Tsubj.expcond, expConditions{g});
    vals = Tsubj.inhFraction(idx);
    col  = groupColors(g,:);
    n    = numel(vals);

    bh = boxplot(ax, vals, 'Positions', g, 'Widths', 0.4, ...
        'Colors', col, 'Symbol', '');
    set(bh, 'LineWidth', 1.5);

    xJit = g + (rand(n,1)-0.5)*2*jitterAmt;
    scatter(ax, xJit, vals, 40, col, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
end

ax.XTick      = 1:nGroups;
ax.XTickLabel = groupLabels;
ax.XLim       = [0.4, nGroups+0.6];
ax.YLim       = [0, 1];
ax.YTick      = 0:0.25:1;
ax.TickDir    = 'out';
ax.Box        = 'off';
ax.FontSize   = 10;
ax.FontName   = 'Arial';
ylabel(ax, 'INH fraction  (N_{INH} / N_{total})', 'FontSize', 11);
title(ax,  'Inhibitory cell fraction per subject', 'FontSize', 11, 'FontWeight', 'normal');

% n labels
yBase = -0.06;
for g = 1:nGroups
    idx = strcmp(Tsubj.expcond, expConditions{g});
    text(ax, g, yBase, sprintf('n=%d', sum(idx)), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4]);
end

saveFig(fig, outputDir, 'celltyperation_inh_fraction');
end % makeFractionFigure


function makeStackedBarFigure(Tsubj, expConditions, groupLabels, groupColors, outputDir)
% Stacked bar: mean N_EXC (bottom) and N_INH (top) per group.
% Individual subject points overlaid on INH counts for transparency.

nGroups = numel(expConditions);
meanEXC = zeros(nGroups,1);
meanINH = zeros(nGroups,1);

for g = 1:nGroups
    idx = strcmp(Tsubj.expcond, expConditions{g});
    meanEXC(g) = mean(Tsubj.nEXC(idx));
    meanINH(g) = mean(Tsubj.nINH(idx));
end

fig = figure('Visible','off','Color','w','Units','centimeters','Position',[0 0 10 12]);
ax  = axes(fig); hold(ax,'on');

excColor = [0.75 0.75 0.75];   % neutral grey for EXC
jitterAmt = 0.12;

bh = bar(ax, 1:nGroups, [meanEXC, meanINH], 'stacked', 'BarWidth', 0.5);
bh(1).FaceColor = excColor;
bh(1).EdgeColor = 'none';

% Color each INH bar segment by group color
for g = 1:nGroups
    bh(2).FaceColor = 'flat';
    bh(2).CData(g,:) = groupColors(g,:);
end
bh(2).EdgeColor = 'none';

% Overlay individual INH counts (stacked on top of EXC bar)
for g = 1:nGroups
    idx  = strcmp(Tsubj.expcond, expConditions{g});
    nINH = Tsubj.nINH(idx);
    nEXC = Tsubj.nEXC(idx);
    yPos = nEXC + nINH;   % top of stacked bar for each subject
    xJit = g + (rand(numel(yPos),1)-0.5)*2*jitterAmt;
    scatter(ax, xJit, yPos, 35, groupColors(g,:), 'filled', ...
        'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
end

ax.XTick      = 1:nGroups;
ax.XTickLabel = groupLabels;
ax.XLim       = [0.4, nGroups+0.6];
ax.TickDir    = 'out';
ax.Box        = 'off';
ax.FontSize   = 10;
ax.FontName   = 'Arial';
ylabel(ax, 'Mean units per subject', 'FontSize', 11);
title(ax,  'EXC / INH unit counts per group', 'FontSize', 11, 'FontWeight', 'normal');
legend(ax, [bh(1), bh(2)], {'EXC','INH'}, 'Location','northwest', ...
    'Box','off','FontSize',9);

saveFig(fig, outputDir, 'celltyperation_stacked_counts');
end % makeStackedBarFigure


function saveFig(fig, outputDir, fname)
print(fig, fullfile(outputDir, [fname '.svg']),  '-dsvg',  '-r0');
print(fig, fullfile(outputDir, [fname '.jpg']),  '-djpeg', '-r300');
close(fig);
end


function exportRatioStats(stats, expConditions, groupLabels, outputDir)

fpath = fullfile(outputDir, 'celltyperation_stats.txt');
fid   = fopen(fpath, 'w');

fprintf(fid, '================================================================\n');
fprintf(fid, '  INH/EXC RATIO ANALYSIS — Statistical Report\n');
fprintf(fid, '  Generated: %s\n', datestr(now));
fprintf(fid, '================================================================\n\n');
fprintf(fid, '  Cell type mapping:\n');
fprintf(fid, '    Excitatory : CAMK2_DEEP, CAMK2_SUP\n');
fprintf(fid, '    Inhibitory : PV+, SST+, ID2+\n\n');

levels = {stats.inh_fraction, stats.inh_ratio};

for L = 1:numel(levels)
    s = levels{L};
    fprintf(fid, '----------------------------------------------------------------\n');
    fprintf(fid, '  %s\n', s.label);
    fprintf(fid, '----------------------------------------------------------------\n\n');

    % Descriptives
    fprintf(fid, '  Descriptives (raw values):\n');
    fprintf(fid, '  %-18s %5s %8s %8s %8s %8s %8s %8s\n', ...
        'Group','N','Median','IQR','Mean','SD','Min','Max');
    for g = 1:numel(s.groups)
        fprintf(fid, '  %-18s %5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            s.groups(g).label, s.groups(g).n, ...
            s.groups(g).median, s.groups(g).iqr, ...
            s.groups(g).mean,   s.groups(g).sd, ...
            s.groups(g).min,    s.groups(g).max);
    end

    % PRIMARY: Permutation test
    fprintf(fid, '\n  PRIMARY — Permutation test on logit-transformed fractions:\n');
    fprintf(fid, '    (Agresti-Coull correction applied; n = %d permutations)\n', ...
        s.perm.nPerm);
    fprintf(fid, '    Omnibus: F_obs = %.4f,  p_perm = %.4f', ...
        s.perm.F_obs, s.perm.p);
    if     s.perm.p < 0.001, fprintf(fid, ' ***');
    elseif s.perm.p < 0.01,  fprintf(fid, ' **');
    elseif s.perm.p < 0.05,  fprintf(fid, ' *');
    end
    fprintf(fid, '\n\n');

    fprintf(fid, '  Pairwise permutation tests (Bonferroni-corrected):\n');
    fprintf(fid, '  %-28s %10s %10s %10s\n', 'Pair','p_raw','p_bonf','d (logit)');
    for p = 1:numel(s.perm.pairLabels)
        pv  = s.perm.p_bonf(p);
        sig = '';
        if pv < 0.001, sig = '***'; elseif pv < 0.01, sig = '**'; elseif pv < 0.05, sig = '*'; end
        fprintf(fid, '  %-28s %10.4f %10.4f %10.3f  %s\n', ...
            s.perm.pairLabels{p}, s.perm.p_raw(p), pv, s.perm.d(p), sig);
    end

    % SECONDARY: KW
    fprintf(fid, '\n  SECONDARY — Kruskal-Wallis on raw values:\n');
    fprintf(fid, '    H(%d) = %.3f,  p = %.4f\n', s.kw.df, s.kw.H, s.kw.p);
    fprintf(fid, '    Effect size: epsilon² = %.3f (%s)\n', ...
        s.kw.epsilon2, s.kw.epsilon2_interp);

    if isfield(s.kw, 'p_bonf')
        fprintf(fid, '  Dunn post-hoc (Bonferroni):\n');
        fprintf(fid, '  %-28s %10s %10s\n', 'Pair','p (adj)','r (RBC)');
        for p = 1:numel(s.kw.pairLabels)
            pv  = s.kw.p_bonf(p);
            sig = '';
            if pv < 0.001, sig = '***'; elseif pv < 0.01, sig = '**'; elseif pv < 0.05, sig = '*'; end
            fprintf(fid, '  %-28s %10.4f %10.3f  %s\n', ...
                s.kw.pairLabels{p}, pv, s.kw.r(p), sig);
        end
    else
        fprintf(fid, '  Mann-Whitney: p = %.4f,  r (RBC) = %.3f\n', ...
            s.kw.p_pair, s.kw.r);
    end

    fprintf(fid, '\n  Significance: * p<0.05  ** p<0.01  *** p<0.001\n');
    fprintf(fid, '  d = Cohen''s d on logit scale; RBC = rank-biserial correlation\n\n');
end

fprintf(fid, '================================================================\n');
fprintf(fid, '  NOTES\n');
fprintf(fid, '================================================================\n');
fprintf(fid, '  PRIMARY DV: INH fraction = N_INH / (N_INH + N_EXC) per subject.\n');
fprintf(fid, '  Counts summed across sessions per subject before computing fraction.\n\n');
fprintf(fid, '  RATIONALE FOR PERMUTATION TEST:\n');
fprintf(fid, '  INH fraction is bounded [0,1] with non-constant variance.\n');
fprintf(fid, '  With n=4-5 subjects/group, asymptotic tests (KW, beta regression)\n');
fprintf(fid, '  cannot rely on large-sample approximations. Permutation testing\n');
fprintf(fid, '  is exact under the null and requires no distributional assumptions.\n');
fprintf(fid, '  Logit transform stabilises variance before computing F statistic.\n\n');
fprintf(fid, '  KW reported as secondary for consistency with convention.\n\n');
fprintf(fid, '  Groups analyzed: ');
for g = 1:numel(expConditions)
    fprintf(fid, '%s (%s)', expConditions{g}, groupLabels{g});
    if g < numel(expConditions), fprintf(fid, ', '); end
end
fprintf(fid, '\n');
fclose(fid);

fprintf('Stats report saved to: %s\n', fpath);
end % exportRatioStats