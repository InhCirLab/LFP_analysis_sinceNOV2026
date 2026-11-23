function stats = analyzeCellCounts(T, groupColors, expConditions, groupLabels, outputDir)
% analyzeCellCounts  Compare total unit counts across experimental groups
%                    at both the session level and subject level.
%
% USAGE:
%   stats = analyzeCellCounts(T, groupColors, expConditions, groupLabels, outputDir)
%
% INPUTS:
%   T             - Table with at minimum columns: 'expcond', 'subject',
%                   'spikech' (or any column uniquely identifying a session;
%                   see SESSION IDENTIFICATION note below). Each row = one unit.
%   groupColors   - [nGroups x 3] matrix of RGB colors (0-1 scale), in the
%                   same order as expConditions.
%   expConditions - Cell array of condition codes to include, matching values
%                   in T.expcond. E.g. {'SS','OS','OD'}.
%                   Rows in T with expcond not in this list are ignored.
%   groupLabels   - Cell array of display labels, same order as expConditions.
%                   E.g. {'Ctrl','ORX','ORX+DHT'}.
%   outputDir     - String. Directory where figures and stats file are saved.
%                   Created if it does not exist.
%
% OUTPUT:
%   stats         - Struct with fields:
%                     .session  : per-session counts and test results
%                     .subject  : per-subject counts (sum + mean) and test results
%
% SESSION IDENTIFICATION:
%   A session is defined as a unique combination of (subject, spikech).
%   If your table uses a different column to identify sessions (e.g.,
%   'sessionID', 'date'), replace 'spikech' in the findgroups line below.
%
% STATISTICS:
%   - Kruskal-Wallis test across groups (non-parametric; count data,
%     typically non-normal, unequal N).
%   - Dunn's post-hoc pairwise comparisons with Bonferroni correction
%     (via multcompare on kruskalwallis output).
%   - Effect size: epsilon-squared (ε²) for Kruskal-Wallis
%     ε² = (H - k + 1) / (n - k), where H = KW statistic, k = nGroups,
%     n = total N. Interpreted as: small≥0.01, medium≥0.06, large≥0.14.
%   - Descriptives: median, IQR, mean, SD, N per group.
%
% FIGURES:
%   - Grouped box/scatter plots (session-level and subject-level panels).
%   - Saved as SVG and JPEG in outputDir.
%
% NOTES ON ANALYSIS LEVELS:
%   Session-level: treats each recording session as an independent
%     observation. Inflates N and violates independence (multiple sessions
%     per subject), but useful for visualizing raw data density and as a
%     supplementary analysis.
%   Subject-level (PRIMARY): aggregates per subject (sum of units across
%     sessions; mean units per session). This is the appropriate unit of
%     analysis for between-subjects designs. Both sum and mean are reported;
%     sum reflects total recorded yield, mean reflects per-session yield
%     controlling for unequal session counts.

%% ---- 0. Setup ----------------------------------------------------------

if nargin < 5 || isempty(outputDir)
    outputDir = fullfile(pwd, 'cellcount_output');
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

nGroups = numel(expConditions);
assert(numel(groupLabels) == nGroups, ...
    'groupLabels and expConditions must have the same number of elements.');
assert(size(groupColors,1) == nGroups, ...
    'groupColors must have one row per group.');

% Filter table to requested conditions only
keepRows = ismember(T.expcond, expConditions);
T = T(keepRows, :);

if height(T) == 0
    error('No rows remain after filtering for requested expConditions.');
end

% Force expcond to cell array of strings — guards against categorical arrays
% where strcmp against plain string inputs silently produces all zeros
if iscategorical(T.expcond) || ischar(T.expcond)
    T.expcond = cellstr(T.expcond);
end

%% ---- 1. Session-level counts -------------------------------------------
% A session = unique (subject, spikech) combination

[sessionGroups, sessionSubject, sessionChan] = findgroups(T.subject, T.spikech);
sessionCond  = splitapply(@(x) x(1), T.expcond,  sessionGroups);
sessionCount = splitapply(@numel,     T.expcond,  sessionGroups); % units per session

% Build session table
Tsess = table(sessionSubject, sessionChan, sessionCond, sessionCount, ...
    'VariableNames', {'subject','spikech','expcond','unitCount'});

%% ---- 2. Subject-level counts -------------------------------------------
% Aggregate per subject: sum and mean across sessions

[subjGroups, subjID] = findgroups(Tsess.subject);
subjCond     = splitapply(@(x) x(1),   Tsess.expcond,    subjGroups);
subjSum      = splitapply(@sum,         Tsess.unitCount,  subjGroups);
subjMean     = splitapply(@mean,        Tsess.unitCount,  subjGroups);
subjNSess    = splitapply(@numel,       Tsess.unitCount,  subjGroups);

Tsubj = table(subjID, subjCond, subjSum, subjMean, subjNSess, ...
    'VariableNames', {'subject','expcond','unitSum','unitMean','nSessions'});

%% ---- 3. Statistics -----------------------------------------------------

stats.session = runGroupStats(Tsess, 'unitCount', expConditions, groupLabels, 'Session-level');
stats.subject_sum  = runGroupStats(Tsubj, 'unitSum',  expConditions, groupLabels, 'Subject-level (sum)');
stats.subject_mean = runGroupStats(Tsubj, 'unitMean', expConditions, groupLabels, 'Subject-level (mean per session)');

%% ---- 4. Figures --------------------------------------------------------

makeFigure(Tsess,  'unitCount', expConditions, groupLabels, groupColors, ...
    'Session-level unit counts', outputDir, 'cellcounts_session');

makeFigure(Tsubj, 'unitSum', expConditions, groupLabels, groupColors, ...
    'Subject-level unit counts (sum)', outputDir, 'cellcounts_subject_sum');

makeFigure(Tsubj, 'unitMean', expConditions, groupLabels, groupColors, ...
    'Subject-level unit counts (mean per session)', outputDir, 'cellcounts_subject_mean');

%% ---- 5. Export stats to text -------------------------------------------

exportStats(stats, expConditions, groupLabels, outputDir);

fprintf('\nDone. Figures and stats saved to: %s\n', outputDir);

end % main function


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function s = runGroupStats(T, countVar, expConditions, groupLabels, label)
% Compute descriptives + Kruskal-Wallis + Dunn post-hoc for one data level.

nGroups = numel(expConditions);
s.label = label;
s.countVar = countVar;

% --- Descriptives per group
for g = 1:nGroups
    idx = strcmp(T.expcond, expConditions{g});
    vals = T.(countVar)(idx);
    s.groups(g).condition = expConditions{g};
    s.groups(g).label     = groupLabels{g};
    s.groups(g).n         = sum(idx);
    s.groups(g).median    = median(vals);
    s.groups(g).iqr       = iqr(vals);
    s.groups(g).mean      = mean(vals);
    s.groups(g).sd        = std(vals);
    s.groups(g).min       = min(vals);
    s.groups(g).max       = max(vals);
    s.groups(g).values    = vals;
end

% --- Kruskal-Wallis
allVals  = T.(countVar);
% Map expcond to integer group index in the order of expConditions
groupIdx = zeros(height(T), 1);
for g = 1:nGroups
    groupIdx(strcmp(T.expcond, expConditions{g})) = g;
end
assert(all(groupIdx > 0), ...
    'Some rows have expcond values not found in expConditions. Check table filtering.');

n = numel(allVals);
k = nGroups;

% Pass group labels as cell array of strings — required for multcompare to work
groupLabelsForKW = groupLabels(groupIdx)';   % groupLabels is a cell array: use () not {}
[s.kw.p, s.kw.tbl, s.kw.stats] = kruskalwallis(allVals, groupLabelsForKW, 'off');
H = s.kw.tbl{2,5}; % chi-sq statistic
s.kw.H = H;
s.kw.df = k - 1;

% Epsilon-squared effect size
s.kw.epsilon2 = (H - k + 1) / (n - k);
s.kw.epsilon2 = max(s.kw.epsilon2, 0); % floor at 0

% Effect size interpretation
if s.kw.epsilon2 < 0.01
    s.kw.epsilon2_interp = 'negligible';
elseif s.kw.epsilon2 < 0.06
    s.kw.epsilon2_interp = 'small';
elseif s.kw.epsilon2 < 0.14
    s.kw.epsilon2_interp = 'medium';
else
    s.kw.epsilon2_interp = 'large';
end

% --- Dunn post-hoc with Bonferroni correction
if k > 2
    [s.posthoc.c, ~, ~, s.posthoc.gnames] = ...
        multcompare(s.kw.stats, 'CType', 'bonferroni', 'Display', 'off');
    s.posthoc.pairs  = s.posthoc.c(:,1:2);
    s.posthoc.p_bonf = s.posthoc.c(:,6);

    nPairs = size(s.posthoc.c,1);
    s.posthoc.pairLabels = cell(nPairs,1);
    s.posthoc.r          = zeros(nPairs,1);
    for p = 1:nPairs
        gi1   = s.posthoc.c(p,1);
        gi2   = s.posthoc.c(p,2);
        name1 = s.posthoc.gnames{gi1};
        name2 = s.posthoc.gnames{gi2};
        s.posthoc.pairLabels{p} = sprintf('%s vs %s', name1, name2);
        idx1 = find(strcmp(groupLabels, name1), 1);
        idx2 = find(strcmp(groupLabels, name2), 1);
        if ~isempty(idx1) && ~isempty(idx2)
            s.posthoc.r(p) = rankBiserialR(s.groups(idx1).values, s.groups(idx2).values);
        else
            s.posthoc.r(p) = NaN;
        end
    end
else
    % Two groups: Mann-Whitney U instead
    v1 = s.groups(1).values;
    v2 = s.groups(2).values;
    [s.posthoc.p, ~, mwstat] = ranksum(v1, v2);
    s.posthoc.U = mwstat.ranksum;
    s.posthoc.r = rankBiserialR(v1, v2);
    s.posthoc.pairLabels = {sprintf('%s vs %s', groupLabels{1}, groupLabels{2})};
end

end % runGroupStats


function r = rankBiserialR(v1, v2)
% Rank-biserial correlation as effect size for Mann-Whitney U.
% r = 1 - (2U) / (n1*n2), ranges [-1, 1].
n1 = numel(v1); n2 = numel(v2);
if n1 == 0 || n2 == 0; r = NaN; return; end
U = 0;
for i = 1:n1
    U = U + sum(v1(i) > v2) + 0.5*sum(v1(i) == v2);
end
r = 1 - (2*U)/(n1*n2);
end


function makeFigure(T, countVar, expConditions, groupLabels, groupColors, titleStr, outputDir, fname)
% Boxplot with overlaid individual data points, one panel.

nGroups = numel(expConditions);
fig = figure('Visible','off','Color','w','Units','centimeters','Position',[0 0 10 12]);
ax = axes(fig); hold(ax,'on');

jitterAmt = 0.12;
allData = cell(nGroups,1);

for g = 1:nGroups
    idx  = strcmp(T.expcond, expConditions{g});
    vals = T.(countVar)(idx);
    allData{g} = vals;
    n = numel(vals);

    col = groupColors(g,:);

    % Box
    bh = boxplot(ax, vals, 'Positions', g, 'Widths', 0.4, ...
        'Colors', col, 'Symbol', '');
    set(bh, 'LineWidth', 1.5);

    % Jittered scatter
    xJit = g + (rand(n,1)-0.5)*2*jitterAmt;
    scatter(ax, xJit, vals, 28, col, 'filled', 'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeColor','none');
end

% Aesthetics
ax.XTick = 1:nGroups;
ax.XTickLabel = groupLabels;
ax.XLim = [0.4, nGroups+0.6];
ax.TickDir = 'out';
ax.Box = 'off';
ax.FontSize = 10;
ax.FontName = 'Arial';
ylabel(ax, strrep(countVar,'unit','Units'), 'FontSize',11);
title(ax, titleStr, 'FontSize',11,'FontWeight','normal');

% N labels below x-axis
for g = 1:nGroups
    text(ax, g, ax.YLim(1) - 0.05*range(ax.YLim), ...
        sprintf('n=%d', numel(allData{g})), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4]);
end

% Save
svgPath  = fullfile(outputDir, [fname '.svg']);
jpegPath = fullfile(outputDir, [fname '.jpg']);
print(fig, svgPath,  '-dsvg',  '-r0');
print(fig, jpegPath, '-djpeg', '-r300');
close(fig);

end % makeFigure


function exportStats(stats, expConditions, groupLabels, outputDir)
% Write a plain-text stats report suitable for methods/results reporting.

fpath = fullfile(outputDir, 'cellcounts_stats.txt');
fid   = fopen(fpath, 'w');

fprintf(fid, '================================================================\n');
fprintf(fid, '  UNIT COUNT ANALYSIS — Statistical Report\n');
fprintf(fid, '  Generated: %s\n', datestr(now));
fprintf(fid, '================================================================\n\n');

levels = {stats.session, stats.subject_sum, stats.subject_mean};

for L = 1:numel(levels)
    s = levels{L};
    fprintf(fid, '----------------------------------------------------------------\n');
    fprintf(fid, '  %s\n', s.label);
    fprintf(fid, '----------------------------------------------------------------\n\n');

    % Descriptives table
    fprintf(fid, '  Descriptives:\n');
    fprintf(fid, '  %-18s %5s %8s %8s %8s %8s %8s %8s\n', ...
        'Group','N','Median','IQR','Mean','SD','Min','Max');
    for g = 1:numel(s.groups)
        fprintf(fid, '  %-18s %5d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n', ...
            s.groups(g).label, ...
            s.groups(g).n, ...
            s.groups(g).median, ...
            s.groups(g).iqr, ...
            s.groups(g).mean, ...
            s.groups(g).sd, ...
            s.groups(g).min, ...
            s.groups(g).max);
    end

    % KW test
    fprintf(fid, '\n  Kruskal-Wallis test:\n');
    fprintf(fid, '    H(%d) = %.3f,  p = %.4f\n', s.kw.df, s.kw.H, s.kw.p);
    fprintf(fid, '    Effect size: epsilon² = %.3f (%s)\n', ...
        s.kw.epsilon2, s.kw.epsilon2_interp);

    % Post-hoc
    fprintf(fid, '\n  Post-hoc comparisons (Dunn / Bonferroni-corrected):\n');
    if isfield(s.posthoc, 'p_bonf')
        fprintf(fid, '  %-28s %10s %10s\n', 'Pair', 'p (adj)', 'r (RBC)');
        for p = 1:numel(s.posthoc.pairLabels)
            sig = '';
            pv  = s.posthoc.p_bonf(p);
            if pv < 0.001, sig = '***';
            elseif pv < 0.01, sig = '**';
            elseif pv < 0.05, sig = '*';
            end
            fprintf(fid, '  %-28s %10.4f %10.3f  %s\n', ...
                s.posthoc.pairLabels{p}, pv, s.posthoc.r(p), sig);
        end
    else
        fprintf(fid, '  %-28s %10s %10s\n', 'Pair', 'p', 'r (RBC)');
        fprintf(fid, '  %-28s %10.4f %10.3f\n', ...
            s.posthoc.pairLabels{1}, s.posthoc.p, s.posthoc.r);
    end

    fprintf(fid, '\n  Significance: * p<0.05, ** p<0.01, *** p<0.001\n');
    fprintf(fid, '  RBC = rank-biserial correlation (effect size for pairwise comparison)\n');
    fprintf(fid, '\n');
end

fprintf(fid, '================================================================\n');
fprintf(fid, '  NOTES\n');
fprintf(fid, '================================================================\n');
fprintf(fid, '  Session-level analysis treats each recording session as\n');
fprintf(fid, '  independent. This inflates N and violates independence\n');
fprintf(fid, '  assumptions. Treat as exploratory / supplementary only.\n\n');
fprintf(fid, '  Subject-level (sum): total units recorded per subject across\n');
fprintf(fid, '  all sessions. Primary recommended analysis.\n\n');
fprintf(fid, '  Subject-level (mean): mean units per session per subject.\n');
fprintf(fid, '  Controls for unequal session counts across subjects.\n\n');
fprintf(fid, '  Groups analyzed: ');
for g = 1:numel(expConditions)
    fprintf(fid, '%s (%s)', expConditions{g}, groupLabels{g});
    if g < numel(expConditions), fprintf(fid, ', '); end
end
fprintf(fid, '\n');
fclose(fid);

fprintf('Stats report saved to: %s\n', fpath);
end % exportStats