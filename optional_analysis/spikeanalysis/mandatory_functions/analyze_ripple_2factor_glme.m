function results = analyze_ripple_2factor_glme(data, measure_name, distribution, ...
    factor1_name, factor1_levels, ...
    factor2_name, factor2_levels, ...
    cell_colors)
% ANALYZE_RIPPLE_2FACTOR_GLME
% Two-factor generalized linear mixed-effects models (main effects + interaction)
% for gamma or weibull distributions, with nested structure:
%   ripples within recordings (ID) within subjects
%
% Contexts:
%   Context A: ID nested in subject  --> y ~ F1*F2 + (1|subject) + (1|subject:ID)
%   Context B: ID independent        --> y ~ F1*F2 + (1|ID)
%
% Results text, figure (JPEG + SVG), and .mat are auto-saved to pwd.
%
% INPUT:
%   data           - table with columns: ID, subject, <factor1_name>, <factor2_name>, <measure_name>
%   measure_name   - string, name of the measure column
%   distribution   - 'gamma' or 'weibull'
%
%   factor1_name   - string, name of first factor column (e.g. 'expcond')
%   factor1_levels - cell array of level names IN ORDER (e.g. {'XX','XY'})
%
%   factor2_name   - string, name of second factor column (e.g. 'sex')
%   factor2_levels - cell array of level names IN ORDER (e.g. {'male','female'})
%
%   cell_colors    - (nF2*nF1) x 3 RGB matrix, one row per cell, in plot order
%                    (factor2 outer/left-right groups, factor1 inner/within group):
%                      Row 1: factor2_levels{1} + factor1_levels{1}
%                      Row 2: factor2_levels{1} + factor1_levels{2}
%                      Row 3: factor2_levels{2} + factor1_levels{1}
%                      Row 4: factor2_levels{2} + factor1_levels{2}  ... etc.
%                    If omitted or empty, lines() colormap is used.
%
% OUTPUT:
%   results - structure with contextA, contextB, and descriptive fields
%
% FILES SAVED (all in pwd):
%   GLME2way_<dist>_<measure>_<f1>_<f2>_results.txt
%   GLME2way_<dist>_<measure>_<f1>_<f2>_contextA.jpeg / .svg
%   GLME2way_<dist>_<measure>_<f1>_<f2>_results.mat
%
% EXAMPLE:
%   cell_colors = [
%       0.30 0.50 0.90;   % male   + XX
%       0.10 0.20 0.70;   % male   + XY
%       0.90 0.40 0.60;   % female + XX
%       0.70 0.10 0.40;   % female + XY
%   ];
%   results = analyze_ripple_2factor_glme(data, 'peak_freq', 'gamma', ...
%       'expcond', {'XX','XY'}, ...
%       'sex',     {'male','female'}, ...
%       cell_colors);

    %% ----------------------------------------------------------------
    %  Validate inputs
    %% ----------------------------------------------------------------
    if ~ismember(lower(distribution), {'gamma', 'weibull'})
        error('Distribution must be ''gamma'' or ''weibull''');
    end
    distribution = lower(distribution);
    dist_name    = 'Gamma';   % fitglme uses Gamma + log link for both gamma and weibull
    link         = 'log';

    nF1 = numel(factor1_levels);
    nF2 = numel(factor2_levels);

    if nargin < 8 || isempty(cell_colors)
        cell_colors = lines(nF1 * nF2);
    else
        assert(size(cell_colors,1) == nF1*nF2 && size(cell_colors,2) == 3, ...
            'cell_colors must be a %dx3 RGB matrix (%d rows = nF2 x nF1).', nF1*nF2, nF1*nF2);
    end

    %% ----------------------------------------------------------------
    %  Extract and clean data (y > 0 required for gamma/weibull)
    %% ----------------------------------------------------------------
    y       = data.(measure_name);
    f1_raw  = data.(factor1_name);
    f2_raw  = data.(factor2_name);
    subject = categorical(data.subject);
    ID      = categorical(data.ID);

    f1 = categorical(f1_raw, factor1_levels);
    f2 = categorical(f2_raw, factor2_levels);

    valid   = ~isnan(y) & isfinite(y) & y > 0 & ~isundefined(f1) & ~isundefined(f2);
    y       = y(valid);
    f1      = f1(valid);
    f2      = f2(valid);
    subject = subject(valid);
    ID      = ID(valid);

    tbl = table(y, f1, f2, subject, ID);

    %% ----------------------------------------------------------------
    %  Build omnibus contrast matrix for the interaction model
    %  Fixed effects order: intercept, f1 dummies, f2 dummies, interaction dummies
    %  Omnibus tests each factor and interaction separately via coefTest
    %% ----------------------------------------------------------------

    %% ----------------------------------------------------------------
    %  CONTEXT A: ID nested within subject
    %% ----------------------------------------------------------------
    formulaA = 'y ~ f1*f2 + (1|subject) + (1|subject:ID)';
    modelA   = fitglme(tbl, formulaA, 'Distribution', dist_name, 'Link', link);

    [pA_f1, pA_f2, pA_int] = omnibus_pvalues(modelA, nF1, nF2);

    results.contextA.model         = modelA;
    results.contextA.AIC           = modelA.ModelCriterion.AIC;
    results.contextA.BIC           = modelA.ModelCriterion.BIC;
    results.contextA.dispersion    = modelA.Dispersion;
    results.contextA.p_f1          = pA_f1;
    results.contextA.p_f2          = pA_f2;
    results.contextA.p_interaction = pA_int;

    %% ----------------------------------------------------------------
    %  CONTEXT B: ID independent
    %% ----------------------------------------------------------------
    formulaB = 'y ~ f1*f2 + (1|ID)';
    modelB   = fitglme(tbl, formulaB, 'Distribution', dist_name, 'Link', link);

    [pB_f1, pB_f2, pB_int] = omnibus_pvalues(modelB, nF1, nF2);

    results.contextB.model         = modelB;
    results.contextB.AIC           = modelB.ModelCriterion.AIC;
    results.contextB.BIC           = modelB.ModelCriterion.BIC;
    results.contextB.dispersion    = modelB.Dispersion;
    results.contextB.p_f1          = pB_f1;
    results.contextB.p_f2          = pB_f2;
    results.contextB.p_interaction = pB_int;

    %% ----------------------------------------------------------------
    %  Descriptive statistics (per cell)
    %% ----------------------------------------------------------------
    for i = 1:nF1
        for j = 1:nF2
            idx = (f1 == factor1_levels{i}) & (f2 == factor2_levels{j});
            gd  = y(idx);
            tag = matlab.lang.makeValidName([factor1_levels{i} '_' factor2_levels{j}]);
            results.descriptive.(tag).mean         = mean(gd);
            results.descriptive.(tag).std          = std(gd);
            results.descriptive.(tag).median       = median(gd);
            results.descriptive.(tag).n            = numel(gd);
            results.descriptive.(tag).n_subjects   = numel(unique(subject(idx)));
            results.descriptive.(tag).n_recordings = numel(unique(ID(idx)));
        end
    end

    %% ----------------------------------------------------------------
    %  Post-hoc: all pairwise cell contrasts (both contexts)
    %  Cell order matches plot layout: factor2 outer, factor1 inner
    %% ----------------------------------------------------------------
    cell_list = {};
    for j = 1:nF2
        for i = 1:nF1
            cell_list{end+1} = {factor1_levels{i}, factor2_levels{j}}; %#ok<AGROW>
        end
    end
    nCells   = numel(cell_list);
    pair_idx = nchoosek(1:nCells, 2);
    nPairs   = size(pair_idx, 1);

    pair_labels = cell(nPairs, 1);
    for k = 1:nPairs
        a = pair_idx(k,1); b = pair_idx(k,2);
        pair_labels{k} = sprintf('%s/%s vs %s/%s', ...
            cell_list{b}{1}, cell_list{b}{2}, ...
            cell_list{a}{1}, cell_list{a}{2});
    end

    for ctx = {'contextA','contextB'}
        context   = ctx{1};
        model_ctx = results.(context).model;
        coefs     = model_ctx.Coefficients;

        % Run post-hoc if any of the three effects is significant
        if results.(context).p_f1 < 0.05 || ...
           results.(context).p_f2 < 0.05 || ...
           results.(context).p_interaction < 0.05

            p_raw = nan(nPairs, 1);
            ests  = nan(nPairs, 1);
            SEs   = nan(nPairs, 1);

            for k = 1:nPairs
                ca = cell_list{pair_idx(k,1)};
                cb = cell_list{pair_idx(k,2)};
                cv = build_contrast(model_ctx, cb{1}, cb{2}, factor1_levels, factor2_levels) - ...
                     build_contrast(model_ctx, ca{1}, ca{2}, factor1_levels, factor2_levels);

                p_raw(k) = coefTest(model_ctx, cv);
                ests(k)  = cv * coefs.Estimate;
                SEs(k)   = sqrt(cv * model_ctx.CoefficientCovariance * cv');
            end

            p_holm = holm_bonferroni(p_raw);
            results.(context).posthoc.labels    = pair_labels;
            results.(context).posthoc.p_raw     = p_raw;
            results.(context).posthoc.p_holm    = p_holm;
            results.(context).posthoc.estimates = ests;
            results.(context).posthoc.SE        = SEs;
        end
    end

    %% ----------------------------------------------------------------
    %  Print to screen + save to text file
    %% ----------------------------------------------------------------
    save_dir = pwd;
    txt_path = fullfile(save_dir, sprintf('GLME2way_%s_%s_%s_%s_results.txt', ...
        distribution, measure_name, factor1_name, factor2_name));
    fid = fopen(txt_path, 'w');
    tee = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid}, 'UniformOutput', false);

    tee('\n========================================\n');
    tee('MEASURE: %s\n', measure_name);
    tee('Distribution: %s (link: %s)\n', upper(distribution), link);
    tee('FACTORS: %s (levels: %s) x %s (levels: %s)\n', ...
        factor1_name, strjoin(factor1_levels,', '), ...
        factor2_name, strjoin(factor2_levels,', '));
    tee('========================================\n');

    tee('\nDescriptive Statistics (per cell):\n');
    for j = 1:nF2
        for i = 1:nF1
            tag = matlab.lang.makeValidName([factor1_levels{i} '_' factor2_levels{j}]);
            d   = results.descriptive.(tag);
            tee('  %s x %s: M=%.3f, SD=%.3f, Mdn=%.3f, n=%d ripples, %d subjects, %d recordings\n', ...
                factor1_levels{i}, factor2_levels{j}, ...
                d.mean, d.std, d.median, d.n, d.n_subjects, d.n_recordings);
        end
    end

    ctx_display = { ...
        'CONTEXT A: Nested (ID within subject)', formulaA, 'contextA', pA_f1, pA_f2, pA_int; ...
        'CONTEXT B: Independent recordings',     formulaB, 'contextB', pB_f1, pB_f2, pB_int};

    for c = 1:2
        lbl     = ctx_display{c,1};
        formula = ctx_display{c,2};
        context = ctx_display{c,3};
        pF1     = ctx_display{c,4};
        pF2     = ctx_display{c,5};
        pI      = ctx_display{c,6};

        tee('\n--- %s ---\n', lbl);
        tee('Formula: %s\n', formula);
        tee('%s (main effect):     p=%.4f%s (omnibus Wald)\n', factor1_name, pF1, sig_star(pF1));
        tee('%s (main effect):     p=%.4f%s (omnibus Wald)\n', factor2_name, pF2, sig_star(pF2));
        tee('Interaction %sx%s: p=%.4f%s (omnibus Wald)\n', ...
            factor1_name, factor2_name, pI, sig_star(pI));
        tee('AIC=%.1f, BIC=%.1f, Dispersion=%.3f\n', ...
            results.(context).AIC, results.(context).BIC, results.(context).dispersion);

        if isfield(results.(context), 'posthoc')
            tee('Post-hoc pairwise (Holm-corrected):\n');
            ph = results.(context).posthoc;
            for k = 1:nPairs
                tee('  %s: Est=%.3f, SE=%.3f, p_raw=%.4f, p_holm=%.4f%s\n', ...
                    ph.labels{k}, ph.estimates(k), ph.SE(k), ...
                    ph.p_raw(k), ph.p_holm(k), sig_star(ph.p_holm(k)));
            end
        end
    end
    tee('========================================\n\n');
    fclose(fid);
    fprintf('Results text saved:\n  %s\n', txt_path);

    %% ----------------------------------------------------------------
    %  Visualization — grouped boxplots (Context A)
    %  factor2 = outer groups | factor1 = inner bars within group
    %% ----------------------------------------------------------------
    gap_within  = 0.7;
    gap_between = 1.8;

    x_pos  = zeros(nF2, nF1);
    cursor = 1;
    for j = 1:nF2
        for i = 1:nF1
            x_pos(j,i) = cursor;
            if i < nF1
                cursor = cursor + gap_within;
            end
        end
        cursor = cursor + gap_between;
    end

    fig = figure('Position', [100 100, max(700, 180*nF1*nF2), 620]);
    ax  = axes('Parent', fig);
    hold(ax, 'on');
    all_y_vals = [];

    % First pass: scatter dots + connecting lines
    for j = 1:nF2
        for i = 1:nF1
            xp       = x_pos(j,i);
            col_cell = cell_colors((j-1)*nF1 + i, :);

            idx_cell   = (f1 == factor1_levels{i}) & (f2 == factor2_levels{j});
            subj_list  = unique(subject(idx_cell));
            subj_means = zeros(numel(subj_list), 1);

            for s = 1:numel(subj_list)
                idx_subj      = idx_cell & (subject == subj_list(s));
                subj_means(s) = mean(y(idx_subj));
                ID_list       = unique(ID(idx_subj));
                rec_means     = zeros(numel(ID_list), 1);
                for r = 1:numel(ID_list)
                    rec_means(r) = mean(y(idx_subj & ID == ID_list(r)));
                end
                all_y_vals = [all_y_vals; rec_means; subj_means(s)]; %#ok<AGROW>
                x_jit = xp + 0.15*(rand(numel(rec_means),1)-0.5);
                scatter(ax, x_jit, rec_means, 30, col_cell, 'filled', 'MarkerFaceAlpha', 0.4);
                for r = 1:numel(rec_means)
                    plot(ax, [x_jit(r), xp], [rec_means(r), subj_means(s)], ...
                        'Color', [col_cell, 0.3], 'LineWidth', 0.5);
                end
            end
            x_sm = xp + 0.15*(rand(numel(subj_means),1)-0.5);
            scatter(ax, x_sm, subj_means, 120, col_cell, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        end
    end

    % Second pass: boxplots on top
    for j = 1:nF2
        for i = 1:nF1
            xp       = x_pos(j,i);
            col_cell = cell_colors((j-1)*nF1 + i, :);

            idx_cell   = (f1 == factor1_levels{i}) & (f2 == factor2_levels{j});
            subj_list  = unique(subject(idx_cell));
            subj_means = zeros(numel(subj_list), 1);
            for s = 1:numel(subj_list)
                subj_means(s) = mean(y(idx_cell & subject == subj_list(s)));
            end
           if numel(subj_means) >= 2
                boxplot(subj_means, 'Positions', xp, 'Width', 0.45, ...
                    'Colors', col_cell, 'Symbol', '', 'Parent', ax);
            end
        end
    end

    % Axes formatting
    if ~isempty(all_y_vals)
        yr = max(all_y_vals) - min(all_y_vals);
        if yr == 0; yr = 1; end
        ylim(ax, [min(all_y_vals) - 0.10*yr, max(all_y_vals) + 0.20*yr]);
    end

    group_centers = mean(x_pos, 2);
    xticks(ax, group_centers);
    xticklabels(ax, factor2_levels);
    xlim(ax, [min(x_pos(:)) - gap_within, max(x_pos(:)) + gap_within]);
    ylabel(ax, strrep(measure_name, '_', ' '));
    title(ax, sprintf('Context A: %s  [%s \times %s]  (%s)', ...
        strrep(measure_name,'_',' '), factor1_name, factor2_name, upper(distribution)));
    set(ax, 'FontSize', 12, 'LineWidth', 1);
    box(ax, 'off');

    % Sub-labels below each bar (factor1 level name, dimmed cell color)
    y_bot    = ax.YLim(1);
    y_offset = 0.04 * (ax.YLim(2) - ax.YLim(1));
    for j = 1:nF2
        for i = 1:nF1
            col_cell = cell_colors((j-1)*nF1 + i, :);
            text(ax, x_pos(j,i), y_bot - y_offset, factor1_levels{i}, ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'Color', col_cell * 0.75);
        end
    end

    % Legend: one entry per cell (factor2 – factor1)
    h_leg      = zeros(nF1*nF2, 1);
    leg_labels = cell(nF1*nF2, 1);
    for j = 1:nF2
        for i = 1:nF1
            k             = (j-1)*nF1 + i;
            col_cell      = cell_colors(k, :);
            h_leg(k)      = scatter(ax, nan, nan, 80, col_cell, 'filled', ...
                'MarkerEdgeColor','k','LineWidth',1.2);
            leg_labels{k} = sprintf('%s \x2013 %s', factor2_levels{j}, factor1_levels{i});
        end
    end
    legend(ax, h_leg, leg_labels, 'Location','northeast', 'FontSize', 10);

    % Significance brackets (Context A, Holm-corrected p < 0.05 only)
    if isfield(results.contextA, 'posthoc')
        add_significance_bars(ax, results.contextA.posthoc, x_pos, factor1_levels, factor2_levels);
    end

    hold(ax, 'off');

    %% ----------------------------------------------------------------
    %  Save figure + .mat
    %% ----------------------------------------------------------------
    base_name = sprintf('GLME2way_%s_%s_%s_%s_contextA', ...
        distribution, measure_name, factor1_name, factor2_name);
    fig_path  = fullfile(save_dir, base_name);

    exportgraphics(fig, [fig_path '.jpeg'], 'Resolution', 300);
    print(fig, [fig_path '.svg'], '-dsvg');
    fprintf('Figure saved:\n  %s.jpeg\n  %s.svg\n', fig_path, fig_path);

    mat_path = fullfile(save_dir, sprintf('GLME2way_%s_%s_%s_%s_results.mat', ...
        distribution, measure_name, factor1_name, factor2_name));
    save(mat_path, 'results');
    fprintf('Results .mat saved:\n  %s\n', mat_path);
end


%% ================================================================
%  Helper: extract omnibus p-values for f1, f2, and f1:f2
%  from a fitted GLME interaction model via Wald coefTest
%% ================================================================
function [p_f1, p_f2, p_int] = omnibus_pvalues(model, nF1, nF2)
    coef_names = model.CoefficientNames;
    nCoef      = numel(coef_names);

    % Identify indices for each term group
    idx_f1  = find(~cellfun(@isempty, regexp(coef_names, '^f1_', 'once')));
    idx_f2  = find(~cellfun(@isempty, regexp(coef_names, '^f2_', 'once')));
    idx_int = find(~cellfun(@isempty, regexp(coef_names, '^f1_.*:f2_', 'once')));

    % Build contrast matrices (one row per dummy coefficient in each term)
    C_f1  = build_term_contrast(nCoef, idx_f1);
    C_f2  = build_term_contrast(nCoef, idx_f2);
    C_int = build_term_contrast(nCoef, idx_int);

    p_f1  = coefTest(model, C_f1);
    p_f2  = coefTest(model, C_f2);
    if ~isempty(idx_int)
        p_int = coefTest(model, C_int);
    else
        p_int = NaN;
    end
end

function C = build_term_contrast(nCoef, term_idx)
    % Builds an identity-style contrast matrix selecting term_idx coefficients
    C = zeros(numel(term_idx), nCoef);
    for k = 1:numel(term_idx)
        C(k, term_idx(k)) = 1;
    end
end


%% ================================================================
%  Helper: build contrast vector for one cell mean
%  Reference levels: factor1_levels{1}, factor2_levels{1}
%% ================================================================
function c = build_contrast(model, f1_val, f2_val, f1_levels, f2_levels)
    coef_names = model.CoefficientNames;
    nCoef      = numel(coef_names);
    c          = zeros(1, nCoef);
    c(1)       = 1;   % intercept

    for k = 2:numel(f1_levels)
        idx = find(strcmp(coef_names, sprintf('f1_%s', f1_levels{k})), 1);
        if ~isempty(idx) && strcmp(f1_val, f1_levels{k}), c(idx) = 1; end
    end

    for k = 2:numel(f2_levels)
        idx = find(strcmp(coef_names, sprintf('f2_%s', f2_levels{k})), 1);
        if ~isempty(idx) && strcmp(f2_val, f2_levels{k}), c(idx) = 1; end
    end

    for k2 = 2:numel(f2_levels)
        for k1 = 2:numel(f1_levels)
            idx = find(strcmp(coef_names, sprintf('f1_%s:f2_%s', f1_levels{k1}, f2_levels{k2})), 1);
            if ~isempty(idx) && strcmp(f1_val, f1_levels{k1}) && strcmp(f2_val, f2_levels{k2})
                c(idx) = 1;
            end
        end
    end
end


%% ================================================================
%  Helper: significance brackets on axes
%% ================================================================
function add_significance_bars(ax, posthoc, x_pos, f1_levels, f2_levels)
    nF1 = numel(f1_levels);
    nF2 = numel(f2_levels);

    cell_x = containers.Map('KeyType','char','ValueType','double');
    for j = 1:nF2
        for i = 1:nF1
            cell_x(sprintf('%s/%s', f1_levels{i}, f2_levels{j})) = x_pos(j,i);
        end
    end

    y_range = ax.YLim(2) - ax.YLim(1);
    bar_h   = 0.04 * y_range;
    y_cur   = ax.YLim(2) + 0.02 * y_range;

    for k = 1:numel(posthoc.labels)
        if posthoc.p_holm(k) < 0.05
            parts = strsplit(posthoc.labels{k}, ' vs ');
            if numel(parts) ~= 2; continue; end
            if ~isKey(cell_x, parts{1}) || ~isKey(cell_x, parts{2}); continue; end
            x1 = cell_x(parts{1});
            x2 = cell_x(parts{2});
            plot(ax, [x1 x1 x2 x2], [y_cur, y_cur+bar_h, y_cur+bar_h, y_cur], 'k-', 'LineWidth', 1);
            text(ax, (x1+x2)/2, y_cur + bar_h*1.1, sig_star(posthoc.p_holm(k)), ...
                'HorizontalAlignment','center','FontSize',11);
            y_cur = y_cur + 2.5*bar_h;
        end
    end
    ylim(ax, [ax.YLim(1), y_cur + bar_h]);
end


%% ================================================================
%  Helper: Holm-Bonferroni correction
%% ================================================================
function p_holm = holm_bonferroni(p_raw)
    n        = numel(p_raw);
    [~, ord] = sort(p_raw);
    p_holm   = p_raw;
    for i = 1:n
        p_holm(ord(i)) = min(1, p_raw(ord(i)) * (n - i + 1));
    end
    for i = 2:n
        p_holm(ord(i)) = max(p_holm(ord(i)), p_holm(ord(i-1)));
    end
end


%% ================================================================
%  Helper: significance star string
%% ================================================================
function s = sig_star(p)
    if     p < 0.001; s = ' ***';
    elseif p < 0.01;  s = ' **';
    elseif p < 0.05;  s = ' *';
    else;             s = '';
    end
end