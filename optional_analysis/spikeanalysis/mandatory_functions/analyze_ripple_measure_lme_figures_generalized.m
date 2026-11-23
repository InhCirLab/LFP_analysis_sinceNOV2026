function results = analyze_ripple_measure_lme_figures_generalized(data, measure_name, group_var, colors, group_order, display_names)
% ANALYZE_RIPPLE_MEASURE_LME_FIGURES_GENERALIZED
% Linear mixed-effects models assuming normal distribution.
% Groups are inferred automatically from data.(group_var).
% Context A: ID nested in subject
% Context B: ID treated as independent
%
% Figures and results are saved to pwd automatically.
%
% INPUT:
%   data         - table with columns: ID, subject, <group_var>, measure_name
%   measure_name - string, name of the measure column
%   group_var    - string, name of the grouping variable (default: 'expcond')
%   colors       - (optional) nGroups x 3 RGB matrix of group colors
%   group_order  - (optional) cell/string array specifying left-to-right plot order of groups

    if nargin < 3, group_var = 'expcond'; end
    if nargin < 6, display_names = [];        end  % add this line
    %% Extract and clean data
    y       = data.(measure_name);
    groups  = unique(string(data.(group_var)), 'stable');
    nGroups = numel(groups);

    if nargin >= 5 && ~isempty(group_order)
        group_order = string(group_order);
        assert(numel(group_order) == nGroups && all(ismember(group_order, groups)), ...
            'group_order must contain exactly the same group names as in %s.', group_var);
        groups = group_order;
    end
    % Resolve display labels (for x-axis only; internal logic uses groups{})
    if isempty(display_names)
        plot_labels = groups;
    else
        display_names = string(display_names);
        assert(numel(display_names) == nGroups, ...
            'display_names must have the same number of entries as groups (%d).', nGroups);
        plot_labels = display_names;
    end
    expcond = categorical(data.(group_var), groups);
    subject = categorical(data.subject);
    ID      = categorical(data.ID);

    valid   = ~isnan(y) & isfinite(y);
    y       = y(valid);
    expcond = expcond(valid);
    subject = subject(valid);
    ID      = ID(valid);

    tbl = table(y, expcond, subject, ID);

    %% CONTEXT A: ID nested within subject
    modelA = fitlme(tbl, 'y ~ expcond + (1|subject) + (1|subject:ID)');
    anovaA = anova(modelA);
    pA     = anovaA.pValue(2);

    results.contextA.model   = modelA;
    results.contextA.p_value = pA;
    results.contextA.anova   = anovaA;
    results.contextA.AIC     = modelA.ModelCriterion.AIC;
    results.contextA.BIC     = modelA.ModelCriterion.BIC;

    %% CONTEXT B: ID independent
    modelB = fitlme(tbl, 'y ~ expcond + (1|ID)');
    anovaB = anova(modelB);
    pB     = anovaB.pValue(2);

    results.contextB.model   = modelB;
    results.contextB.p_value = pB;
    results.contextB.anova   = anovaB;
    results.contextB.AIC     = modelB.ModelCriterion.AIC;
    results.contextB.BIC     = modelB.ModelCriterion.BIC;

    %% Descriptive statistics
    for i = 1:nGroups
        idx      = expcond == groups{i};
        grp_data = y(idx);
        results.descriptive.(groups{i}).mean         = mean(grp_data);
        results.descriptive.(groups{i}).std          = std(grp_data);
        results.descriptive.(groups{i}).n            = numel(grp_data);
        results.descriptive.(groups{i}).n_subjects   = numel(unique(subject(idx)));
        results.descriptive.(groups{i}).n_recordings = numel(unique(ID(idx)));
    end

    %% Post-hoc contrasts (all pairwise, reference = groups{1})
    pair_idx        = nchoosek(1:nGroups, 2);
    nPairs          = size(pair_idx, 1);
    contrast_labels = cell(nPairs, 1);
    for k = 1:nPairs
        contrast_labels{k} = sprintf('%s vs %s', groups{pair_idx(k,2)}, groups{pair_idx(k,1)});
    end

    coefNames = modelA.CoefficientNames;

    for ctx = {'contextA', 'contextB'}
        context   = ctx{1};
        model_ctx = results.(context).model;
        coefs     = model_ctx.Coefficients;

        if results.(context).p_value < 0.05
            p_raw = nan(nPairs, 1);
            ests  = nan(nPairs, 1);
            SEs   = nan(nPairs, 1);

            for k = 1:nPairs
                i1 = pair_idx(k, 1);
                i2 = pair_idx(k, 2);

                c = zeros(1, numel(coefs.Estimate));
                if i1 == 1
                    coef2_name = sprintf('expcond_%s', groups{i2});
                    ci2 = find(strcmp(coefNames, coef2_name), 1);
                    if ~isempty(ci2), c(ci2) = 1; end
                    ests(k) = coefs.Estimate(ci2);
                    SEs(k)  = coefs.SE(ci2);
                else
                    coef1_name = sprintf('expcond_%s', groups{i1});
                    coef2_name = sprintf('expcond_%s', groups{i2});
                    ci1 = find(strcmp(coefNames, coef1_name), 1);
                    ci2 = find(strcmp(coefNames, coef2_name), 1);
                    if ~isempty(ci1), c(ci1) = -1; end
                    if ~isempty(ci2), c(ci2) =  1; end
                    ests(k) = coefs.Estimate(ci2) - coefs.Estimate(ci1);
                    SEs(k)  = sqrt(coefs.SE(ci1)^2 + coefs.SE(ci2)^2);
                end
                p_raw(k) = coefTest(model_ctx, c);
            end

            p_holm = holm_bonferroni(p_raw);
            results.(context).posthoc.contrasts = contrast_labels;
            results.(context).posthoc.p_raw     = p_raw;
            results.(context).posthoc.p_holm    = p_holm;
            results.(context).posthoc.estimates = ests;
            results.(context).posthoc.SE        = SEs;
        end
    end

    %% Display results + save to text file
    save_dir  = pwd;
    txt_path  = fullfile(save_dir, sprintf('LME_%s_%s_results.txt', measure_name, group_var));
    fid       = fopen(txt_path, 'w');
    tee       = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid}, 'UniformOutput', false);

    tee('\n========================================\n');
    tee('MEASURE: %s\n', measure_name);
    tee('========================================\n');
    tee('\nDescriptive Statistics:\n');
    for i = 1:nGroups
        tee('  %s: M=%.3f, SD=%.3f, n=%d (ripples), %d subjects, %d recordings\n', ...
            groups{i}, ...
            results.descriptive.(groups{i}).mean, ...
            results.descriptive.(groups{i}).std, ...
            results.descriptive.(groups{i}).n, ...
            results.descriptive.(groups{i}).n_subjects, ...
            results.descriptive.(groups{i}).n_recordings);
    end

    for ctx = {'contextA', 'contextB'}
        context = ctx{1};
        if strcmp(context, 'contextA')
            tee('\n--- CONTEXT A: Nested structure (ID within subject) ---\n');
            tee('Formula: y ~ %s + (1|subject) + (1|subject:ID)\n', group_var);
            thisAnova = anovaA; thisP = pA;
        else
            tee('\n--- CONTEXT B: Independent recordings ---\n');
            tee('Formula: y ~ %s + (1|ID)\n', group_var);
            thisAnova = anovaB; thisP = pB;
        end
        tee('Fixed effect (%s): F(%d,%d)=%.2f, p=%.4f\n', group_var, ...
            thisAnova.DF1(2), thisAnova.DF2(2), thisAnova.FStat(2), thisP);
        tee('AIC=%.1f, BIC=%.1f\n', results.(context).AIC, results.(context).BIC);
        if thisP < 0.05 && isfield(results.(context), 'posthoc')
            tee('Post-hoc (Holm-corrected):\n');
            for k = 1:nPairs
                tee('  %s: Est=%.3f, SE=%.3f, p_raw=%.4f, p_holm=%.4f %s\n', ...
                    contrast_labels{k}, ...
                    results.(context).posthoc.estimates(k), ...
                    results.(context).posthoc.SE(k), ...
                    results.(context).posthoc.p_raw(k), ...
                    results.(context).posthoc.p_holm(k), ...
                    iff(results.(context).posthoc.p_holm(k) < 0.05, '*', ''));
            end
        end
    end
    tee('========================================\n\n');
    fclose(fid);
    fprintf('Results text saved:\n  %s\n', txt_path);

    %% Visualization (Context A)
    if nargin < 4 || isempty(colors)
        colors = lines(nGroups);
    else
        assert(size(colors,1) == nGroups && size(colors,2) == 3, ...
            'colors must be an nGroups x 3 RGB matrix (%d rows expected).', nGroups);
    end
    x_positions = 1:nGroups;

    fig = figure('Position', [100 100 200*nGroups + 200, 600]);
    ax  = axes('Parent', fig);
    hold(ax, 'on');
    all_y = [];

    % First pass: scatter + lines
    for i = 1:nGroups
        idx_grp    = expcond == groups{i};
        subj_list  = unique(subject(idx_grp));
        x_pos      = x_positions(i);
        subj_means = zeros(numel(subj_list), 1);

        for s = 1:numel(subj_list)
            idx_subj      = idx_grp & subject == subj_list(s);
            subj_means(s) = mean(y(idx_subj));
            ID_list       = unique(ID(idx_subj));
            rec_means     = zeros(numel(ID_list), 1);
            for r = 1:numel(ID_list)
                rec_means(r) = mean(y(idx_subj & ID == ID_list(r)));
            end
            all_y    = [all_y; rec_means; subj_means(s)]; %#ok<AGROW>
            x_jitter = x_pos + 0.15 * (rand(numel(rec_means), 1) - 0.5);
            scatter(ax, x_jitter, rec_means, 30, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.4);
            for r = 1:numel(rec_means)
                plot(ax, [x_jitter(r), x_pos], [rec_means(r), subj_means(s)], ...
                    'Color', [colors(i,:) 0.3], 'LineWidth', 0.5);
            end
        end
        scatter(ax, x_pos + 0.15*(rand(numel(subj_means),1)-0.5), subj_means, 120, colors(i,:), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end

    % Second pass: boxplots
    for i = 1:nGroups
        idx_grp    = expcond == groups{i};
        subj_list  = unique(subject(idx_grp));
        subj_means = zeros(numel(subj_list), 1);
        for s = 1:numel(subj_list)
            subj_means(s) = mean(y(idx_grp & subject == subj_list(s)));
        end
        if ~isempty(subj_means)
            boxplot(subj_means, 'Positions', x_positions(i), 'Width', 0.3, ...
                'Colors', colors(i,:), 'Symbol', '', 'Parent', ax);
        end
    end

    y_range = max(all_y) - min(all_y);
    ylim(ax, [min(all_y) - 0.10*y_range, max(all_y) + 0.15*y_range]);
    xlim(ax, [0.5, nGroups + 0.5]);
    xticks(ax, x_positions);
    % xticklabels(ax, groups);
    xticklabels(ax, plot_labels);
    ylabel(ax, strrep(measure_name, '_', ' '));
    title(ax, sprintf('Context A: %s', strrep(measure_name, '_', ' ')));
    set(ax, 'FontSize', 12, 'LineWidth', 1);
    box(ax, 'off');
    hold(ax, 'off');

    %% Save figure and results
    base_name = sprintf('LME_%s_%s_contextA', measure_name, group_var);
    fig_path  = fullfile(save_dir, base_name);

    exportgraphics(fig, [fig_path '.jpeg'], 'Resolution', 300);
    print(fig, [fig_path '.svg'], '-dsvg');
    fprintf('Figure saved:\n  %s.jpeg\n  %s.svg\n', fig_path, fig_path);

    mat_path = fullfile(save_dir, sprintf('LME_%s_%s_results.mat', measure_name, group_var));
    save(mat_path, 'results');
    fprintf('Results saved:\n  %s\n', mat_path);
end

function p_holm = holm_bonferroni(p_raw)
    [n, ~] = size(p_raw);
    [sorted_p, idx] = sort(p_raw);
    p_holm = zeros(n, 1);
    for i = 1:n
        p_holm(idx(i)) = min(sorted_p(i) * (n - i + 1), 1);
    end
    for i = 2:n
        p_holm(idx(i)) = max(p_holm(idx(i)), p_holm(idx(i-1)));
    end
end

function iff_result = iff(condition, true_val, false_val)
    if condition
        iff_result = true_val;
    else
        iff_result = false_val;
    end
end