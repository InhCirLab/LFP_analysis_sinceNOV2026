function results = analyze_psd_lme_2way_lineplot(data, measure_name, group_var, group_order, error_type, colors_input)

    if nargin < 3, group_var    = 'expcond'; end
    if nargin < 4, group_order  = [];        end
    if nargin < 5, error_type   = 'sem';     end
    if nargin < 6, colors_input = [];        end

    layer_order = {'or','pyr','rad','slm'};
    nLayers     = numel(layer_order);

    groups_data = unique(string(data.(group_var)), 'stable');
    if isempty(group_order)
        groups = groups_data;
    else
        groups = string(group_order(:)');
    end
    nGroups = numel(groups);

    y       = data.(measure_name);
    expcond = categorical(data.(group_var), groups);
    layer   = categorical(data.layer, layer_order);
    subject = categorical(data.subject);
    ID      = categorical(data.ID);

    valid   = ~isnan(y) & isfinite(y);
    y       = y(valid); expcond = expcond(valid);
    layer   = layer(valid); subject = subject(valid);
    ID      = ID(valid);

    tbl    = table(y, expcond, layer, subject, ID);
    fmlaA  = 'y ~ expcond * layer + (1|subject) + (1|subject:ID)';
    modelA = fitlme(tbl, fmlaA);
    anovaA = anova(modelA);

    termsA = string(anovaA.Term);
    get_pA = @(t) anovaA.pValue(strcmp(termsA, t));
    pA_grp = get_pA(group_var);
    pA_lay = get_pA('layer');
    pA_int = get_pA([group_var ':layer']);

    gstats = struct();
    for g = 1:nGroups
        for l = 1:nLayers
            idx       = expcond == groups{g} & layer == layer_order{l};
            subj_list = unique(subject(idx));
            smeans    = arrayfun(@(s) mean(y(idx & subject==s)), subj_list);
            n         = numel(smeans);
            sem       = std(smeans) / sqrt(n);
            t95       = tinv(0.975, n-1);
            gstats(g,l).mean      = mean(smeans);
            gstats(g,l).sem       = sem;
            gstats(g,l).ci95      = t95 * sem;
            gstats(g,l).sd        = std(smeans);
            gstats(g,l).n         = n;
            gstats(g,l).smeans    = smeans;
        end
    end

    pair_idx  = nchoosek(1:nGroups, 2);
    nPairs    = size(pair_idx, 1);
    coefNames = modelA.CoefficientNames;
    coefs     = modelA.Coefficients;
    posthoc   = struct();

    for l = 1:nLayers
        p_raw = nan(nPairs,1);
        ests  = nan(nPairs,1);
        for k = 1:nPairs
            i1 = pair_idx(k,1); i2 = pair_idx(k,2);
            c  = zeros(1, numel(coefs.Estimate));
            n1  = sprintf('expcond_%s', groups{i1});
            n2  = sprintf('expcond_%s', groups{i2});
            nl1 = sprintf('expcond_%s:layer_%s', groups{i1}, layer_order{l});
            nl2 = sprintf('expcond_%s:layer_%s', groups{i2}, layer_order{l});
            ci  = @(nm) find(strcmp(coefNames, nm), 1);
            if i1 == 1
                if ~isempty(ci(n2)),  c(ci(n2))  =  1; end
                if ~isempty(ci(nl2)), c(ci(nl2)) =  1; end
            else
                if ~isempty(ci(n1)),  c(ci(n1))  = -1; end
                if ~isempty(ci(n2)),  c(ci(n2))  =  1; end
                if ~isempty(ci(nl1)), c(ci(nl1)) = -1; end
                if ~isempty(ci(nl2)), c(ci(nl2)) =  1; end
            end
            p_raw(k) = coefTest(modelA, c);
            ests(k)  = gstats(i2,l).mean - gstats(i1,l).mean;
        end
        p_holm = holm_bonferroni(p_raw);
        posthoc.(layer_order{l}).p_raw     = p_raw;
        posthoc.(layer_order{l}).p_holm    = p_holm;
        posthoc.(layer_order{l}).estimates = ests;
        posthoc.(layer_order{l}).contrasts = arrayfun(@(k) ...
            sprintf('%s vs %s', groups{pair_idx(k,2)}, groups{pair_idx(k,1)}), ...
            (1:nPairs)', 'UniformOutput', false);
    end

    results.model   = modelA;
    results.anova   = anovaA;
    results.posthoc = posthoc;

    %% --- Figure ---
    if isempty(colors_input)
        base_colors = [
            0.18, 0.38, 0.78;
            0.62, 0.18, 0.78;
            0.05, 0.62, 0.85;
            0.82, 0.12, 0.55;
            0.08, 0.48, 0.60;
        ];
        colors = base_colors(1:nGroups, :);
    else
        colors = colors_input;
    end

    x_pos = linspace(1.3, nLayers*0.85, nLayers);

    fig = figure('Position', [100 100 700 500], 'Color', 'w');
    ax  = axes('Parent', fig);
    hold(ax, 'on');

    for g = 1:nGroups
        col   = colors(g,:);
        means = arrayfun(@(l) gstats(g,l).mean, 1:nLayers);

        if strcmp(error_type, 'ci95')
            errs = arrayfun(@(l) gstats(g,l).ci95, 1:nLayers);
        else
            errs = arrayfun(@(l) gstats(g,l).sem,  1:nLayers);
        end

        plot(ax, x_pos, means, '-', 'Color', col, 'LineWidth', 1.8);

        for l = 1:nLayers
            plot(ax, [x_pos(l) x_pos(l)], ...
                [means(l)-errs(l) means(l)+errs(l)], ...
                '-', 'Color', col, 'LineWidth', 1.5);
            plot(ax, x_pos(l), means(l), 'o', 'MarkerSize', 8, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
        end
    end

    % if pA_int < 0.05
    %     y_range = ylim(ax);
    %     y_top   = y_range(2);
    %     step    = 0.06 * diff(y_range);
    %     for l = 1:nLayers
    %         ph    = posthoc.(layer_order{l});
    %         sig_k = find(ph.p_holm < 0.05);
    %         for ki = 1:numel(sig_k)
    %             k = sig_k(ki);
    %             text(ax, x_pos(l), y_top + step*ki, ...
    %                 sprintf('%s %s', ph.contrasts{k}, sig_stars(ph.p_holm(k))), ...
    %                 'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2]);
    %         end
    %     end
    % end

    set(ax, 'XLim', [x_pos(1)-0.3, x_pos(end)+0.3], ...
        'XTick', x_pos, 'XTickLabel', layer_order, 'FontSize', 12, ...
        'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);
    xlabel(ax, 'Layer',                      'FontSize', 13, 'FontWeight', 'bold');
    ylabel(ax, strrep(measure_name,'_',' '), 'FontSize', 13);

    if strcmp(error_type, 'ci95'), err_label = '95% CI'; else, err_label = 'SEM'; end

    % leg_h = gobjects(nGroups,1);
    % for g = 1:nGroups
    %     leg_h(g) = plot(ax, nan, nan, 'o-', 'Color', colors(g,:), ...
    %         'MarkerFaceColor', colors(g,:), 'LineWidth', 1.8);
    % end
    % legend(ax, leg_h, groups, 'Location', 'best', 'Box', 'off', 'FontSize', 11);

    % title(ax, sprintf('%s: p=%.3f%s  |  layer: p=%.3f%s  |  interaction: p=%.3f%s  |  error: %s', ...
    %     group_var, pA_grp, sig_stars(pA_grp), pA_lay, sig_stars(pA_lay), ...
    %     pA_int, sig_stars(pA_int), err_label), 'FontSize', 9, 'FontWeight', 'normal');

    hold(ax, 'off');

    %% --- Save text results ---
    save_dir  = pwd;
    base_name = sprintf('LME2way_line_%s_%s', measure_name, group_var);

    lines = {};
    lines{end+1} = '========================================';
    lines{end+1} = sprintf('LME RESULTS: %s', measure_name);
    lines{end+1} = sprintf('Formula:     %s', fmlaA);
    lines{end+1} = sprintf('Group order: %s', strjoin(groups, ' > '));
    lines{end+1} = sprintf('Date:        %s', datestr(now));
    lines{end+1} = '========================================';

    lines{end+1} = '';
    lines{end+1} = 'MODEL FIT';
    lines{end+1} = sprintf('  AIC = %.4f', modelA.ModelCriterion.AIC);
    lines{end+1} = sprintf('  BIC = %.4f', modelA.ModelCriterion.BIC);
    lines{end+1} = sprintf('  LogLikelihood = %.4f', modelA.LogLikelihood);

    lines{end+1} = '';
    lines{end+1} = 'FIXED EFFECTS (ANOVA)';
    lines{end+1} = sprintf('  %-30s  %6s  %8s  %8s  %8s  %s', ...
        'Term', 'df1', 'df2', 'F', 'p', 'Sig');
    lines{end+1} = ['  ' repmat('-', 1, 68)];
    for r = 1:height(anovaA)
        lines{end+1} = sprintf('  %-30s  %6d  %8.2f  %8.4f  %8.4f  %s', ...
            string(anovaA.Term(r)), anovaA.DF1(r), anovaA.DF2(r), ...
            anovaA.FStat(r), anovaA.pValue(r), sig_stars(anovaA.pValue(r)));
    end

    lines{end+1} = '';
    lines{end+1} = 'RANDOM EFFECTS';
    re = randomEffects(modelA);
    reNames = modelA.VariableInfo.Properties.RowNames;
    reStats = modelA.covarianceParameters;
    for r = 1:numel(reStats)
        lines{end+1} = sprintf('  Group %d variance: %.4f', r, reStats{r}(1));
    end
    lines{end+1} = sprintf('  Residual variance: %.4f', modelA.MSE);

    lines{end+1} = '';
    lines{end+1} = 'DESCRIPTIVE STATISTICS (subject-level means)';
    lines{end+1} = sprintf('  %-10s  %-5s  %8s  %8s  %8s  %8s  %5s', ...
        'Group', 'Layer', 'Mean', 'SD', 'SEM', '95% CI', 'N');
    lines{end+1} = ['  ' repmat('-', 1, 62)];
    for g = 1:nGroups
        for l = 1:nLayers
            d = gstats(g,l);
            lines{end+1} = sprintf('  %-10s  %-5s  %8.4f  %8.4f  %8.4f  %8.4f  %5d', ...
                groups{g}, layer_order{l}, d.mean, d.sd, d.sem, d.ci95, d.n);
        end
    end

    lines{end+1} = '';
    if pA_int < 0.05
        lines{end+1} = 'POST-HOC: within-layer between-group contrasts (Holm-corrected) [interaction SIGNIFICANT]';
    else
        lines{end+1} = 'POST-HOC: within-layer between-group contrasts (Holm-corrected) [interaction n.s. — exploratory]';
    end
    lines{end+1} = sprintf('  %-5s  %-25s  %10s  %10s  %10s  %s', ...
        'Layer', 'Contrast', 'Estimate', 'p_raw', 'p_holm', 'Sig');
    lines{end+1} = ['  ' repmat('-', 1, 68)];
    for l = 1:nLayers
        ph = posthoc.(layer_order{l});
        for k = 1:nPairs
            lines{end+1} = sprintf('  %-5s  %-25s  %10.4f  %10.4f  %10.4f  %s', ...
                layer_order{l}, ph.contrasts{k}, ph.estimates(k), ...
                ph.p_raw(k), ph.p_holm(k), sig_stars(ph.p_holm(k)));
        end
    end

    lines{end+1} = '';
    lines{end+1} = 'FIXED EFFECT COEFFICIENTS';
    lines{end+1} = sprintf('  %-45s  %10s  %10s  %10s  %10s', ...
        'Name', 'Estimate', 'SE', 't', 'p');
    lines{end+1} = ['  ' repmat('-', 1, 85)];
    for r = 1:height(coefs)
        lines{end+1} = sprintf('  %-45s  %10.4f  %10.4f  %10.4f  %10.4f  %s', ...
            coefs.Name{r}, coefs.Estimate(r), coefs.SE(r), ...
            coefs.tStat(r), coefs.pValue(r), sig_stars(coefs.pValue(r)));
    end

    lines{end+1} = '';
    lines{end+1} = '========================================';

    txt_path = fullfile(save_dir, [base_name '_results.txt']);
    fid = fopen(txt_path, 'w');
    for i = 1:numel(lines), fprintf(fid, '%s\n', lines{i}); end
    fclose(fid);

    exportgraphics(fig, fullfile(save_dir, [base_name '.jpeg']), 'Resolution', 300);
    print(fig, fullfile(save_dir, [base_name '.svg']), '-dsvg');
    save(fullfile(save_dir, [base_name '_results.mat']), 'results');
    fprintf('\nSaved: %s (.jpeg / .svg / _results.mat / _results.txt)\n', fullfile(save_dir, base_name));
end

function p_holm = holm_bonferroni(p_raw)
    n = numel(p_raw);
    [sp, idx] = sort(p_raw);
    p_holm = zeros(n,1);
    for i = 1:n, p_holm(idx(i)) = min(sp(i)*(n-i+1), 1); end
    for i = 2:n, p_holm(idx(i)) = max(p_holm(idx(i)), p_holm(idx(i-1))); end
end

function s = sig_stars(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = '';
    end
end