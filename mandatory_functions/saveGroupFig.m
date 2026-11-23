function saveGroupFig(tableRun, tableNonRun, groupLabels, figFolder)
    bands      = {'Theta (4-12 Hz)', 'Low Gamma (20-60 Hz)', 'High Gamma (60-100 Hz)'};
    dbCols     = {"medianTheta (db)",  "medianLgamma (db)",  "medianHgamma (db)"};
    zmodCols   = {"zmodTheta (au)",    "zmodLgamma (au)",    "zmodHgamma (au)"};
    yLabels    = {'Median power (dB)', 'Robust z-score (au)'};
    metricCols = {dbCols, zmodCols};
    figNames   = {'group_comparison_dB', 'group_comparison_zmod'};
    figTitles  = {'Theta & Gamma Power — dB', 'Theta & Gamma Power — Robust z-score'};
    condTables = {tableRun, tableNonRun};
    condLabels = {'Run', 'Non-run'};
    nGroups    = numel(groupLabels);
    palette    = lines(nGroups);
    figW       = max(520, 240 * nGroups * 2);

    for m = 1:2
        cols = metricCols{m};
        fig  = figure('Visible','off','Color','w','Position',[50 50 figW 820]);
        for b = 1:3
            for cond = 1:2
                ax = subplot(3, 2, (b-1)*2 + cond);
                hold(ax, 'on');
                for g = 1:nGroups
                    T    = condTables{cond};
                    mask = strcmp(T.ExpCond, groupLabels{g});
                    vals = T.(cols{b})(mask);
                    nAn  = numel(vals);
                    mu   = mean(vals,'omitnan');
                    sem  = std(vals,'omitnan') / sqrt(nAn);
                    bar(ax, g, mu, 0.55, 'FaceColor',palette(g,:), 'FaceAlpha',0.60, 'EdgeColor','none');
                    errorbar(ax, g, mu, sem, 'k-', 'LineWidth',1.6, 'CapSize',7);
                    scatter(ax, g+(rand(nAn,1)-0.5)*0.22, vals, 32, palette(g,:), 'filled', ...
                            'MarkerFaceAlpha',0.90, 'MarkerEdgeColor','k', 'LineWidth',0.4);
                end
                set(ax,'XTick',1:nGroups,'XTickLabel',groupLabels,'FontSize',9,...
                       'TickDir','out','XLim',[0.4, nGroups+0.6]);
                grid(ax,'on'); box(ax,'off');
                if b    == 1, title(ax, condLabels{cond}, 'FontSize',11,'FontWeight','bold'); end
                if cond == 1, ylabel(ax, sprintf('%s\n%s', bands{b}, yLabels{m}), 'FontSize',9); end
            end
        end
        sgtitle(figTitles{m}, 'FontSize',13,'FontWeight','bold');
        exportgraphics(fig, fullfile(figFolder, [figNames{m} '.png']), 'Resolution',200);
        savefig(fig, fullfile(figFolder, [figNames{m} '.fig']));
        close(fig);
        fprintf('  Figura grupal guardada: %s.png\n', figNames{m});
    end
end

