%% ========================================================================
%  FUNCIONES AUXILIARES — PARTE II
%% ========================================================================

function savePSFig(T, titleStr, yUnit, layerNames, groupLabels, figFolder, fname)
    nL      = numel(layerNames);
    nG      = numel(groupLabels);
    palette = lines(nG);
    figH    = max(300, 200 * nL);
    figW    = max(400, 180 * nG);

    fig = figure('Visible','off','Color','w','Position',[50 50 figW figH]);

    for L = 1:nL
        ax = subplot(nL, 1, L);
        hold(ax, 'on');

        for g = 1:nG
            mask = strcmp(T.Group, groupLabels{g});
            vals = T.(layerNames{L})(mask);
            nAn  = sum(mask);
            mu   = mean(vals,'omitnan');
            sem  = std(vals,'omitnan') / sqrt(nAn);

            bar(ax, g, mu, 0.55, 'FaceColor',palette(g,:), 'FaceAlpha',0.60, 'EdgeColor','none');
            errorbar(ax, g, mu, sem, 'k-', 'LineWidth',1.5, 'CapSize',6);
            scatter(ax, g+(rand(nAn,1)-0.5)*0.2, vals, 28, palette(g,:), 'filled', ...
                    'MarkerFaceAlpha',0.85, 'MarkerEdgeColor','k', 'LineWidth',0.4);
        end

        set(ax,'XTick',1:nG,'XTickLabel',groupLabels,'FontSize',9,...
               'TickDir','out','XLim',[0.4, nG+0.6]);
        ylabel(ax, sprintf('%s\n(%s)', upper(layerNames{L}), yUnit), 'FontSize',9);
        grid(ax,'on'); box(ax,'off');

        if L == 1
            title(ax, strrep(titleStr,'_','\_'), 'FontSize',10,'FontWeight','bold');
        end
    end

    exportgraphics(fig, fullfile(figFolder, [fname '.png']), 'Resolution',150);
    close(fig);
end