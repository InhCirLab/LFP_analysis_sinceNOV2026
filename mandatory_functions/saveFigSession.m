function saveFigSession(recName, expCond, runPow, nonrunPow, figFolder)
    bands     = {'Theta (4-12)', 'LGamma (20-60)', 'HGamma (60-100)'};
    runMed    = [median(runPow.dbT),    median(runPow.dbL),    median(runPow.dbH)];
    nonrunMed = [median(nonrunPow.dbT), median(nonrunPow.dbL), median(nonrunPow.dbH)];
    runAll    = [runPow.dbT(:),    runPow.dbL(:),    runPow.dbH(:)];
    nonrunAll = [nonrunPow.dbT(:), nonrunPow.dbL(:), nonrunPow.dbH(:)];
    x   = 1:3;
    fig = figure('Visible','off','Color','w','Position',[100 100 560 380]);
    h1  = bar(x-0.2, runMed,    0.35, 'FaceColor',[0.20 0.50 0.80], 'FaceAlpha',0.75);
    hold on;
    h2  = bar(x+0.2, nonrunMed, 0.35, 'FaceColor',[0.80 0.30 0.20], 'FaceAlpha',0.75);
    for b = 1:3
        scatter((b-0.2)+(rand(size(runAll,1),1)-0.5)*0.08,    runAll(:,b),    18, [0.10 0.30 0.65], 'filled', 'MarkerFaceAlpha',0.5);
        scatter((b+0.2)+(rand(size(nonrunAll,1),1)-0.5)*0.08, nonrunAll(:,b), 18, [0.65 0.10 0.10], 'filled', 'MarkerFaceAlpha',0.5);
    end
    set(gca,'XTick',x,'XTickLabel',bands,'FontSize',9);
    ylabel('Median power (dB)');
    title(sprintf('[%s] %s', expCond, strrep(recName,'_','\_')), 'FontSize',9);
    legend([h1 h2],{'Run','Non-run'},'Location','best','FontSize',8);
    grid on; box off;
    exportgraphics(fig, fullfile(figFolder, [expCond '_' recName '_summary.png']), 'Resolution',150);
    close(fig);
end