
%% SCRIPT FOR SIMPLIFIED CSD ANALYSIS

%% 1. Analyze data based on rad-sink vs slm-sink (obtained with criteria) 

% 6.1. One-way ANOVA (for all ripples)

radsinkdata = all_ripples(all_ripples.sinktype == "rad", :);
slmsinkdata = all_ripples(all_ripples.sinktype == "slm", :);
othersinkdata = all_ripples(all_ripples.sinktype == "other", :);

% Colors for males (SS, OS, OD)
% 
% groupcolors = [0, 0, 128;...
%               132,197,247;...
%               7, 126, 149]/255;

% Colors for females (P, E, M, D)
% 
% groupcolors = [229,159,50;...
%                 50,86,23;...
%                 31,58,133;...
%                 103,0,43]/255;

% Colors for XX vs XY:
% % 
% groupcolors = [184,73,26;...
%                 255,200,169]/255;


% FCG males:
% 
% groupcolors = [70, 129, 49;...
%         77, 225, 18]/255;

% FCG females OVX:
% 
% groupcolors = [224, 101, 18;...
%         255, 216, 190]/255;

% Combined analysis:

% XX vs XY (both males and females):
% groupcolors = [0, 0, 0;...
%     162, 160, 156]/255;

% Males vs females:

% groupcolors = [70, 129, 49;...
%         184,73,26]/255;

% sexchr vs sexdif
groupcolors = [
    114, 153, 93;   % male-XX    (blue)
    60,67, 14;   % male-XY    (red)
    124, 19, 194;   % female-XX  
    234, 135, 173;   % female-XY  
]/255;
% Now analysis are performed excluding outliers:

%% Measures for rad-sink
all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'PeakNormPow', 'expcond', groupcolors, ["XX", "XY"]);
all_peakamp = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'PeakAmp', 'expcond', groupcolors,  ["XX", "XY"]);
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'PeakFq', 'expcond', groupcolors,  ["XX", "XY"]);
all_duration = analyze_ripple_measure_glme_generalized(radsinkdata, 'Duration', 'gamma', 'expcond', groupcolors,  ["XX", "XY"]);
all_entropy = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'spectralEntropy', 'expcond', groupcolors,  ["XX", "XY"]);
all_interriplefq = analyze_ripple_measure_glme_generalized(radsinkdata, 'interrippleFq', 'weibull', 'expcond', groupcolors, ["XX", "XY"]);

%% Measures for slm-sink
all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'PeakNormPow', 'expcond', groupcolors, ["XX", "XY"]);
all_peakamp = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'PeakAmp', 'expcond', groupcolors,  ["XX", "XY"]);
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'PeakFq', 'expcond', groupcolors,  ["XX", "XY"]);
all_duration = analyze_ripple_measure_glme_generalized(slmsinkdata, 'Duration', 'gamma', 'expcond', groupcolors,  ["XX", "XY"]);
all_entropy = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'spectralEntropy', 'expcond', groupcolors, ["XX", "XY"]);
all_interriplefq = analyze_ripple_measure_glme_generalized(slmsinkdata, 'interrippleFq', 'weibull', 'expcond', groupcolors,  ["XX", "XY"]);

%% Measures for other-sink
all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'PeakNormPow', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_peakamp = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'PeakAmp', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'PeakFq', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_duration = analyze_ripple_measure_glme_generalized(othersinkdata, 'Duration', 'gamma', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_entropy = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'spectralEntropy', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_interriplefq = analyze_ripple_measure_glme_generalized(othersinkdata, 'interrippleFq', 'weibull', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});

%% MEasures for all ripples

all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(all_ripples, 'PeakNormPow', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_peakamp = analyze_ripple_measure_lme_figures_generalized(all_ripples, 'PeakAmp', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(all_ripples, 'PeakFq', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_duration = analyze_ripple_measure_glme_generalized(all_ripples, 'Duration', 'gamma', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_entropy = analyze_ripple_measure_lme_figures_generalized(all_ripples, 'spectralEntropy', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});
all_interriplefq = analyze_ripple_measure_glme_generalized(all_ripples, 'interrippleFq', 'weibull', 'expcond', groupcolors, ["SS", "OS", "OD"], {"Ctrl", "ORX", "ORX+DHT"});

%%
%% Measures for rad-sink
all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'PeakNormPow', 'expcond', groupcolors, ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_peakamp = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'PeakAmp', 'expcond', groupcolors, ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'PeakFq', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_duration = analyze_ripple_measure_glme_generalized(radsinkdata, 'Duration', 'gamma', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_entropy = analyze_ripple_measure_lme_figures_generalized(radsinkdata, 'spectralEntropy', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_interriplefq = analyze_ripple_measure_glme_generalized(radsinkdata, 'interrippleFq', 'weibull', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);

%% Measures for slm-sink
all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'PeakNormPow', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Diestrus"]);
all_peakamp = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'PeakAmp', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Diestrus"]);
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'PeakFq', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Diestrus"]);
all_duration = analyze_ripple_measure_glme_generalized(slmsinkdata, 'Duration', 'gamma', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Diestrus"]);
all_entropy = analyze_ripple_measure_lme_figures_generalized(slmsinkdata, 'spectralEntropy', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Diestrus"]);
all_interriplefq = analyze_ripple_measure_glme_generalized(slmsinkdata, 'interrippleFq', 'weibull', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Diestrus"]);

%% Measures for other-sink
all_peaknormpow = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'PeakNormPow', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_peakamp = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'PeakAmp', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_peakfreq = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'PeakFq', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_duration = analyze_ripple_measure_glme_generalized(othersinkdata, 'Duration', 'gamma', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_entropy = analyze_ripple_measure_lme_figures_generalized(othersinkdata, 'spectralEntropy', 'expcond', groupcolors,  ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);
all_interriplefq = analyze_ripple_measure_glme_generalized(othersinkdata, 'interrippleFq', 'weibull', 'expcond', groupcolors, ["Proestrus", "Estrus", "Metestrus", "Diestrus"]);

%% 2-way ANOVA for FCG combinations

groupcolors = [
    70, 129, 49;... % male   + XX
    77, 225, 18;... % male   + XY   
    184,73,26;...   % female + XX
    255,200,169]/255;% female + XY

% 
% results = analyze_ripple_2factor_lme(data, 'peak_freq', ...
%     'expcond', {'XX','XY'}, ...
%     'sex',     {'male','female'}, ...
%     cell_colors);


cyclingripples = all_ripples;
cyclingripples.expcond = repmat("cycling", size(all_ripples, 1), 1);

ovxripples = all_ripples;
ovxripples.expcond = repmat("OVX", size(ovxripples, 1), 1);

femaleripples.sex = repmat("female", size(femaleripples,1), 1);
femaleripples.expcond = femaleripples.sexchr;
femaleripples = removevars(femaleripples, "sexchr");

maleripples.sex = repmat("male", size(maleripples, 1), 1);
maleripples.ID = categorical(maleripples.ID);
maleripples.subject = categorical(maleripples.subject);

all_ripples = [femaleripples; maleripples];

%%  Rearrange data for females cycling vs females ovx

ovxripples.ID = categorical(ovxripples.ID);
ovxripples.subject = categorical(ovxripples.subject);
femaleripples.expcond = categorical(femaleripples.expcond);
all_ripples = [cyclingripples; ovxripples];

%%

radsinkdata = all_ripples(all_ripples.sinktype == "rad", :);
slmsinkdata = all_ripples(all_ripples.sinktype == "slm", :);
othersinkdata = all_ripples(all_ripples.sinktype == "other", :);

%% 
groupcolors = [184,73,26;...
                255,200,169;...
                224, 101, 18;...
             255, 216, 190]/255;

%% All ripples
all_peaknormpow = analyze_ripple_2factor_lme(all_ripples, 'PeakNormPow', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_peakamp = analyze_ripple_2factor_lme(all_ripples, 'PeakAmp', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_peakfreq = analyze_ripple_2factor_lme(all_ripples, 'PeakFq', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_duration = analyze_ripple_2factor_glme(all_ripples, 'Duration', 'gamma', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_entropy = analyze_ripple_2factor_lme(all_ripples, 'spectralEntropy', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_interriplefq = analyze_ripple_2factor_glme(all_ripples, 'interrippleFq', 'weibull', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);

%% Rad ripples

all_peaknormpow = analyze_ripple_2factor_lme(radsinkdata, 'PeakNormPow', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_peakamp = analyze_ripple_2factor_lme(radsinkdata, 'PeakAmp', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_peakfreq = analyze_ripple_2factor_lme(radsinkdata, 'PeakFq', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_duration = analyze_ripple_2factor_glme(radsinkdata, 'Duration', 'gamma', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_entropy = analyze_ripple_2factor_lme(radsinkdata, 'spectralEntropy', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_interriplefq = analyze_ripple_2factor_glme(radsinkdata, 'interrippleFq', 'weibull', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);

%% SLM ripples

all_peaknormpow = analyze_ripple_2factor_lme(slmsinkdata, 'PeakNormPow', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_peakamp = analyze_ripple_2factor_lme(slmsinkdata, 'PeakAmp', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_peakfreq = analyze_ripple_2factor_lme(slmsinkdata, 'PeakFq', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_duration = analyze_ripple_2factor_glme(slmsinkdata, 'Duration', 'gamma', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_entropy = analyze_ripple_2factor_lme(slmsinkdata, 'spectralEntropy', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
all_interriplefq = analyze_ripple_2factor_glme(slmsinkdata, 'interrippleFq', 'weibull', 'sexchr', {'XX','XY'}, 'expcond',{'male','female'}, groupcolors);
%% otherripples

all_peaknormpow = analyze_ripple_2factor_lme(othersinkdata, 'PeakNormPow', 'sexchr', {'XX','XY'}, 'expcond',{'cycling','OVX'}, groupcolors);
all_peakamp = analyze_ripple_2factor_lme(othersinkdata, 'PeakAmp', 'sexchr', {'XX','XY'}, 'expcond',{'cycling','OVX'}, groupcolors);
all_peakfreq = analyze_ripple_2factor_lme(othersinkdata, 'PeakFq', 'sexchr', {'XX','XY'}, 'expcond',{'cycling','OVX'}, groupcolors);
all_duration = analyze_ripple_2factor_glme(othersinkdata, 'Duration', 'gamma', 'sexchr', {'XX','XY'}, 'expcond',{'cycling','OVX'}, groupcolors);
all_entropy = analyze_ripple_2factor_lme(othersinkdata, 'spectralEntropy', 'sexchr', {'XX','XY'}, 'expcond',{'cycling','OVX'}, groupcolors);
all_interriplefq = analyze_ripple_2factor_glme(othersinkdata, 'interrippleFq', 'weibull', 'sexchr', {'XX','XY'}, 'expcond',{'cycling','OVX'}, groupcolors);


%% =========================================================================
%  Legend figure
%% =========================================================================


%% Standalone Legend Script
% Configure these two variables depending on your experiment:

% Colors for males (SS, OS, OD)
% 
% groupcolors = [0, 0, 128;...
%               132,197,247;...
%               7, 126, 149]/255;

% Colors for females (P, E, M, D)
% 
% groupcolors = [229,159,50;...
%                 50,86,23;...
%                 31,58,133;...
%                 103,0,43]/255;

% Colors for XX vs XY:
% % 
% groupcolors = [184,73,26;...
%                 255,200,169]/255;


% FCG males:
% 
% groupcolors = [70, 129, 49;...
%         77, 225, 18]/255;

% FCG females OVX:
% % 
% groupcolors = [224, 101, 18;...
%         255, 216, 190]/255;

% Combined analysis:

% XX vs XY (both males and females):
% groupcolors = [0, 0, 0;...
%     162, 160, 156]/255;

% Males vs females:

% groupcolors = [70, 129, 49;...
%         184,73,26]/255;

% sexchr vs sexdif
% groupcolors = [
%     114, 153, 93;   % male-XX    (blue)
%     60,67, 14;   % male-XY    (red)
%     124, 19, 194;   % female-XX  
%     234, 135, 173;   % female-XY  
% ]/255;

group_labels = {'XX', 'XY'};  % One label per group

% --- Generate legend figure (no changes needed below) ---

figure;
hold on;

for i = 1:numel(group_labels)
    plot(nan, nan, 'o', ...
        'MarkerFaceColor', groupcolors(i,:), ...
        'MarkerEdgeColor', groupcolors(i,:), ...
        'MarkerSize', 8);
end

lgd = legend(group_labels, 'Location', 'best');
lgd.Box = 'off';
lgd.FontSize = 11;

axis off;

%% Save
exportgraphics(gcf, 'legend.svg', 'ContentType', 'vector');
exportgraphics(gcf, 'legend.jpeg', 'Resolution', 300);

%% Extract UMAP

% filteredripples = filteredripples(~rowidx, :);
filteredRipples = table2array(filteredripples);
if size(filteredRipples, 2) == 187
    midpoint = round(size(filteredRipples, 2)/2);
    filteredRipples_centered = filteredRipples(:, midpoint-63:midpoint+63); % 127 samples
else
    filteredRipples_centered = filteredRipples;
end

% 8.2. Apply UMAP with optimal dimensionality

% [reduction, umap] = run_umap(filteredRipples_centered, ...
%     'n_components', 4, 'metric', 'euclidean', 'min_dist', 0.1);


%% Repeat figures but with clim constrained by LME results:

% % 8.3. Figure per groups
% data = all_ripples;
% % Extract group labels from data table
% if istable(data) && ismember('sinktype', data.Properties.VariableNames)
%     groupLabels = data.sinktype;
% elseif isstruct(data) && isfield(data, 'sinktype')
%     groupLabels = data.sinktype;
% else
%     error('Variable "expcond" not found in data. Please ensure data contains group labels.');
% end
% 
% groupLabels = categorical(groupLabels, {'rad','slm','other'}, 'Ordinal', true);
% groups = unique(groupLabels);
% % groups = categorical(groups, {'SS','OS','OD'}, 'Ordinal', true);
% colors = lines(numel(groups));
% 
% figure(1); hold on;
% for i = 1:3
%     idx = groupLabels == groups(i);
%     scatter(reduction(idx,1), reduction(idx,2), 10, colors(i,:), 'filled','MarkerFaceAlpha', 0.2);
% end
% xlabel('UMAP1');
% ylabel('UMAP2');
% title('UMAP Embedding of Sharp-Wave Ripple Measures');
% legend(groups, 'Location', 'best');
% % grid on;
% hold off;


%% Plot 3: Group Distribution per Sink Type
fig3 = figure('Position', [100, 100, 1400, 800]);

% --- Setup with FIXED ORDER ---
expcond  = categorical(all_ripples.expcond);
sinktype = all_ripples.sinktype;

% Enforce explicit ordering as categoricals
% group_order = categorical({'SS', 'OS', 'OD'});
group_order = categorical({'XX', 'XY'});
sink_order  = categorical({'rad', 'slm', 'other'});
n_groups    = length(group_order);
n_clusters  = length(sink_order);

% Verify all expected levels are present
for g = 1:n_groups
    if ~any(expcond == group_order(g))
        warning('Group "%s" not found in expcond!', char(group_order(g)));
    end
end
for k = 1:n_clusters
    if ~any(sinktype == sink_order(k))
        warning('Sink type "%s" not found in sinktype!', char(sink_order(k)));
    end
end

% Nice colors: one per GROUP (SS, OS, OD)
% group_colors = [
%     0.22, 0.50, 0.72;   % muted blue   – SS
%     0.85, 0.37, 0.01;   % burnt orange – OS
%     0.13, 0.62, 0.38;   % teal green   – OD
%     0.10, 0.27, 0.08];

% Colors for XX vs XY:

group_colors = groupcolors;

% --- Count matrix: rows = sinktypes, cols = groups (both in fixed order) ---
group_cluster_counts = zeros(n_clusters, n_groups);
for k = 1:n_clusters
    for g = 1:n_groups
        group_cluster_counts(k, g) = sum(sinktype == sink_order(k) & expcond == group_order(g));
    end
end

% Proportion within each sinktype (rows sum to 1)
cluster_totals = sum(group_cluster_counts, 2);
cluster_totals(cluster_totals == 0) = NaN;    % avoid div by zero
cluster_props  = group_cluster_counts ./ cluster_totals;

% X tick labels as plain strings for plotting
sink_labels  = cellstr(sink_order);
group_labels = cellstr(group_order);

% --- Subplot 1: Grouped bar – proportion within each sinktype ---
subplot(2, 2, 1);
b = bar(cluster_props, 'grouped');
for g = 1:n_groups
    b(g).FaceColor = group_colors(g,:);
end
set(gca, 'XTick', 1:n_clusters, 'XTickLabel', sink_labels);
xlabel('Sink Type',                       'FontSize', 11);
ylabel('Proportion within Sink Type',     'FontSize', 11);
title('Group Proportions within Each Sink Type', 'FontSize', 12);
legend(group_labels, 'Location', 'northeastoutside');
max_prop = max(cluster_props(:), [], 'omitnan');
if ~isnan(max_prop) && max_prop > 0
    ylim([0, max_prop * 1.15]);
end

% --- Subplot 2: Stacked bar – proportion ---
subplot(2, 2, 2);
b = bar(cluster_props, 'stacked');
for g = 1:n_groups
    b(g).FaceColor = group_colors(g,:);
end
set(gca, 'XTick', 1:n_clusters, 'XTickLabel', sink_labels);
xlabel('Sink Type',  'FontSize', 11);
ylabel('Proportion', 'FontSize', 11);
title('Group Distribution per Sink Type (Stacked)', 'FontSize', 12);
legend(group_labels, 'Location', 'northeastoutside');
ylim([0, 1]);

% --- Subplot 3: Stacked bar – absolute counts ---
subplot(2, 2, 3);
b = bar(group_cluster_counts, 'stacked');
for g = 1:n_groups
    b(g).FaceColor = group_colors(g,:);
end
set(gca, 'XTick', 1:n_clusters, 'XTickLabel', sink_labels);
xlabel('Sink Type',      'FontSize', 11);
ylabel('Absolute Count', 'FontSize', 11);
title('Absolute Ripple Counts per Sink Type', 'FontSize', 12);
legend(group_labels, 'Location', 'northeastoutside');

% --- Subplot 4: Heatmap of observed counts ---
subplot(2, 2, 4);
imagesc(group_cluster_counts);
colormap(gca, parula);
colorbar;

% Chi-square on the actual per-ripple categorical vectors
[~, chi2_val, p_chi] = crosstab(sinktype, expcond);

set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sink_labels, ...
         'XTick', 1:n_groups,   'XTickLabel', group_labels);
ylabel('Sink Type', 'FontSize', 11);
xlabel('Group',     'FontSize', 11);
title(sprintf('Observed Counts (χ² p=%.4f)', p_chi), 'FontSize', 12);

% Overlay count values
for k = 1:n_clusters
    for g = 1:n_groups
        text(g, k, num2str(group_cluster_counts(k, g)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, ...
            'Color', 'w', 'FontWeight', 'bold');
    end
end

sgtitle('Group Distribution per Sink Type', 'FontSize', 14, 'FontWeight', 'bold');
savepath = pwd;
saveas(fig3, fullfile(savepath, 'Fig3_Group_Cluster_Distribution.svg'));
saveas(fig3, fullfile(savepath, 'Fig3_Group_Cluster_Distribution.jpg'));


%% Plot 3b: Cluster proportions within each group (pie charts)
fig3b = figure('Position', [100, 100, 1400, 550]);

% Colors: one per sink type (rad, sink, other) — consistent with fig3
sink_colors = [
    0.22, 0.50, 0.72;   % muted blue   – rad
    0.85, 0.37, 0.01;   % burnt orange – sink
    0.13, 0.62, 0.38;   % teal green   – other
];

for g = 1:n_groups
    subplot(1, n_groups, g);

    % Proportion of each sinktype within this group
    % group_cluster_counts: rows=sinktypes, cols=groups
    pie_data = group_cluster_counts(:, g);          % absolute counts, n_clusters x 1
    pie_data_safe = pie_data;
    pie_data_safe(pie_data == 0) = [];              % pie() crashes on zero slices

    % Track which sinktypes are non-zero for correct color assignment
    nonzero_mask = pie_data > 0;

    h = pie(pie_data_safe);

    % Apply consistent colors per sink type
    patch_idx = 1;
    color_idx  = find(nonzero_mask);
    for k = 1:length(color_idx)
        h(patch_idx).FaceColor = sink_colors(color_idx(k), :);
        h(patch_idx).EdgeColor = 'w';
        h(patch_idx+1).FontSize = 10;               % percentage label
        patch_idx = patch_idx + 2;
    end

    % Group total = sum across sinktypes for this group
    group_n = sum(group_cluster_counts(:, g));
    title(sprintf('%s  (n=%d)', char(group_order(g)), group_n), ...
        'FontSize', 13, 'FontWeight', 'bold');
end

% Single shared legend at the bottom, one entry per sink type
lgd = legend(sink_labels, ...
    'Orientation', 'horizontal', 'FontSize', 11);
lgd.Position(1) = 0.5 - lgd.Position(3)/2;   % center horizontally
lgd.Position(2) = 0.02;                        % near bottom

sgtitle('Sink Type Proportions within Each Group', 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig3b, fullfile(savepath, 'Fig3b_Cluster_Proportions_Pretty.svg'));
saveas(fig3b, fullfile(savepath, 'Fig3b_Cluster_Proportions_Pretty.jpg'));


%%

%% Plot 4: Representative CSD Patterns per Sink Type (using all_csd_allch)
all_csd_allch = all_csd;
Fs = 1250;
fig4 = figure('Units', 'normalized', ...
              'Position', [100, 100, 1600, 400]);

[n_time_full, n_channels_all, ~] = size(all_csd_allch);
time_vec_full    = ((1:n_time_full) - ceil(n_time_full/2)) / Fs * 1000;
depth_vec_all    = 1:n_channels_all;
global_csd_limit = prctile(abs(all_csd_allch(:)), 95);

% Red=positive (source), White=zero, Blue=negative (sink)
n_colors = 256;
half     = n_colors / 2;
rwb      = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1);   % blue -> white
            ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];   % white -> red

for k = 1:n_clusters

    cluster_idx = find(sinktype == sink_order(k));
    avg_pattern = mean(all_csd_allch(:, :, cluster_idx), 3)';

    subplot(1, n_clusters, k);
    contourf(time_vec_full, depth_vec_all, avg_pattern, 40, 'LineColor', 'none');
    colormap(gca, rwb);
    colorbar;
    set(gca, 'YDir', 'reverse');
    ylim([1, n_channels_all]);
    yticks(1:n_channels_all);
    caxis([-global_csd_limit, global_csd_limit]);

    xlabel('Time (ms)', 'FontSize', 11);
    ylabel('Channel',   'FontSize', 11);
    title(sprintf('%s  (n=%d)', char(sink_order(k)), length(cluster_idx)), ...
        'FontSize', 12, 'FontWeight', 'bold');

    if exist('lfpwaveform', 'var') && ~isempty(cluster_idx)
        mean_lfp = squeeze(mean(lfpwaveform(:, :, cluster_idx), 3));
        n_lfp_t  = size(mean_lfp, 1);
        ch_lfp   = size(mean_lfp, 2);

        if n_lfp_t ~= n_time_full
            t_lfp    = linspace(time_vec_full(1), time_vec_full(end), n_lfp_t);
            mean_lfp = interp1(t_lfp, mean_lfp, time_vec_full);
        end
        if ch_lfp ~= n_channels_all
            mean_lfp = interp1(1:ch_lfp, mean_lfp', ...
                linspace(1, ch_lfp, n_channels_all), 'linear')';
        end

        mean_lfp = reshape(mean_lfp, n_time_full, n_channels_all);
        gmax     = max(abs(mean_lfp(:))) + eps;
        scale    = (max(depth_vec_all) - min(depth_vec_all)) / 4;
        overlay  = repmat(depth_vec_all(:)', n_time_full, 1) + (mean_lfp ./ gmax) * scale;

        pad      = 0.03 * (max(depth_vec_all) - min(depth_vec_all));
        new_ylim = [min(1, min(overlay(:))) - pad, ...
                    max(n_channels_all, max(overlay(:))) + pad];
        hold on;
        for ch = 1:n_channels_all
            plot(time_vec_full, overlay(:, ch), 'k-', 'LineWidth', 1);
        end
        hold off;
        ylim(new_ylim);
    end
end

sgtitle('CSD Patterns per Sink Type', 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig4, fullfile(savepath, 'Fig4_CSD_Patterns_AllChannels.svg'));
saveas(fig4, fullfile(savepath, 'Fig4_CSD_Patterns_AllChannels.jpg'));

%% Shared colormap (red=positive/source, blue=negative/sink)
n_colors = 256; half = n_colors/2;
rwb = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1);
       ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];

%% Plot 5: CSD Patterns + Waveforms per Sink Type
% lfpwaveform = allLFPwaveforms;
Fs = 1250;
if exist('filteredripples', 'var') && exist('rawripples', 'var')
    fig5 = figure('Position', [100, 100, 1800, 900]);

    if istable(filteredripples); filteredripples = table2array(filteredripples); end
    if istable(rawripples);      rawripples      = table2array(rawripples);      end

    n_samples_wave_filt = size(filteredripples, 2);
    n_samples_wave_raw  = size(rawripples, 2);
    time_vec_wave_filt  = ((1:n_samples_wave_filt) - ceil(n_samples_wave_filt/2)) / Fs * 1000;
    time_vec_wave_raw   = ((1:n_samples_wave_raw)  - ceil(n_samples_wave_raw/2))  / Fs * 1000;

    all_filt_avgs = zeros(n_clusters, n_samples_wave_filt);
    all_raw_avgs  = zeros(n_clusters, n_samples_wave_raw);
    for k = 1:n_clusters
        cluster_idx = find(sinktype == sink_order(k));   % FIX: categorical
        all_filt_avgs(k,:) = mean(filteredripples(cluster_idx,:), 1);
        all_raw_avgs(k,:)  = mean(rawripples(cluster_idx,:),      1);
    end
    ylim_filt = [min(all_filt_avgs(:)), max(all_filt_avgs(:))];
    ylim_raw  = [min(all_raw_avgs(:)),  max(all_raw_avgs(:))];

    for k = 1:n_clusters
        cluster_idx = find(sinktype == sink_order(k));   % FIX: categorical
        sink_label  = char(sink_order(k));

        % Row 1: CSD
        subplot(3, n_clusters, k);
        avg_pattern = mean(all_csd_allch(:,:,cluster_idx), 3)';
        contourf(time_vec_full, depth_vec_all, avg_pattern, 40, 'LineColor', 'none');
        colormap(gca, rwb); colorbar;
        set(gca, 'YDir', 'reverse'); ylim([1, n_channels_all]); yticks(1:n_channels_all);
        caxis([-global_csd_limit, global_csd_limit]);
        ylabel('Channel');
        title(sprintf('%s CSD (n=%d)', sink_label, length(cluster_idx)));  % FIX: label

        if exist('lfpwaveform','var') && ~isempty(cluster_idx)
            mean_lfp = squeeze(mean(lfpwaveform(:,:,cluster_idx), 3));
            n_lfp_t = size(mean_lfp,1); ch_lfp = size(mean_lfp,2);
            if n_lfp_t ~= n_time_full
                t_lfp = linspace(time_vec_full(1), time_vec_full(end), n_lfp_t);
                mean_lfp = interp1(t_lfp, mean_lfp, time_vec_full);
            end
            if ch_lfp ~= n_channels_all
                mean_lfp = interp1(1:ch_lfp, mean_lfp', linspace(1,ch_lfp,n_channels_all), 'linear')';
            end
            mean_lfp = reshape(mean_lfp, n_time_full, n_channels_all);
            gmax  = max(abs(mean_lfp(:))) + eps;
            scale = (max(depth_vec_all) - min(depth_vec_all)) / 4;
            overlay = repmat(depth_vec_all(:)', n_time_full, 1) + (mean_lfp ./ gmax) * scale;
            pad = 0.03 * (max(depth_vec_all) - min(depth_vec_all));
            hold on;
            for ch = 1:n_channels_all; plot(time_vec_full, overlay(:,ch), 'k-', 'LineWidth', 1); end
            hold off;
            ylim([min(1,min(overlay(:)))-pad, max(n_channels_all,max(overlay(:)))+pad]);
        end

        % Row 2: filtered waveform
        subplot(3, n_clusters, n_clusters+k);
        plot(time_vec_wave_filt, all_filt_avgs(k,:), 'k-', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Filtered LFP'); title('Mean Filtered Waveform');
        xlim([time_vec_wave_filt(1) time_vec_wave_filt(end)]); ylim(ylim_filt);

        % Row 3: raw waveform
        subplot(3, n_clusters, 2*n_clusters+k);
        plot(time_vec_wave_raw, all_raw_avgs(k,:), 'k-', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Raw LFP'); title('Mean Raw Waveform');
        xlim([time_vec_wave_raw(1) time_vec_wave_raw(end)]); ylim(ylim_raw);
    end

    saveas(fig5, fullfile(savepath, 'Fig5_CSD_Waveforms.svg'));
    saveas(fig5, fullfile(savepath, 'Fig5_CSD_Waveforms.jpg'));
end

%% Plot 6: CSD Patterns per Group
fig6 = figure('Position', [100, 100, 1400, 400]);
for g = 1:n_groups
    group_idx = find(expcond == group_order(g));   % FIX: categorical, not integer index
    avg_pattern = mean(all_csd_allch(:,:,group_idx), 3)';

    subplot(1, n_groups, g);
    contourf(time_vec_full, depth_vec_all, avg_pattern, 40, 'LineColor', 'none');
    colormap(gca, rwb); colorbar;
    set(gca, 'YDir', 'reverse'); ylim([1, n_channels_all]); yticks(1:n_channels_all);
    caxis([-global_csd_limit, global_csd_limit]);
    xlabel('Time (ms)'); ylabel('Channel');
    title(sprintf('%s (n=%d)', char(group_order(g)), length(group_idx)));  % FIX: label

    if exist('lfpwaveform','var') && ~isempty(group_idx)
        mean_lfp = squeeze(mean(lfpwaveform(:,:,group_idx), 3));
        n_lfp_t = size(mean_lfp,1); ch_lfp = size(mean_lfp,2);
        if n_lfp_t ~= n_time_full
            t_lfp = linspace(time_vec_full(1), time_vec_full(end), n_lfp_t);
            mean_lfp = interp1(t_lfp, mean_lfp, time_vec_full);
        end
        if ch_lfp ~= n_channels_all
            mean_lfp = interp1(1:ch_lfp, mean_lfp', linspace(1,ch_lfp,n_channels_all), 'linear')';
        end
        mean_lfp = reshape(mean_lfp, n_time_full, n_channels_all);
        gmax  = max(abs(mean_lfp(:))) + eps;
        scale = (max(depth_vec_all) - min(depth_vec_all)) / 4;
        overlay = repmat(depth_vec_all(:)', n_time_full, 1) + (mean_lfp ./ gmax) * scale;
        pad = 0.03 * (max(depth_vec_all) - min(depth_vec_all));
        hold on;
        for ch = 1:n_channels_all; plot(time_vec_full, overlay(:,ch), 'k-', 'LineWidth', 1); end
        hold off;
        ylim([min(1,min(overlay(:)))-pad, max(n_channels_all,max(overlay(:)))+pad]);
    end
end

saveas(fig6, fullfile(savepath, 'Fig6_CSD_Patterns_ByGroup.svg'));
saveas(fig6, fullfile(savepath, 'Fig6_CSD_Patterns_ByGroup.jpg'));

%% Rearrange 4 layers data

% LOAD CSD DATA WITH 4 COLUMNS: ONE FOR EACH LAYER 

% Average over mid 50ms (-25ms to +25ms) centered at peak (assuming uniform sampling)
nT = size(all_csd_layers, 1);
mid_idx = round(nT/2 - nT*25/200) : round(nT/2 + nT*25/200 - 1);
mean_csd = squeeze(mean(all_csd_layers(mid_idx,:,:), 1))';

% Layer names in order: or, pyr, rad, slm
layer_names = {'or', 'pyr', 'rad', 'slm'};
nripples = size(mean_csd, 1);

% Replicate metadata 4x (once per layer) and stack layer data
finaltable = table(...
    repmat(all_ripples.subject, 4, 1), ...
    repmat((1:nripples)',        4, 1), ...
    repmat(all_ripples.ID,      4, 1), ...
    repmat(all_ripples.expcond, 4, 1), ...
    repmat(all_ripples.sinktype,4, 1), ...
    categorical(repelem(layer_names', nripples)), ...
    vertcat(mean_csd(:,1), mean_csd(:,2), mean_csd(:,3), mean_csd(:,4)), ...
    'VariableNames', {'subject','rippleID','ID','expcond','sinktype','layer','CSD'});

%% STATS

% Pre-aggregate to subject x layer means before fitting

radtable = finaltable(finaltable.sinktype == "rad", :);
slmtable = finaltable(finaltable.sinktype == "slm", :);
othertable = finaltable(finaltable.sinktype == "other", :);

radmeanstbl = groupsummary(radtable, {'subject','ID','expcond','layer'}, 'mean', 'CSD');
radmeanstbl.Properties.VariableNames{'mean_CSD'} = 'CSD';
results_rad = analyze_csd_lme_2way(radmeanstbl, 'CSD', 'expcond', ["SS", "OS", "OD"]);

slmmeanstbl = groupsummary(slmtable, {'subject','ID','expcond','layer'}, 'mean', 'CSD');
slmmeanstbl.Properties.VariableNames{'mean_CSD'} = 'CSD';
results_slm = analyze_csd_lme_2way(slmmeanstbl, 'CSD', 'expcond',  ["SS", "OS", "OD"]);

othermeanstbl = groupsummary(othertable, {'subject','ID','expcond','layer'}, 'mean', 'CSD');
othermeanstbl.Properties.VariableNames{'mean_CSD'} = 'CSD';
results_other = analyze_csd_lme_2way(othermeanstbl, 'CSD', 'expcond',  ["SS", "OS", "OD"]);


% results = analyze_csd_lme_2way(finaltable, 'CSD', 'expcond'); >>> running
% this is too slow, better to summarize data

%%

%% Plot 9: CSD Patterns per Sink Type × Group
% Requires: all_csd_allch, sinktype, expcond (from main script)
% sink_order, group_order, n_clusters, n_groups, time_vec_full,
% depth_vec_all, n_channels_all, n_time_full, global_csd_limit, savepath

% Define custom colormaps
redwhiteblue = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
                ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];
bluewhitered = flipud(redwhiteblue);

fig9 = figure('Units', 'normalized', ...
              'Position', [0.02, 0.05, 0.96, 0.88], ...
              'Color', 'w');

n_groups = length(group_labels);
n_clusters = 3;

for g = 1:n_groups
    for k = 1:n_clusters

        % --- Find ripples belonging to this group × sink type combination ---
        combined_idx = find(sinktype == sink_order(k) & expcond == group_order(g));

        ax = subplot(n_groups, n_clusters, (g-1)*n_clusters + k);

        if ~isempty(combined_idx)

            avg_pattern = mean(all_csd_allch(:, :, combined_idx), 3)';

            contourf(time_vec_full, depth_vec_all, avg_pattern, 40, 'LineColor', 'none');
            colormap(ax, redwhiteblue);
            colorbar(ax);
            set(ax, 'YDir', 'reverse');
            ylim(ax, [1, n_channels_all]);
            yticks(ax, 1:n_channels_all);
            caxis(ax, [-global_csd_limit, global_csd_limit]);

            % --- Optional LFP overlay ---
            if exist('lfpwaveform', 'var') && ~isempty(lfpwaveform)

                mean_lfp = squeeze(mean(lfpwaveform(:, :, combined_idx), 3));
                n_lfp_t  = size(mean_lfp, 1);
                ch_lfp   = size(mean_lfp, 2);

                % Resample time axis if needed
                if n_lfp_t ~= n_time_full
                    t_lfp    = linspace(time_vec_full(1), time_vec_full(end), n_lfp_t);
                    mean_lfp = interp1(t_lfp, mean_lfp, time_vec_full);
                end

                % Resample channels if needed
                if ch_lfp ~= n_channels_all
                    mean_lfp = interp1(1:ch_lfp, mean_lfp', ...
                                       linspace(1, ch_lfp, n_channels_all), 'linear')';
                end

                mean_lfp = reshape(mean_lfp, n_time_full, n_channels_all);
                gmax     = max(abs(mean_lfp(:))) + eps;
                scale    = (max(depth_vec_all) - min(depth_vec_all)) / 4;
                overlay  = repmat(depth_vec_all(:)', n_time_full, 1) + ...
                           (mean_lfp ./ gmax) * scale;

                pad     = 0.03 * (max(depth_vec_all) - min(depth_vec_all));
                new_ylim = [min(min(depth_vec_all), min(overlay(:))) - pad, ...
                             max(max(depth_vec_all), max(overlay(:))) + pad];

                hold(ax, 'on');
                for ch = 1:n_channels_all
                    plot(ax, time_vec_full, overlay(:, ch), 'k-', 'LineWidth', 0.8);
                end
                hold(ax, 'off');
                ylim(ax, new_ylim);
            end

        else
            % Empty panel
            text(0.5, 0.5, 'No data', ...
                 'HorizontalAlignment', 'center', ...
                 'Units', 'normalized', ...
                 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
            set(ax, 'XTick', [], 'YTick', []);
        end

        % --- Labels: only on border subplots to reduce clutter ---
        % Column headers (sink type): top row only
        if g == 1
            title(ax, char(sink_order(k)), 'FontSize', 12, 'FontWeight', 'bold');
        end

        % Row labels (group): left column only
        if k == 1
            ylabel(ax, char(group_order(g)), 'FontSize', 12, 'FontWeight', 'bold');
        else
            ylabel(ax, '');
        end

        % X-axis label: bottom row only
        if g == n_groups
            xlabel(ax, 'Time (ms)', 'FontSize', 10);
        else
            xlabel(ax, '');
        end

        set(ax, 'FontSize', 9, 'Box', 'off');
    end
end

sgtitle('CSD Patterns per Sink Type × Group', 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig9, fullfile(savepath, 'Fig9_CSD_SinkType_by_Group.svg'));
saveas(fig9, fullfile(savepath, 'Fig9_CSD_SinkType_by_Group.jpg'));

%% =========================================================================
%  SINK TYPE PROPORTION ANALYSIS — HIERARCHICAL VERSION
%  Question: Does the proportion of rad / slm / other sink-type ripples
%            differ between experimental groups?
%
%  Statistical unit : recording session (ID)
%  Random effect    : subject (animal), to account for multiple sessions
%                     per animal
%  Model            : logit(proportion) ~ group + (1|subject)
%                     fitted with fitglme using binomial distribution, or
%                     fitlme on logit-transformed proportions
% =========================================================================

%% ── 0. Setup ─────────────────────────────────────────────────────────────
% Set these explicitly for the current experiment
group_order  = categorical({'XX','XY'});   % <-- change per experiment
sink_order   = categorical({'rad','slm','other'});
n_groups     = numel(group_order);
n_sinktypes  = numel(sink_order);
sink_vars    = {'prop_rad','prop_slm','prop_other'};
sink_names   = cellstr(sink_order);
N_PERM       = 5000;    % reduced from 10000 — LME permutations are slower
EPSILON      = 0.001;
ALPHA        = 0.05;

% Helper functions
p_to_star = @(p) repmat('*', 1, (p<0.05) + (p<0.01) + (p<0.001));
logit_fn  = @(p) log((p + EPSILON) ./ (1 - p + EPSILON));

% Sanity check
fprintf('\n=== DATA CHECK ===\n');
fprintf('%-14s  %8s  %8s  %10s\n','Group','Sessions','Subjects','Ripples');
fprintf('%s\n', repmat('-',1,44));
for g = 1:n_groups
    grp_mask = all_ripples.expcond == group_order(g);
    fprintf('%-14s  %8d  %8d  %10d\n', char(group_order(g)), ...
        numel(unique(all_ripples.ID(grp_mask))), ...
        numel(unique(all_ripples.subject(grp_mask))), ...
        sum(grp_mask));
end
fprintf('\n');


%% ── 1. Build per-session proportion table ────────────────────────────────
sessions   = unique(all_ripples.ID);
n_sessions = numel(sessions);

prop_table = table('Size', [n_sessions, 4 + n_sinktypes], ...
    'VariableTypes', [{'categorical','categorical','categorical','double'}, ...
                       repmat({'double'}, 1, n_sinktypes)], ...
    'VariableNames', [{'session','subject','group','n_ripples'}, sink_vars]);

valid_rows = 0;
for s = 1:n_sessions
    sess      = sessions(s);
    sess_mask = all_ripples.ID == sess;
    sess_data = all_ripples(sess_mask, :);
    total     = height(sess_data);
    if total == 0, continue; end

    valid_rows = valid_rows + 1;
    prop_table.session(valid_rows)   = sess;
    prop_table.subject(valid_rows)   = sess_data.subject(1);
    prop_table.group(valid_rows)     = sess_data.expcond(1);
    prop_table.n_ripples(valid_rows) = total;

    for k = 1:n_sinktypes
        prop_table.(sink_vars{k})(valid_rows) = ...
            sum(sess_data.sinktype == sink_order(k)) / total;
    end
end
prop_table = prop_table(1:valid_rows, :);

% Print session/subject breakdown
fprintf('=== Per-session proportion table (%d sessions, %d subjects) ===\n', ...
    valid_rows, numel(unique(prop_table.subject)));
fprintf('Sessions per subject:\n');
subj_list = unique(prop_table.subject);
for s = 1:numel(subj_list)
    subj_mask = prop_table.subject == subj_list(s);
    fprintf('  %s (%s): %d sessions\n', ...
        char(subj_list(s)), ...
        char(prop_table.group(find(subj_mask,1))), ...
        sum(subj_mask));
end
fprintf('\n');
disp(prop_table);


%% ── 2. Logit-transform ───────────────────────────────────────────────────
logit_vars = {'logit_rad','logit_slm','logit_other'};
for k = 1:n_sinktypes
    prop_table.(logit_vars{k}) = logit_fn(prop_table.(sink_vars{k}));
end


%% ── 3. Linear Mixed Model: logit(prop) ~ group + (1|subject) ────────────
%
%  This is the primary statistical test. The random intercept for subject
%  absorbs between-animal variance, so the fixed effect of group is tested
%  only on within- and between-subject variation that is not explained by
%  individual differences.
%
%  For 2 groups: the group coefficient t-test from fitlme is the result.
%  For 3+ groups: we compare the full model vs an intercept-only model
%  using a likelihood ratio test to get an omnibus p-value, then run
%  pairwise contrasts for post-hoc comparisons.
%
fprintf('\n=== LINEAR MIXED MODEL: logit(prop) ~ group + (1|subject) ===\n');

lme_results = struct('model',[],'stat',[],'df',[],'p_raw',[],'p_fdr',[], ...
                     'posthoc',[],'posthoc_names',[],'coefTable',[]);
lme_results = repmat(lme_results, 1, n_sinktypes);

for k = 1:n_sinktypes
    sname = sink_names{k};
    lvar  = logit_vars{k};

    % Fit full model and null model (intercept + random effect only)
    formula_full = sprintf('%s ~ group + (1|subject)', lvar);
    formula_null = sprintf('%s ~ 1 + (1|subject)', lvar);

    lme_full = fitlme(prop_table, formula_full, ...
        'FitMethod','REML', 'DummyVarCoding','reference');
    lme_null = fitlme(prop_table, formula_null, ...
        'FitMethod','REML', 'DummyVarCoding','reference');

    lme_results(k).model     = lme_full;
    lme_results(k).coefTable = lme_full.Coefficients;

    fprintf('\n--- %s ---\n', sname);
    fprintf('Random effect variance (subject): %.4f\n', ...
        lme_full.covarianceParameters{1}(1));
    fprintf('Residual variance               : %.4f\n\n', ...
        lme_full.MSE);
    disp(lme_full.Coefficients);

    if n_groups == 2
        % For 2 groups: use the t-test on the group coefficient directly
        coef_idx = 2;   % row 1 = intercept, row 2 = group effect
        t_val    = lme_full.Coefficients.tStat(coef_idx);
        df_val   = lme_full.Coefficients.DF(coef_idx);
        p_val    = lme_full.Coefficients.pValue(coef_idx);

        lme_results(k).stat   = t_val;
        lme_results(k).df     = df_val;
        lme_results(k).p_raw  = p_val;

        fprintf('%s  |  LME t(%.1f)=%.3f,  p=%.4f\n', sname, df_val, t_val, p_val);

    else
        % For 3+ groups: likelihood ratio test (full vs null)
        lr_results = compare(lme_null, lme_full, 'CheckNesting', true);
        chi2_val   = lr_results.LRStat(2);
        df_lr      = lr_results.deltaDF(2);
        p_lr       = lr_results.pValue(2);

        lme_results(k).stat  = chi2_val;
        lme_results(k).df    = df_lr;
        lme_results(k).p_raw = p_lr;

        fprintf('%s  |  LRT chi2(%d)=%.3f,  p=%.4f\n', sname, df_lr, chi2_val, p_lr);

        % Post-hoc pairwise contrasts using coefTest
        fprintf('  Post-hoc pairwise contrasts (Bonferroni-corrected):\n');
        coef_names = lme_full.CoefficientNames;
        group_coef_idx = find(~strcmp(coef_names, '(Intercept)'));
        group_levels   = group_order;   % all group levels
        ref_group      = group_levels(1);

        % Build all pairs
        pairs     = nchoosek(1:n_groups, 2);
        n_pairs   = size(pairs, 1);
        posthoc   = zeros(n_pairs, 4);   % [g1 g2 t p]
        ph_names  = cell(n_pairs, 2);

        for p_idx = 1:n_pairs
            g1 = pairs(p_idx, 1);
            g2 = pairs(p_idx, 2);
            ph_names{p_idx,1} = char(group_levels(g1));
            ph_names{p_idx,2} = char(group_levels(g2));

            % Build contrast vector
            % Reference group (g=1) has no explicit coefficient — it is
            % captured by the intercept. All other groups have a coefficient.
            c = zeros(1, numel(coef_names));
            if g1 == 1
                % contrast: 0 - coef(g2)
                coef_g2 = find(strcmp(coef_names, ...
                    sprintf('group_%s', char(group_levels(g2)))));
                if isempty(coef_g2)
                    % MATLAB names dummies as 'group_XY' etc — try alternative
                    coef_g2 = group_coef_idx(g2-1);
                end
                c(coef_g2) = -1;
            elseif g2 == 1
                coef_g1 = find(strcmp(coef_names, ...
                    sprintf('group_%s', char(group_levels(g1)))));
                if isempty(coef_g1)
                    coef_g1 = group_coef_idx(g1-1);
                end
                c(coef_g1) = 1;
            else
                % Both non-reference: contrast is coef(g1) - coef(g2)
                coef_g1 = group_coef_idx(g1-1);
                coef_g2 = group_coef_idx(g2-1);
                c(coef_g1) =  1;
                c(coef_g2) = -1;
            end

            [p_c, F_c, df1_c, df2_c] = coefTest(lme_full, c);
            t_c = sign(c * lme_full.Coefficients.Estimate) * sqrt(F_c);
            posthoc(p_idx,:) = [g1, g2, t_c, p_c];
            fprintf('    %-12s vs %-12s  t=%.3f  p=%.4f  (Bonf p=%.4f)\n', ...
                ph_names{p_idx,1}, ph_names{p_idx,2}, t_c, p_c, ...
                min(p_c * n_pairs, 1));
        end

        lme_results(k).posthoc       = posthoc;
        lme_results(k).posthoc_names = ph_names;
    end
end


%% ── 4. Permutation test (subject-aware, fast) ────────────────────────────
%
%  Statistic: F from one-way ANOVA on per-subject means (one value per
%             animal). This respects the nested structure without needing
%             to refit the LME, and matches what the LME is testing.
%  Permutation: group labels shuffled at subject level.
%  N_PERM = 10000 runs in seconds with this approach.
%
N_PERM = 10000;

fprintf('\n=== SUBJECT-LEVEL PERMUTATION TEST (%d permutations) ===\n', N_PERM);

perm_results = struct('obs_stat',[],'null_dist',[],'p_raw',[],'p_fdr',[]);
perm_results = repmat(perm_results, 1, n_sinktypes);

% Build subject-level summary table (one row per subject)
subj_list   = unique(prop_table.subject);
n_subj      = numel(subj_list);
subj_grp    = arrayfun(@(s) prop_table.group(find(prop_table.subject==s,1)), ...
                        subj_list);   % n_subj x 1 categorical

rng(42);

for k = 1:n_sinktypes
    lvar = logit_vars{k};

    % Compute per-subject mean of logit-transformed proportion
    subj_means = zeros(n_subj, 1);
    for s = 1:n_subj
        subj_mask     = prop_table.subject == subj_list(s);
        subj_means(s) = mean(prop_table.(lvar)(subj_mask));
    end

    % Observed F on subject means
    [~, tbl_obs, ~] = anova1(subj_means, subj_grp, 'off');
    obs_stat = tbl_obs{2,5};

    % Permutation null distribution
    null_dist = zeros(N_PERM, 1);
    for i = 1:N_PERM
        perm_grp      = subj_grp(randperm(n_subj));
        [~, tbl_p, ~] = anova1(subj_means, perm_grp, 'off');
        null_dist(i)  = tbl_p{2,5};
    end

    p_perm = mean(null_dist >= obs_stat);
    perm_results(k).obs_stat  = obs_stat;
    perm_results(k).null_dist = null_dist;
    perm_results(k).p_raw     = p_perm;

    fprintf('%s  |  obs F=%.4f,  p=%.4f\n', sink_names{k}, obs_stat, p_perm);
end

%% ── 5. FDR correction ────────────────────────────────────────────────────
p_lme  = zeros(1, n_sinktypes);
p_perm = zeros(1, n_sinktypes);
for k = 1:n_sinktypes
    p_lme(k)  = lme_results(k).p_raw;
    p_perm(k) = perm_results(k).p_raw;
end

p_lme_fdr  = deal_bh(p_lme);
p_perm_fdr = deal_bh(p_perm);

fprintf('\n=== FDR-CORRECTED p-VALUES ===\n');
fprintf('%-10s  %10s  %10s  %10s  %10s\n', ...
    'SinkType','p_LME','p_LME_FDR','p_perm','p_perm_FDR');
fprintf('%s\n', repmat('-',1,55));
for k = 1:n_sinktypes
    fprintf('%-10s  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
        sink_names{k}, p_lme(k), p_lme_fdr(k), p_perm(k), p_perm_fdr(k));
    lme_results(k).p_fdr  = p_lme_fdr(k);
    perm_results(k).p_fdr = p_perm_fdr(k);
end


%% ── 6. Summary table ─────────────────────────────────────────────────────
stat_label = 'LME_t';
if n_groups > 2, stat_label = 'LRT_chi2'; end

lme_stats   = zeros(n_sinktypes,1);
lme_p_raw   = zeros(n_sinktypes,1);
lme_p_fdr   = zeros(n_sinktypes,1);
perm_stats  = zeros(n_sinktypes,1);
perm_p_raw  = zeros(n_sinktypes,1);
perm_p_fdr  = zeros(n_sinktypes,1);

for k = 1:n_sinktypes
    lme_stats(k)  = lme_results(k).stat;
    lme_p_raw(k)  = lme_results(k).p_raw;
    lme_p_fdr(k)  = lme_results(k).p_fdr;
    perm_stats(k) = perm_results(k).obs_stat;
    perm_p_raw(k) = perm_results(k).p_raw;
    perm_p_fdr(k) = perm_results(k).p_fdr;
end

summary_table = table(sink_names(:), lme_stats, lme_p_raw, lme_p_fdr, ...
    perm_stats, perm_p_raw, perm_p_fdr, ...
    'VariableNames', {'SinkType', stat_label, 'p_LME', 'p_LME_FDR', ...
                      'perm_F', 'p_permutation', 'p_perm_FDR'});

fprintf('\n=== SUMMARY TABLE ===\n');
disp(summary_table);


%% ── 6b. Save statistical report ─────────────────────────────────────────
savepath = pwd;
report_path = fullfile(savepath, 'SinkType_Proportion_Stats.txt');
fid = fopen(report_path, 'w');

fprintf(fid, '=========================================================\n');
fprintf(fid, '  SINK TYPE PROPORTION ANALYSIS — STATISTICAL REPORT\n');
fprintf(fid, '  Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '=========================================================\n\n');

fprintf(fid, '--- EXPERIMENTAL DESIGN ---\n');
fprintf(fid, 'Statistical unit  : recording session (ID)\n');
fprintf(fid, 'Hierarchy         : sessions nested within subjects (animals)\n');
fprintf(fid, 'Primary test      : LME  logit(prop) ~ group + (1|subject)\n');
fprintf(fid, 'Backup test       : subject-level permutation test\n');
fprintf(fid, 'Groups            : %s\n', strjoin(cellstr(group_order), ', '));
fprintf(fid, 'Sink types        : %s\n', strjoin(sink_names, ', '));
fprintf(fid, 'Total sessions    : %d\n', valid_rows);
fprintf(fid, 'Total subjects    : %d\n', numel(unique(prop_table.subject)));
fprintf(fid, 'Total ripples     : %d\n\n', sum(prop_table.n_ripples));

fprintf(fid, '--- SAMPLE SIZES ---\n');
fprintf(fid, '%-14s  %8s  %8s  %10s\n','Group','Sessions','Subjects','Ripples');
fprintf(fid, '%s\n', repmat('-',1,44));
for g = 1:n_groups
    grp_mask = prop_table.group == group_order(g);
    fprintf(fid, '%-14s  %8d  %8d  %10d\n', char(group_order(g)), ...
        sum(grp_mask), ...
        numel(unique(prop_table.subject(grp_mask))), ...
        sum(prop_table.n_ripples(grp_mask)));
end
fprintf(fid, '\n');

fprintf(fid, '--- SESSIONS PER SUBJECT ---\n');
for s = 1:numel(subj_list)
    subj_mask = prop_table.subject == subj_list(s);
    fprintf(fid, '  %-10s  group=%-12s  %d sessions\n', ...
        char(subj_list(s)), ...
        char(prop_table.group(find(subj_mask,1))), ...
        sum(subj_mask));
end
fprintf(fid, '\n');

fprintf(fid, '--- DESCRIPTIVE STATISTICS (proportions per session) ---\n');
for k = 1:n_sinktypes
    fprintf(fid, '\nSink type: %s\n', sink_names{k});
    fprintf(fid, '  %-14s  %6s  %8s  %8s  %8s  %8s\n', ...
        'Group','N','Mean','SD','SEM','Median');
    fprintf(fid, '  %s\n', repmat('-',1,58));
    for g = 1:n_groups
        y = prop_table.(sink_vars{k})(prop_table.group == group_order(g));
        if isempty(y), continue; end
        fprintf(fid, '  %-14s  %6d  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
            char(group_order(g)), numel(y), mean(y), std(y), ...
            std(y)/sqrt(numel(y)), median(y));
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- LME RESULTS: logit(prop) ~ group + (1|subject) ---\n\n');
for k = 1:n_sinktypes
    fprintf(fid, 'Sink type: %s\n', sink_names{k});
    fprintf(fid, '  Random effect variance (subject): %.4f\n', ...
        lme_results(k).model.covarianceParameters{1}(1));
    fprintf(fid, '  Residual variance               : %.4f\n', ...
        lme_results(k).model.MSE);
    if n_groups == 2
        fprintf(fid, '  t(%.1f) = %.4f,  p_raw = %.4f,  p_FDR = %.4f  %s\n\n', ...
            lme_results(k).df, lme_results(k).stat, ...
            lme_results(k).p_raw, lme_results(k).p_fdr, ...
            yesno(lme_results(k).p_fdr < ALPHA));
    else
        fprintf(fid, '  LRT chi2(%d) = %.4f,  p_raw = %.4f,  p_FDR = %.4f  %s\n', ...
            lme_results(k).df, lme_results(k).stat, ...
            lme_results(k).p_raw, lme_results(k).p_fdr, ...
            yesno(lme_results(k).p_fdr < ALPHA));
        if ~isempty(lme_results(k).posthoc)
            fprintf(fid, '  Post-hoc pairwise contrasts (Bonferroni):\n');
            ph = lme_results(k).posthoc;
            phn = lme_results(k).posthoc_names;
            n_pairs = size(ph,1);
            for r = 1:n_pairs
                fprintf(fid, '    %-12s vs %-12s  t=%.3f  p=%.4f  Bonf_p=%.4f  %s\n', ...
                    phn{r,1}, phn{r,2}, ph(r,3), ph(r,4), ...
                    min(ph(r,4)*n_pairs,1), yesno(ph(r,4)*n_pairs < ALPHA));
            end
        end
        fprintf(fid, '\n');
    end
end

fprintf(fid, '--- SUBJECT-LEVEL PERMUTATION TEST ---\n');
fprintf(fid, 'Permutations : %d (seed = 42)\n', N_PERM);
fprintf(fid, 'Statistic    : F (ANOVA on logit-transformed proportions)\n');
fprintf(fid, 'Labels permuted at subject level to preserve session nesting\n\n');
fprintf(fid, '%-10s  %10s  %10s  %10s  %15s\n', ...
    'SinkType','obs_F','p_raw','p_FDR','Significant');
fprintf(fid, '%s\n', repmat('-',1,60));
for k = 1:n_sinktypes
    fprintf(fid, '%-10s  %10.4f  %10.4f  %10.4f  %15s\n', ...
        sink_names{k}, perm_results(k).obs_stat, ...
        perm_results(k).p_raw, perm_results(k).p_fdr, ...
        yesno(perm_results(k).p_fdr < ALPHA));
end
fprintf(fid, '\n');

fprintf(fid, '--- MULTIPLE COMPARISONS CORRECTION ---\n');
fprintf(fid, 'Method : Benjamini-Hochberg FDR\n');
fprintf(fid, 'Alpha  : %.2f\n', ALPHA);
fprintf(fid, 'N tests: %d (one per sink type; applied separately to LME\n', n_sinktypes);
fprintf(fid, '         and permutation p-values)\n\n');

fprintf(fid, '--- METHODOLOGICAL NOTES ---\n');
fprintf(fid, '1. Sessions are nested within subjects. fitlme with a random\n');
fprintf(fid, '   intercept per subject accounts for between-animal variance,\n');
fprintf(fid, '   so the fixed group effect is tested appropriately.\n');
fprintf(fid, '2. Proportions were logit-transformed prior to LME fitting.\n');
fprintf(fid, '   Epsilon = %.4f added to avoid logit(0)/logit(1).\n', EPSILON);
fprintf(fid, '3. The permutation test shuffles group labels at the subject\n');
fprintf(fid, '   level (not session level), preserving the nested structure.\n');
fprintf(fid, '4. FDR correction applied across %d sink types separately for\n', n_sinktypes);
fprintf(fid, '   LME and permutation p-values.\n\n');

fprintf(fid, '--- RAW PER-SESSION DATA ---\n');
fprintf(fid, '%-12s  %-10s  %-14s  %10s  %10s  %10s  %10s\n', ...
    'Session','Subject','Group','N_ripples','prop_rad','prop_slm','prop_other');
fprintf(fid, '%s\n', repmat('-',1,82));
for s = 1:valid_rows
    fprintf(fid, '%-12s  %-10s  %-14s  %10d  %10.4f  %10.4f  %10.4f\n', ...
        char(prop_table.session(s)), char(prop_table.subject(s)), ...
        char(prop_table.group(s)),   prop_table.n_ripples(s), ...
        prop_table.prop_rad(s), prop_table.prop_slm(s), prop_table.prop_other(s));
end

fclose(fid);
fprintf('\nReport saved to: %s\n', report_path);


%% ── 7. Figure: subject-mean dots + session scatter ───────────────────────
%
%  Large filled dot  = per-subject mean across sessions  (the LME unit)
%  Small faded dot   = individual session value          (for transparency)
%  Error bar         = SEM across subjects within group
%
figA = figure('Color','w','Position',[50 50 350*n_sinktypes 520]);

for k = 1:n_sinktypes
    ax = subplot(1, n_sinktypes, k);
    hold(ax, 'on');
    all_subj_means = [];

    for g = 1:n_groups
        grp_mask    = prop_table.group == group_order(g);
        subj_in_grp = unique(prop_table.subject(grp_mask));
        n_s         = numel(subj_in_grp);
        subj_means  = zeros(n_s, 1);

        for s = 1:n_s
            sess_mask     = prop_table.subject == subj_in_grp(s) & grp_mask;
            sess_vals     = prop_table.(sink_vars{k})(sess_mask);
            subj_means(s) = mean(sess_vals);

            % Small faded dots: individual sessions
            x_sess = g + (rand(size(sess_vals)) - 0.5) * 0.18;
            scatter(ax, x_sess, sess_vals, 25, groupcolors(g,:), 'filled', ...
                'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none');

            % Thin line connecting sessions of the same subject
            if numel(sess_vals) > 1
                plot(ax, x_sess, sess_vals, '-', ...
                    'Color', [groupcolors(g,:), 0.2], 'LineWidth', 0.6);
            end
        end

        all_subj_means = [all_subj_means; subj_means]; %#ok<AGROW>

        % Large filled dot: per-subject mean, jittered
        x_subj = g + (rand(n_s,1) - 0.5) * 0.18;
        scatter(ax, x_subj, subj_means, 80, groupcolors(g,:), 'filled', ...
            'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor','w', 'LineWidth', 0.8);

        % Group mean ± SEM across subjects (the quantity the LME tests)
        gm  = mean(subj_means);
        sem = std(subj_means) / sqrt(n_s);
        plot(ax, [g-0.22, g+0.22], [gm, gm], '-', ...
            'Color', groupcolors(g,:), 'LineWidth', 3);
        plot(ax, [g, g], [gm-sem, gm+sem], '-', ...
            'Color', groupcolors(g,:), 'LineWidth', 2);
    end

    % Significance annotation (based on LME p-value)
    if ~isempty(all_subj_means)
        y_range = range(all_subj_means);
        if y_range == 0, y_range = 0.1; end
        y_top = max(all_subj_means) + y_range * 0.15;

        if n_groups == 2
            p_annot = lme_results(k).p_fdr;
            star    = p_to_star(p_annot);
            if isempty(star), star = 'ns'; end
            tick_h  = y_range * 0.04;
            plot(ax, [1,2], [y_top,y_top],          'k-','LineWidth',1.2);
            plot(ax, [1,1], [y_top-tick_h,y_top],   'k-','LineWidth',1.2);
            plot(ax, [2,2], [y_top-tick_h,y_top],   'k-','LineWidth',1.2);
            text(ax, 1.5, y_top+tick_h*0.5, star, ...
                'HorizontalAlignment','center','VerticalAlignment','bottom',...
                'FontSize',13,'FontWeight','bold');
            ylim(ax, [0, max(0.05, y_top + y_range*0.2)]);

        else
            bracket_level = 0;
            if ~isempty(lme_results(k).posthoc)
                ph      = lme_results(k).posthoc;
                n_pairs = size(ph,1);
                [~, sort_ord] = sort(abs(ph(:,1)-ph(:,2)));
                for r = 1:n_pairs
                    ri     = sort_ord(r);
                    p_bonf = min(ph(ri,4)*n_pairs, 1);
                    if p_bonf >= ALPHA, continue; end
                    bracket_level = bracket_level + 1;
                    x1 = ph(ri,1); x2 = ph(ri,2);
                    yb = y_top + y_range*0.12*bracket_level;
                    tick_h = y_range*0.03;
                    star   = p_to_star(p_bonf);
                    plot(ax,[x1,x2],[yb,yb],        'k-','LineWidth',1.2);
                    plot(ax,[x1,x1],[yb-tick_h,yb], 'k-','LineWidth',1.2);
                    plot(ax,[x2,x2],[yb-tick_h,yb], 'k-','LineWidth',1.2);
                    text(ax, mean([x1,x2]), yb+tick_h*0.5, star, ...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom',...
                        'FontSize',12,'FontWeight','bold');
                end
            end
            ylim(ax,[0, max(0.05, y_top+y_range*0.12*(bracket_level+2))]);
        end
    end

    set(ax,'XTick',1:n_groups,'XTickLabel',cellstr(group_order),...
           'Box','off','FontSize',11);
    ylabel(ax,'Proportion','FontSize',11);
    title(ax, sink_names{k},'FontSize',13,'FontWeight','bold');
    xlim(ax,[0.5, n_groups+0.5]);
    xlabel(ax, sprintf('LME p_{FDR}=%.4f  |  perm p_{FDR}=%.4f', ...
        lme_results(k).p_fdr, perm_results(k).p_fdr), ...
        'FontSize',9,'Color',[0.4 0.4 0.4]);
    hold(ax,'off');
end

sgtitle('Sink type proportions  —  large dots = subject means  |  small dots = sessions', ...
    'FontSize',12,'FontWeight','bold');
saveas(figA, fullfile(savepath,'FigA_Proportion_DotPlot.svg'));
saveas(figA, fullfile(savepath,'FigA_Proportion_DotPlot.jpg'));

%% ── Local functions ──────────────────────────────────────────────────────
function p_adj = deal_bh(p_raw)
    p_raw  = p_raw(:);
    m      = numel(p_raw);
    [p_sorted, sort_idx] = sort(p_raw);
    p_adj_sorted    = zeros(m,1);
    p_adj_sorted(m) = p_sorted(m);
    for i = (m-1):-1:1
        p_adj_sorted(i) = min(p_sorted(i)*m/i, p_adj_sorted(i+1));
    end
    p_adj_sorted       = min(p_adj_sorted,1);
    p_adj(sort_idx,1)  = p_adj_sorted;
    p_adj              = reshape(p_adj, size(p_raw));
end

function out = yesno(tf)
    if tf, out = 'YES *'; else, out = 'no'; end
end