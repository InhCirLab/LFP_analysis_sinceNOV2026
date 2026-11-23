

%% CELL-METRICS PROCESSING FOR DATA ANALYSIS:

%% 
mainFolder = uigetdir(); % Change this to your main folder path
folderInfo = dir(mainFolder);

% Filter only directories and exclude '.' and '..'
subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

% Get the names of the subfolders
lfppaths = {subfolders.name};

%% Extract and organize measures

varlist = ["ab_ratio" "burstIndex_Mizuseki2012" "burstIndex_Royer2012" "burstIndex_WAKEtheta" "cv2" "cv2_WAKEtheta",...
    "firingRate" "firingRateCV" "firingRate_WAKEtheta" "population_modIndex" "ripples_modulationIndex" "ripples_modulationRatio"...
    "thetaModulationIndex" "troughToPeak"];

%% SS
% varnames = ["ID", "expcond", "celltype", "spikech","layer","abratio", "burstidx_Mizu", "burstidx_Royer", "burstidx_wakeTheta", "cv2", "cv2_wakeTheta", "firingRate", ...
%     "firingRateCV", "firingRatetheta" ,"population_modidx", "ripples_modix", "ripples_modratio", "theta_modidx", "troughtopeak"];

layersnames = ["oriens", "pyr", "rad", "slm"];
SScellmeasures = [];
for ii = 1:length(lfppaths)
    try
    thisfilename = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".cell_metrics.cellinfo.mat");
    thisspikename = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".spikes.cellinfo.mat");
    thislayersname = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".hippocampalLayers.channelinfo.mat");
    if ~exist(thisfilename)
        continue
    else
        load(thisfilename)
        load(thisspikename)
        load(thislayersname)
    end
    thispeakvoltage = spikes.peakVoltage_sorted;
    thisfileclass = cell_metrics.ground_truth_classification.cell_subtypes;
    fieldidx = isfield(cell_metrics, varlist);
    thismeasures = zeros(size(thisfileclass, 1), size(varlist, 2));
    for jj = 1:size(varlist, 2)
        if fieldidx(jj)
        thismeasures(:, jj) = cell_metrics.(varlist(jj));
        else
        thismeasures(:, jj) = nan(size(thismeasures, 1), 1);
        end
    end
    spikech = spikes.phy_maxWaveformCh1.';
    layers = hippocampalLayers.layers{1,1}.all;
    layerchannels = strings([16 1]);
    spikelayers = zeros(size(spikech, 1), 1);
    for kk = 1:size(spikech, 1)
        [~, spikelayers(kk, 1)] = min(abs(layers - spikech(kk))); 
    end
    oriensidx = spikelayers == 1;pyridx = spikelayers == 2;
    radidx = spikelayers == 3;slmidx = spikelayers == 4;
    spikelayersnames = strings([size(spikech, 1), 1]);
    spikelayersnames(oriensidx) = layersnames(1);
    spikelayersnames(pyridx) = layersnames(2);
    spikelayersnames(radidx) = layersnames(3);
    spikelayersnames(slmidx) = layersnames(4);

    % ID
    thisID = repmat(lfppaths(ii), size(thismeasures, 1), 1);
    thisexpcond = repmat("SS", size(thismeasures,1), 1);
    thismeasurestable = array2table(thismeasures, 'VariableNames', varlist);
    thisIDtable = table(thisID, thisexpcond, thisfileclass,spikech, spikelayersnames, 'VariableNames',["ID", "expcond", "celltype","spikech","layer",]);
    thisrectable = [thisIDtable thismeasurestable];
    SScellmeasures = [SScellmeasures; thisrectable];
    end

end

%%  OS
varlist = ["ab_ratio" "burstIndex_Mizuseki2012" "burstIndex_Royer2012" "burstIndex_WAKEtheta" "cv2" "cv2_WAKEtheta",...
    "firingRate" "firingRateCV" "firingRate_WAKEtheta" "population_modIndex" "ripples_modulationIndex" "ripples_modulationRatio"...
    "thetaModulationIndex" "troughToPeak"];
layersnames = ["oriens", "pyr", "rad", "slm"];

OScellmeasures = [];
for ii = 1:length(lfppaths)
    thisfilename = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".cell_metrics.cellinfo.mat");
    thisspikename = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".spikes.cellinfo.mat");
    thislayersname = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".hippocampalLayers.channelinfo.mat");
    if ~exist(thisfilename)
        continue
    else
        load(thisfilename)
        load(thisspikename)
        load(thislayersname)
    end
    thispeakvoltage = spikes.peakVoltage_sorted;
    thisfileclass = cell_metrics.ground_truth_classification.cell_subtypes;
    fieldidx = isfield(cell_metrics, varlist);
    thismeasures = zeros(size(thisfileclass, 1), size(varlist, 2));
    for jj = 1:size(varlist, 2)
        if fieldidx(jj)
        thismeasures(:, jj) = cell_metrics.(varlist(jj));
        else
        thismeasures(:, jj) = nan(size(thismeasures, 1), 1);
        end
    end
    spikech = spikes.phy_maxWaveformCh1.';
    layers = hippocampalLayers.layers{1,1}.all;
    layerchannels = strings([16 1]);
    spikelayers = zeros(size(spikech, 1), 1);
    for kk = 1:size(spikech, 1)
        [~, spikelayers(kk, 1)] = min(abs(layers - spikech(kk))); 
    end
    oriensidx = spikelayers == 1;pyridx = spikelayers == 2;
    radidx = spikelayers == 3;slmidx = spikelayers == 4;
    spikelayersnames = strings([size(spikech, 1), 1]);
    spikelayersnames(oriensidx) = layersnames(1);
    spikelayersnames(pyridx) = layersnames(2);
    spikelayersnames(radidx) = layersnames(3);
    spikelayersnames(slmidx) = layersnames(4);

    % ID
    thisID = repmat(lfppaths(ii), size(thismeasures, 1), 1);
    thisexpcond = repmat("OS", size(thismeasures,1), 1);
    thismeasurestable = array2table(thismeasures, 'VariableNames', varlist);
    thisIDtable = table(thisID, thisexpcond, thisfileclass,spikech, spikelayersnames, 'VariableNames',["ID", "expcond", "celltype","spikech","layer",]);
    thisrectable = [thisIDtable thismeasurestable];
    OScellmeasures = [OScellmeasures; thisrectable];

end

%%  OD
varlist = ["ab_ratio" "burstIndex_Mizuseki2012" "burstIndex_Royer2012" "burstIndex_WAKEtheta" "cv2" "cv2_WAKEtheta",...
    "firingRate" "firingRateCV" "firingRate_WAKEtheta" "population_modIndex" "ripples_modulationIndex" "ripples_modulationRatio"...
    "thetaModulationIndex" "troughToPeak"];
layersnames = ["oriens", "pyr", "rad", "slm"];

ODcellmeasures = [];
for ii = 1:length(lfppaths)
    thisfilename = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".cell_metrics.cellinfo.mat");
    thisspikename = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".spikes.cellinfo.mat");
    thislayersname = strcat(mainFolder, "\", lfppaths(ii), "\", lfppaths(ii), ".hippocampalLayers.channelinfo.mat");
    if ~exist(thisfilename)
        continue
    else
        load(thisfilename)
        load(thisspikename)
        load(thislayersname)
    end
    thispeakvoltage = spikes.peakVoltage_sorted;
    thisfileclass = cell_metrics.ground_truth_classification.cell_subtypes;
    fieldidx = isfield(cell_metrics, varlist);
    thismeasures = zeros(size(thisfileclass, 1), size(varlist, 2));
    for jj = 1:size(varlist, 2)
        if fieldidx(jj)
        thismeasures(:, jj) = cell_metrics.(varlist(jj));
        else
        thismeasures(:, jj) = nan(size(thismeasures, 1), 1);
        end
    end
    spikech = spikes.phy_maxWaveformCh1.';
    layers = hippocampalLayers.layers{1,1}.all;
    layerchannels = strings([16 1]);
    spikelayers = zeros(size(spikech, 1), 1);
    for kk = 1:size(spikech, 1)
        [~, spikelayers(kk, 1)] = min(abs(layers - spikech(kk))); 
    end
    oriensidx = spikelayers == 1;pyridx = spikelayers == 2;
    radidx = spikelayers == 3;slmidx = spikelayers == 4;
    spikelayersnames = strings([size(spikech, 1), 1]);
    spikelayersnames(oriensidx) = layersnames(1);
    spikelayersnames(pyridx) = layersnames(2);
    spikelayersnames(radidx) = layersnames(3);
    spikelayersnames(slmidx) = layersnames(4);

    % ID
    thisID = repmat(lfppaths(ii), size(thismeasures, 1), 1);
    thisexpcond = repmat("OD", size(thismeasures,1), 1);
    thismeasurestable = array2table(thismeasures, 'VariableNames', varlist);
    thisIDtable = table(thisID, thisexpcond, thisfileclass,spikech, spikelayersnames, 'VariableNames',["ID", "expcond", "celltype","spikech","layer",]);
    thisrectable = [thisIDtable thismeasurestable];
    ODcellmeasures = [ODcellmeasures; thisrectable];

end

%% Join all data:

allrec = [SScellmeasures; OScellmeasures; ODcellmeasures];
datapath = uigetdir();
allfile = strcat(datapath, "\all_cellmetrics_withlayer.xlsx");
writetable(allrec, allfile)

measures = table2array(allrec(:, 4));table2array(allrec(:, 6:end));
nanidx = isnan(measures); nanidx = sum(nanidx, 2) > 0;
measureskeep = measures(~nanidx, :);
% UMAP:

id = string(allrec.expcond);
id = id(~nanidx);

% [umap_embedding, umap_params] = run_umap(measureskeep, 'n_components', 2);
% % 3. Extract group labels and define colors
% groupLabels = id;  % Assuming 'Group' is categorical or a cell array of strings
% groups = unique(groupLabels);
% colors = lines(numel(groups));  % Generate distinct colors
% % 4. Plot the UMAP results
% figure; hold on;
% for i = 1:numel(groups)
%     idx = groupLabels == groups(i);
%     scatter(umap_embedding(idx,1), umap_embedding(idx,2), 10, colors(i,:), 'filled');
% end
% xlabel('UMAP1');
% ylabel('UMAP2');
% title('UMAP Embedding of Sharp-Wave Ripple Measures');
% legend(groups, 'Location', 'best');
% % grid on;
% hold off;

%% Save data:

savefolder = uigetdir();
celldatafile = strcat(savefolder, "\allcellsdata.xlsx");
writetable(allrec, celldatafile)


%% EXCITATORY NEURONS:
% 
% % Separar en fx de excit e inh:
% celltype = unique(allrec.celltype);
% camkidx = string(allrec.celltype) == celltype(1);
% excneurons = allrec(camkidx, :);
% excmeasures = table2array(excneurons(:, 4:end));
% nanidx = isnan(excmeasures); nanidx = sum(nanidx, 2) > 0;
% excmeasureskeep = excmeasures(~nanidx, :);
% 
% id = string(excneurons.expcond);
% id = id(~nanidx);

% [umap_embedding, umap_params] = run_umap(excmeasureskeep, 'n_components', 2);
% % 3. Extract group labels and define colors
% groupLabels = id;  % Assuming 'Group' is categorical or a cell array of strings
% groups = unique(groupLabels);
% colors = lines(numel(groups));  % Generate distinct colors
% % 4. Plot the UMAP results
% figure; hold on;
% for i = 1:numel(groups)
%     idx = groupLabels == groups(i);
%     scatter(umap_embedding(idx,1), umap_embedding(idx,2), 10, colors(i,:), 'filled');
% end
% xlabel('UMAP1');
% ylabel('UMAP2');
% title('UMAP Embedding of Sharp-Wave Ripple Measures');
% legend(groups, 'Location', 'best');
% % grid on;
% hold off;

%% INHIBITORY NEURONS

% % Separar en fx de excit e inh:
% celltype = unique(allrec.celltype);
% inhidx = (string(allrec.celltype) ~= celltype(1)) & (string(allrec.celltype) ~= celltype(3));
% inhneurons = allrec(inhidx, :);
% inhmeasures = table2array(inhneurons(:, 4:end));
% nanidx = isnan(inhmeasures); nanidx = sum(nanidx, 2) > 0;
% inhmeasureskeep = inhmeasures(~nanidx, :);
% 
% id = string(inhneurons.expcond);
% id = id(~nanidx);
% 
% [umap_embedding, umap_params] = run_umap(inhmeasureskeep, 'n_components', 2);
% % 3. Extract group labels and define colors
% groupLabels = id;  % Assuming 'Group' is categorical or a cell array of strings
% groups = unique(groupLabels);
% colors = lines(numel(groups));  % Generate distinct colors
% % 4. Plot the UMAP results
% figure; hold on;
% for i = 1:numel(groups)
%     idx = groupLabels == groups(i);
%     scatter(umap_embedding(idx,1), umap_embedding(idx,2), 10, colors(i,:), 'filled');
% end
% xlabel('UMAP1');
% ylabel('UMAP2');
% title('UMAP Embedding of Sharp-Wave Ripple Measures');
% legend(groups, 'Location', 'best');
% % grid on;
hold off;

%% ANALIZAR Y PLOTEAR MEDIDAS

% % --- Step 3: Run Plotting for Both Types ---
% 
% fprintf('\nRunning analysis for EXCITATORY cells...\n');
% plot_cell_measuresII(excneurons, 'Excitatory');
% 
% fprintf('\nRunning analysis for INHIBITORY cells...\n');
% plot_cell_measuresII(inhneurons, 'Inhibitory');
% 
% pvcellsidx = inhneurons.celltype == "PV+";
% sstcellsidx = inhneurons.celltype == "SST+";
% id2cellsidx = inhneurons.celltype == "ID2+";
% 
% pvcells = inhneurons(pvcellsidx, :);
% sstcells = inhneurons(sstcellsidx, :);
% id2cells = inhneurons(id2cellsidx, :);
% 
% plot_cell_measuresII(pvcells, 'PV+')
% plot_cell_measuresII(sstcells, 'SST+')
% plot_cell_measuresII(id2cells, 'ID2+')

%% IDENTIFICAR CANAL 
celltype = unique(allrec.celltype);
inhidx = (string(allrec.celltype) ~= celltype(1)) & (string(allrec.celltype) ~= celltype(3));
inhneurons = allrec(inhidx, :);

layer = unique(inhneurons.layer);
oridx = string(inhneurons.layer) == layer(1);
pyridx = string(inhneurons.layer) == layer(2);
radidx = string(inhneurons.layer) == layer(3);
slmidx = string(inhneurons.layer) == layer(4);

orneurons = inhneurons(oridx, :);
pyrneurons = inhneurons(pyridx, :);
radneurons = inhneurons(radidx, :);
slmneurons = inhneurons(slmidx, :);

plot_cell_measuresII(orneurons, 'Oriens')
plot_cell_measuresII(pyrneurons, 'Pyr')
plot_cell_measuresII(radneurons, 'Rad')
plot_cell_measuresII(slmneurons, 'Slm')

% Descartando excitatorias que no estén en piramidal..

celltype = unique(allrec.celltype);
camkidx = (string(allrec.celltype) == celltype(1)) & (string(allrec.layer) ==  "pyr");
excneurons = allrec(camkidx, :);

%% Excitatory neurons: deep vs sup


SS = ["C089", "C090", "C091", "C199", "E075"];
OS = ["C093", "C094", "C095", "D094"];
OD = ["D097", "D096", "E078", "E079"];
% all_ripples = SexHormonesMalesripplesoriginal;

subject = repmat("", size(allrec, 1), 1);
allrec.ID = string(allrec.ID);
for ii = 1:length(SS)
    thissubject = contains(allrec.ID, SS(ii));
    subject(thissubject) = SS(ii);
end

for ii = 1:length(OS)
    thissubject = contains(allrec.ID, OS(ii));
    subject(thissubject) = OS(ii);
end

for ii = 1:length(OD)
    thissubject = contains(allrec.ID, OD(ii));
    subject(thissubject) = OD(ii);
end

allrec.subject = subject;
%%
camkneurons = allrec(contains(allrec.celltype, "CAMK2"), :);
deepneurons = camkneurons(camkneurons.celltype == "CAMK2_DEEP", :);
supneurons = camkneurons(camkneurons.celltype == "CAMK2_SUP", :);

varnames = allrec.Properties.VariableNames(:, 6:end-1);
for ii = 1:length(varnames)
    % deepresults{ii} = analyze_ripple_measure_lme(deepneurons, varnames{ii});
    supresults{ii} = analyze_ripple_measure_lme_2way(supneurons, varnames{ii});
    % camkresults{ii} = analyze_ripple_measure_lme(camkneurons, varnames{ii});
    % pause(5)
    % close all
end

savePath = pwd;

figs = findobj('Type','figure');
for i = 1:numel(figs)
    f = figs(i);
    name = sprintf('Figure_%d', f.Number);
    saveas(f, fullfile(savePath, [name '.svg']));
end
close all

%% REMOVE OUTLIERS FOR CELL METRICS

%% 5. Hierarchical robust outlier removal (per-recording -> per-subject -> per-group -> global)
[dataTable, mouseCol, recCol] = deal(camkneurons, 'subject', 'ID');

% user parameters (adjust if needed)
vars = varnames;
logVars = varnames;   % variables to log-transform before MAD
thresh = 3.5;            % robust z threshold
minN_rec = 20;           % min samples to compute per-recording MAD
minN_sub = 20;           % min to compute per-subject
minN_group = 20;         % min to compute per-group (expcond)

n = height(dataTable);
flagsMat = false(n, numel(vars));
out.union = struct();

% helper inline function to compute robust z and flags
robustFlags = @(x,mask) deal( ...
    arrayfun(@(v) (abs((v - median(v))./(max(1.4826*median(abs(v - median(v))), eps)))) > thresh, v2cell(x(mask))), ...
    []); % not used directly, kept for clarity

% loop variables
uniqueRecs = unique(dataTable.(recCol));
uniqueSubs = unique(dataTable.(mouseCol));
uniqueGroups = unique(dataTable.expcond);

for vi = 1:numel(vars)
    var = vars{vi};
    xAll = dataTable.(var);
    flagsVar = false(n,1);

    isLog = ismember(var, logVars);

    % iterate recordings
    for r = 1:numel(uniqueRecs)
        recID = uniqueRecs(r);
        idxRec = ismember(dataTable.(recCol), recID);
        % valid numeric positive values
        valid = idxRec & ~isnan(xAll);
        if isLog
            valid = valid & (xAll > 0);
        end
        x = xAll(valid);
        if isempty(x), continue; end

        % choose level depending on sample size
        if numel(x) >= minN_rec
            % per-recording
            if isLog, xproc = log10(x); else xproc = x; end
        else
            % fallback to subject
            subjID = unique(dataTable.(mouseCol)(idxRec));
            idxSub = ismember(dataTable.(mouseCol), subjID);
            validSub = idxSub & ~isnan(xAll);
            if isLog, validSub = validSub & (xAll > 0); end
            xSub = xAll(validSub);
            if numel(xSub) >= minN_sub
                if isLog, xproc = log10(xSub); else xproc = xSub; end
                valid = validSub;  % we'll mark rows within this subject level
            else
                % fallback to group (expcond)
                grpID = unique(dataTable.expcond(idxRec));
                idxGrp = ismember(dataTable.expcond, grpID);
                validGrp = idxGrp & ~isnan(xAll);
                if isLog, validGrp = validGrp & (xAll > 0); end
                xGrp = xAll(validGrp);
                if numel(xGrp) >= minN_group
                    if isLog, xproc = log10(xGrp); else xproc = xGrp; end
                    valid = validGrp;
                else
                    % fallback to global
                    validGlob = ~isnan(xAll);
                    if isLog, validGlob = validGlob & (xAll > 0); end
                    xGlob = xAll(validGlob);
                    if isempty(xGlob), continue; end
                    if isLog, xproc = log10(xGlob); else xproc = xGlob; end
                    valid = validGlob;
                end
            end
        end

        % compute robust z using median & MAD
        med = median(xproc);
        madv = median(abs(xproc - med));
        scaledMAD = 1.4826 * madv;
        if scaledMAD == 0
            % fallback to std (rare) or very small epsilon
            s = std(xproc);
            if s > 0, scaledMAD = s; else scaledMAD = eps; end
        end

        % compute z only for the selected 'valid' indices and mark outliers
        idxValidList = find(valid);          % absolute indices into table
        xvalid_proc = (isLog) * 0; %#ok<NASGU> % placeholder so editor doesn't complain
        if isLog
            xvalid_proc = log10(xAll(valid));
        else
            xvalid_proc = xAll(valid);
        end
        z = (xvalid_proc - med) / scaledMAD;
        outIdx = idxValidList(abs(z) > thresh);
        flagsVar(outIdx) = true;
    end

    flagsMat(:,vi) = flagsVar;
    out.union.(var) = find(flagsVar);
    fprintf('%s: %d flagged rows\n', var, numel(out.union.(var)));
    if ~isempty(out.union.(var)), fprintf('  first indices: %s\n', mat2str(out.union.(var)(1:min(10,end)))); end
end

% union across variables
unionIdx = any(flagsMat,2);
discardIdx = find(unionIdx);
fprintf('TOTAL unique rows to discard (union across vars): %d\n', numel(discardIdx));

% remove flagged rows from all_ripples
rowidx = false(height(camkneurons),1);
rowidx(discardIdx) = true;
camkneurons = camkneurons(~rowidx, :);
% filteredRipples = all_ripples;  % keep same name semantics as your code

%% 

deepneurons = camkneurons(contains(camkneurons.celltype, "DEEP"), :);
supneurons = camkneurons(contains(camkneurons.celltype, "SUP"), :);

varnames = deepneurons.Properties.VariableNames(:, 6:19);
for ii = 1:length(varnames)
    analyze_ripple_measure_lme_figures_generalized(deepneurons, varnames{ii});
    % pause(5)
end

for ii = 1:length(varnames)
    analyze_ripple_measure_lme_figures_generalized(supneurons, varnames{ii});
    % pause(5)
end

% analyze_ripple_measure_lme_figures(deepneurons, 'firingRate_WAKEtheta');

%% 

pvneurons = allrec(contains(allrec.celltype, "PV"), :);
sstneurons = allrec(contains(allrec.celltype, "SST"), :);

for ii = 1:length(varnames)
    analyze_ripple_measure_lme_figures_generalized(pvneurons, varnames{ii});
    % pause(5)
end

for ii = 1:length(varnames)
    analyze_ripple_measure_lme_figures_generalized(sstneurons, varnames{ii});
    % pause(5)
end

% analyze_ripple_measure_lme_figures(sstneurons, 'ab_ratio')

%% 2 way anova

for ii = 1:length(varnames)
    analyze_ripple_measure_lme_2way_RutORXproject(camkneurons, varnames{ii}, 'celltype', unique(camkneurons.celltype));
    pause(10)
end

% 1y anova for all camkneurons

for ii = 1:length(varnames)
    analyze_ripple_measure_lme_figures(camkneurons, varnames{ii});
    pause(5)
end

% ID2 neurons

id2neurons = allrec(allrec);%que coño es esto jajjajja
for ii = 1:length(varnames)
        analyze_ripple_measure_lme_figures(camkneurons, varnames{ii});
    pause(5)
end