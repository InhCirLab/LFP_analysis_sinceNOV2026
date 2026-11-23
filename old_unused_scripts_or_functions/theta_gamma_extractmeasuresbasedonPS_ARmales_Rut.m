

%% Compute Power Spectrum from 1 to 120 Hz & Compute Theta and Gamma Measures for different events %% 

%% 1 . Load and extract each groups data:

PVdata = ps1to120("PV");
SSTdata = ps1to120("SST");
% ODdata = ps1to120("OD");

orallrec = SSdata.or.allrec;
corrdata = corrcoef(SSdata.or.allrec);

freqs = SSdata.freqs;
% plot_correlation_heatmap(freqs, freqs, corrdata, "Fq (Hz)", "Fq (Hz)")
% imagesc(freqs, freqs, corrdata)

%% 2. Data that I want to extract:

% 2.1. Layers:

datafields = fields(SSdata);
layernames = datafields(2:5);

% 2.2. Measures:
measuresfields = fields(SSdata.or);


%% 3. Compute Peak Frequency of desired band and extract max power

% Modify categorical table (add another column for control vs AR):
% 
% PVctrlidx = contains(PVdata.idtable.ID, "Control");
% SSTctrlidx = contains(SSTdata.idtable.ID, "Control");
% 
% PVgroup = repmat("", length(PVctrlidx),1); PVgroup(PVctrlidx) = "Control";
% PVgroup(~PVctrlidx) = "AR"; PVdata.idtable.group = PVgroup;
% 
% SSTgroup = repmat("", length(SSTctrlidx),1); SSTgroup(SSTctrlidx) = "Control";
% SSTgroup(~SSTctrlidx) = "AR"; SSTdata.idtable.group = SSTgroup;

% 2.1. Layers:

datafields = fields(SSTdata);
layernames = datafields(2:5);

% 2.2. Measures:
measuresfields = fields(SSTdata.or);

freqs = PVdata.freqs;

% deltaidx = freqs >= 1 & freqs <= 3.5;
% thetaidx = freqs >= 6 & freqs <= 12;
% betaidx = freqs >= 12 & freqs <= 30;
% lowgammaidx = freqs >= 30 & freqs <= 60;
% mediumgammaidx = freqs >= 60 & freqs <= 90;
highgammaidx = freqs >= 90 & freqs <= 120;

freqidx = highgammaidx;
freqvector = freqs(freqidx);
thisfreq = freqvector;

for ii = 1:length(measuresfields)

    PVnanRows = any(isnan(PVdata.(layernames{1}).(measuresfields{ii})), 2);
    PVlayerpow = zeros(sum(~PVnanRows), length(layernames));
    PVlayeridx = PVlayerpow; PVlayerpeakfq = PVlayerpow;
    PVgrouplabelmod = PVgroup(~PVnanRows);

    SSTnanRows = any(isnan(SSTdata.(layernames{1}).(measuresfields{ii})), 2);
    SSTlayerpow = zeros(sum(~SSTnanRows), length(layernames));
    SSTlayeridx = SSTlayerpow; SSTlayerpeakfq = SSTlayerpow;
    SSTgrouplabelmod = SSTgroup(~SSTnanRows);

    for jj = 1:length(layernames)
        if ~contains(measuresfields(ii), "db")
            % Check for nans and keep data without nans: 
            PVnanRows = any(isnan(PVdata.(layernames{jj}).(measuresfields{ii})), 2);
            tempPVdata = pow2db(PVdata.(layernames{jj}).(measuresfields{ii})(~PVnanRows, freqidx));
            SSTnanRows = any(isnan(SSTdata.(layernames{jj}).(measuresfields{ii})), 2);
            tempSSTdata = pow2db(SSTdata.(layernames{jj}).(measuresfields{ii})(~SSTnanRows, freqidx));
        else
            % Check for nans and keep data without nans: 
            PVnanRows = any(isnan(PVdata.(layernames{jj}).(measuresfields{ii})), 2);
            tempPVdata = PVdata.(layernames{jj}).(measuresfields{ii})(~PVnanRows, freqidx);
            SSTnanRows = any(isnan(SSTdata.(layernames{jj}).(measuresfields{ii})), 2);
            tempSSTdata = SSTdata.(layernames{jj}).(measuresfields{ii})(~SSTnanRows, freqidx);
        end

        % Compute power and index to find peak frequency:
        [PVlayerpow(:, jj), PVlayeridx(:, jj)] = max(tempPVdata,[], 2);
        [SSTlayerpow(:, jj), SSTlayeridx(:, jj)] = max(tempSSTdata,[], 2);

        PVlayerpeakfq(:, jj) = thisfreq(PVlayeridx(:, jj).').';
        SSTlayerpeakfq(:, jj) = thisfreq(SSTlayeridx(:, jj).').';        
 
    end   
   
    celltypelabel = table([repmat("PV", size(PVlayerpeakfq, 1), 1);...
        repmat("SST", size(SSTlayerpeakfq, 1), 1)], 'VariableNames', "CellType");
    grouplabel = table([PVgrouplabelmod;SSTgrouplabelmod], 'VariableNames', "Group");

    peakfqdata = [PVlayerpeakfq; SSTlayerpeakfq];
    powtopeakdata = [PVlayerpow; SSTlayerpow];

    peakfq = array2table(peakfqdata, 'VariableNames', layernames);
    powtopeak = array2table(powtopeakdata, 'VariableNames', layernames);

    thismeasure_peakfq = [celltypelabel grouplabel peakfq];
    thismeasure_powtopeak = [celltypelabel grouplabel powtopeak];
    % thismeasure_peakfq.Group = categorical(thismeasure_peakfq.Group);
    % thismeasure_powtopeak.Group = categorical(thismeasure_powtopeak.Group);

    freqmeasures.(measuresfields{ii}).peakfq = thismeasure_peakfq;
    freqmeasures.(measuresfields{ii}).powtopeak = thismeasure_powtopeak;

end

% deltameasures = freqmeasures;
% thetameasures = freqmeasures;
% betameasures = freqmeasures;
% lowgammameasures = freqmeasures;
% mediumgammameasures = freqmeasures;
highgammameasures = freqmeasures;


%% 5. Plot and analyze data:


% 2.1. Layers:
% 
% datafields = fields(SSdata);
% layernames = datafields(2:5);
% 
% % 2.2. Measures:
% measuresfields = fields(SSdata.or);

% thismeasure_peakfq.Group = categorical(thismeasure_peakfq.Group);
% savepath = uigetdir();
% hippocampal_permutation_analysisv3(thismeasure_peakfq, "Peak Frequency Loc uv", "Hz", savepath)

for ii = 1:length(measuresfields)
    msg(strcat("Choose path folder for ", measuresfields(ii)))
    savepath = uigetdir();
    % Peak Frequency
    peakfqtitle = strcat("Peak Fq ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(deltameasures.(measuresfields{ii}).peakfq, peakfqtitle, "Hz", savepath)
    % Pow to Peak
    powtopeaktitle = strcat("Power to Peak ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(deltameasures.(measuresfields{ii}).powtopeak, powtopeaktitle, "a.u.", savepath)
end

for ii = 1:length(measuresfields)
    msg(strcat("Choose path folder for ", measuresfields(ii)))
    savepath = uigetdir();
    % Peak Frequency
    peakfqtitle = strcat("Peak Fq ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(thetameasures.(measuresfields{ii}).peakfq, peakfqtitle, "Hz", savepath)
    % Pow to Peak
    powtopeaktitle = strcat("Power to Peak ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(thetameasures.(measuresfields{ii}).powtopeak, powtopeaktitle, "a.u.", savepath)
end

for ii = 1:length(measuresfields)
    msg(strcat("Choose path folder for ", measuresfields(ii)))
    savepath = uigetdir();
    % Peak Frequency
    peakfqtitle = strcat("Peak Fq ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(betameasures.(measuresfields{ii}).peakfq, peakfqtitle, "Hz", savepath)
    % Pow to Peak
    powtopeaktitle = strcat("Power to Peak ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(betameasures.(measuresfields{ii}).powtopeak, powtopeaktitle, "a.u.", savepath)
end

for ii = 1:length(measuresfields)
    msg(strcat("Choose path folder for ", measuresfields(ii)))
    savepath = uigetdir();
    % Peak Frequency
    peakfqtitle = strcat("Peak Fq ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(lowgammameasures.(measuresfields{ii}).peakfq, peakfqtitle, "Hz", savepath)
    % Pow to Peak
    powtopeaktitle = strcat("Power to Peak ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(lowgammameasures.(measuresfields{ii}).powtopeak, powtopeaktitle, "a.u.", savepath)
end


for ii = 1:length(measuresfields)
    msg(strcat("Choose path folder for ", measuresfields(ii)))
    savepath = uigetdir();
    % Peak Frequency
    peakfqtitle = strcat("Peak Fq ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(mediumgammameasures.(measuresfields{ii}).peakfq, peakfqtitle, "Hz", savepath)
    % Pow to Peak
    powtopeaktitle = strcat("Power to Peak ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(mediumgammameasures.(measuresfields{ii}).powtopeak, powtopeaktitle, "a.u.", savepath)
end

for ii = 1:length(measuresfields)
    msg(strcat("Choose path folder for ", measuresfields(ii)))
    savepath = uigetdir();
    % Peak Frequency
    peakfqtitle = strcat("Peak Fq ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(highgammameasures.(measuresfields{ii}).peakfq, peakfqtitle, "Hz", savepath)
    % Pow to Peak
    powtopeaktitle = strcat("Power to Peak ", measuresfields(ii));
    ARgroups_saveandplotsthetagamma(highgammameasures.(measuresfields{ii}).powtopeak, powtopeaktitle, "a.u.", savepath)
end


