
%% %%%%%%%%% UMAP FINAL SCRIPT PRE-PROCESSING AND ANALYSIS %%%%%%%%%%%% %%

%% 1. Get all paths:

for ii = 1:3
    thisFolder = uigetdir(); % Change this to your main folder path
    mainFolder{ii} = thisFolder;
    folderInfo = dir(mainFolder{ii});
    
    % Filter only directories and exclude '.' and '..'
    subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'
    
    % Get the names of the subfolders
    lfppaths{ii} = {subfolders.name};
end
ripplefolder = uigetdir();

%% 2. Recompute ripples

SSlfppaths = lfppaths{1}; OSlfppaths = lfppaths{2}; ODlfppaths = lfppaths{3};
SSslmsize = nan(length(SSlfppaths), 1); SSradsize = nan(length(SSlfppaths), 1);
SSothersize = nan(length(SSlfppaths), 1);
rippleMasterDetector_threshold = [2 5];
passband = [80 240];
SS_ripplesraw = [];
SS_ripplesfilter = [];
SStableraw = [];
SStablefilt = [];
modelpath = 'C:\Users\PC\cnn-matlab-master\cnn';
cnnpath = 'C:\Users\PC\anaconda3\envs\cnn-env\python.exe'; 


for ii = 1:length(SSlfppaths)
    % Load ripple lfp params:
    thisfoldername = strcat(mainFolder{1}, "\", SSlfppaths(ii));
    thisfilename = strcat(thisfoldername, "\", SSlfppaths(ii), ".hippocampalLayers.channelinfo.mat");
    load(thisfilename);
    rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
    radch = hippocampalLayers.layers{1,1}.radiatum;
    SWChannel = hippocampalLayers.layers{1,1}.radiatum;
    slmch = hippocampalLayers.layers{1,1}.slm;
    cd(thisfoldername)
    % To get ripple data excluding lfp from locomotion
    thetafilename = strcat(thisfoldername, "\", SSlfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(thetafilename);
    tsthetaevents = thisthetaEpochs.thetaRun.newintervals;
    % For CSD data:
    lfp = getLFP('all');
    try
        SSripples{ii} = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,...
            'detector', 'customcnn','exec_env',cnnpath, 'model_file', modelpath,...
            'passband',passband,'force',true,'skipStimulationPeriods',true,...
            'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold', false,...
            'excludeIntervals', tsthetaevents);
        [SScsd{ii}, SSlfpAvg{ii}, SScsd_events{ii}, SSlfp_events{ii}] = buzcode_eacheventCSD(lfp, SSripples{ii}.peaks);
        close all
    
        tsthisripple = nan(length(SSripples{ii}.peaks), 2);
        durthisripple = tsthisripple(:, 1);
        rangethisripple = durthisripple;

        for jj = 1:size(SScsd_events{ii},3)
            SScsdsummary{ii}.tsthisripple(jj, :) = SSripples{ii}.timestamps(jj, :);
            SScsdsummary{ii}.durthisripple(jj) = diff(SScsdsummary{ii}.tsthisripple(jj, :))*1000; % paso a mseg para que round fx bien
            SScsdsummary{ii}.rangethisripple(jj) = round(SScsdsummary{ii}.durthisripple(jj)/2);
            SScsdsummary{ii}.csdtsthisrippleidx = SScsd{ii}.timestamps >= -SScsdsummary{ii}.rangethisripple(jj) & SScsd{ii}.timestamps <= SScsdsummary{ii}.rangethisripple(jj);
            SScsdsummary{ii}.csdsummary(jj, :) = mean(SScsd_events{ii}(SScsdsummary{ii}.csdtsthisrippleidx, :, jj), 1);
        end

        % 2.3. Categorize ripples according to their sink (rad or slm)
        [val, idx] = min(SScsdsummary{ii}.csdsummary, [], 2);
        slmpeaks = (idx >= radch + 1); fprintf("SLM ripples:  %4.2f ", sum(slmpeaks));
        SSslmsize(ii) = sum(slmpeaks);  
        radpeaks = idx >= rippleChannel + 1 & idx <= radch; fprintf("RAD ripples:  %4.2f ", sum(radpeaks));
        SSradsize(ii) = sum(radpeaks);
        otherpeaks = ~radpeaks & ~slmpeaks; fprintf("OTHER ripples:  %4.2f ", sum(otherpeaks));
        SSothersize(ii) = sum(otherpeaks);
        if sum(radpeaks) == 0
            warning("No RAD ripples for this file")
            SSradcsd{ii} = [];
        else
            SSradcsd{ii} = buzcode_eventCSD(lfp, SSripples{ii}.peaks(radpeaks));
        end
        if sum(slmpeaks) == 0
            warning("No SLM ripples for this file")
            SSslmcsd{ii} = [];
        else
            SSslmcsd{ii} = buzcode_eventCSD(lfp, SSripples{ii}.peaks(slmpeaks));
        end
        if sum(otherpeaks) == 0
            warning("No other ripples for this file")
            SSothercsd{ii} = [];
        else
            SSothercsd{ii} = buzcode_eventCSD(lfp, SSripples{ii}.peaks(otherpeaks));
        end
        SSripples{ii}.rippleStats.data.radpeaksidx = radpeaks;
        SSripples{ii}.rippleStats.data.slmpeaksidx = slmpeaks;
        SSripples{ii}.rippleStats.data.others = otherpeaks;

                % Get measures:
        filtrip = array2table(SSripples{ii}.rippleStats.maps.ripples_filtered, 'VariableNames', string(SSripples{ii}.rippleStats.maps.timestamps));
        rawrip = array2table(SSripples{ii}.rippleStats.maps.ripples_raw, 'VariableNames', string(SSripples{ii}.rippleStats.maps.timestamps));
        thisid = repmat(string(SSlfppaths(ii)), size(filtrip, 1), 1);
        expcond = repmat("SS", size(filtrip, 1), 1);
        idtable = table(thisid, expcond,'VariableNames',["ID", "expcond"]);
        thisidtable_raw = [idtable, rawrip];
        thisidtable_filt =  [idtable, filtrip];
    
        % Put it all together:
        SS_ripplesfilter = [SS_ripplesfilter; filtrip];
        SS_ripplesraw = [SS_ripplesraw; rawrip];
    
        SStableraw = [SStableraw; thisidtable_raw];
        SStablefilt = [SStablefilt; thisidtable_filt];

    catch
        SSripples{ii} = [];
    end
    close all
end

OSslmsize = nan(length(OSlfppaths), 1); OSradsize = nan(length(OSlfppaths), 1);
OSothersize = nan(length(OSlfppaths), 1);
OS_ripplesraw = [];
OS_ripplesfilter = [];
OStableraw = [];
OStablefilt = [];

for ii = 1:length(OSlfppaths)
    thisfoldername = strcat(mainFolder{2}, "\", OSlfppaths(ii));
    thisfilename = strcat(thisfoldername, "\", OSlfppaths(ii), ".hippocampalLayers.channelinfo.mat");
    load(thisfilename);
    rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
        radch = hippocampalLayers.layers{1,1}.radiatum;
    SWChannel = hippocampalLayers.layers{1,1}.radiatum; % dependerá de si son SR sink o SLM sink, no? PENDIENTE REVISAR
        slmch = hippocampalLayers.layers{1,1}.slm;

    cd(thisfoldername)
    thetafilename = strcat(thisfoldername, "\", OSlfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(thetafilename);
    tsthetaevents = thisthetaEpochs.thetaRun.newintervals;
    % For CSD data:
    lfp = getLFP('all');
    try
        OSripples{ii} = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,...
            'detector', 'customcnn','exec_env',cnnpath, 'model_file', modelpath,...
            'passband',passband,'force',true,'skipStimulationPeriods',true,...
            'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold', false,...
            'excludeIntervals', tsthetaevents);
        [OScsd{ii}, OSlfpAvg{ii}, OScsd_events{ii}, OSlfp_events{ii}] = buzcode_eacheventCSD(lfp, OSripples{ii}.peaks);
                
        close all

        tsthisripple = nan(length(OSripples{ii}.peaks), 2);
        durthisripple = tsthisripple(:, 1);
        rangethisripple = durthisripple;

        for jj = 1:size(OScsd_events{ii},3)
            OScsdsummary{ii}.tsthisripple(jj, :) = OSripples{ii}.timestamps(jj, :);
            OScsdsummary{ii}.durthisripple(jj) = diff(OScsdsummary{ii}.tsthisripple(jj, :))*1000; % paso a mseg para que round fx bien
            OScsdsummary{ii}.rangethisripple(jj) = round(OScsdsummary{ii}.durthisripple(jj)/2);
            OScsdsummary{ii}.csdtsthisrippleidx = OScsd{ii}.timestamps >= -OScsdsummary{ii}.rangethisripple(jj) & OScsd{ii}.timestamps <= OScsdsummary{ii}.rangethisripple(jj);
            OScsdsummary{ii}.csdsummary(jj, :) = mean(OScsd_events{ii}(OScsdsummary{ii}.csdtsthisrippleidx, :, jj), 1);
        end

        % 2.3. Categorize ripples according to their sink (rad or slm)
        [val, idx] = min(OScsdsummary{ii}.csdsummary, [], 2);
        slmpeaks = (idx >= radch + 1); fprintf("SLM ripples:  %4.2f ", sum(slmpeaks));
        OSslmsize(ii) = sum(slmpeaks);  
        radpeaks = idx >= rippleChannel + 1 & idx <= radch; fprintf("RAD ripples:  %4.2f ", sum(radpeaks));
        OSradsize(ii) = sum(radpeaks);
        otherpeaks = ~radpeaks & ~slmpeaks; fprintf("OTHER ripples:  %4.2f ", sum(otherpeaks));
        OSothersize(ii) = sum(otherpeaks);
        if sum(radpeaks) == 0
            warning("No RAD ripples for this file")
            OSradcsd{ii} = [];
        else
            OSradcsd{ii} = buzcode_eventCSD(lfp, OSripples{ii}.peaks(radpeaks));
        end
        if sum(slmpeaks) == 0
            warning("No SLM ripples for this file")
            OSslmcsd{ii} = [];
        else
            OSslmcsd{ii} = buzcode_eventCSD(lfp, OSripples{ii}.peaks(slmpeaks));
        end
        if sum(otherpeaks) == 0
            warning("No other ripples for this file")
            OSothercsd{ii} = [];
        else
            OSothercsd{ii} = buzcode_eventCSD(lfp, OSripples{ii}.peaks(otherpeaks));
        end
        OSripples{ii}.rippleStats.data.radpeaksidx = radpeaks;
        OSripples{ii}.rippleStats.data.slmpeaksidx = slmpeaks;
        OSripples{ii}.rippleStats.data.others = otherpeaks;
                   
        % Get measures:
        filtrip = array2table(OSripples{ii}.rippleStats.maps.ripples_filtered, 'VariableNames', string(OSripples{ii}.rippleStats.maps.timestamps));
        rawrip = array2table(OSripples{ii}.rippleStats.maps.ripples_raw, 'VariableNames', string(OSripples{ii}.rippleStats.maps.timestamps));
        thisid = repmat(string(OSlfppaths(ii)), size(filtrip, 1), 1);
        expcond = repmat("OS", size(filtrip, 1), 1);
        idtable = table(thisid, expcond,'VariableNames',["ID", "expcond"]);
        thisidtable_raw = [idtable, rawrip];
        thisidtable_filt =  [idtable, filtrip];
    
        % Put it all together:
        OS_ripplesfilter = [OS_ripplesfilter; filtrip];
        OS_ripplesraw = [OS_ripplesraw; rawrip];
    
        OStableraw = [OStableraw; thisidtable_raw];
        OStablefilt = [OStablefilt; thisidtable_filt];
    catch
        OSripples{ii} = [];
    end
end


ODslmsize = nan(length(ODlfppaths), 1); ODradsize = nan(length(ODlfppaths), 1);
ODothersize = nan(length(ODlfppaths),1);
OD_ripplesraw = [];
OD_ripplesfilter = [];
ODtableraw = [];
ODtablefilt = [];

for ii = 1:length(ODlfppaths)
    thisfoldername = strcat(mainFolder{3}, "\", ODlfppaths(ii));
    thisfilename = strcat(thisfoldername, "\", ODlfppaths(ii), ".hippocampalLayers.channelinfo.mat");
    load(thisfilename);
    radch = hippocampalLayers.layers{1,1}.radiatum;
    rippleChannel = hippocampalLayers.layers{1,1}.pyramidal;
    SWChannel = hippocampalLayers.layers{1,1}.radiatum;
        slmch = hippocampalLayers.layers{1,1}.slm;

    cd(thisfoldername)
    thetafilename = strcat(thisfoldername, "\", ODlfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(thetafilename);
    tsthetaevents = thisthetaEpochs.thetaRun.newintervals;
    lfp = getLFP('all');
    try
        ODripples{ii} = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,...
            'detector', 'customcnn','exec_env',cnnpath, 'model_file', modelpath,...
            'passband',passband,'force',true,'skipStimulationPeriods',true,...
            'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold', false,...
            'excludeIntervals', tsthetaevents);
        [ODcsd{ii}, ODlfpAvg{ii}, ODcsd_events{ii}, ODlfp_events{ii}] = buzcode_eacheventCSD(lfp, ODripples{ii}.peaks);
                
        close all
    
        tsthisripple = nan(length(ODripples{ii}.peaks), 2);
        durthisripple = tsthisripple(:, 1);
        rangethisripple = durthisripple;

        for jj = 1:size(ODcsd_events{ii},3)
            ODcsdsummary{ii}.tsthisripple(jj, :) = ODripples{ii}.timestamps(jj, :);
            ODcsdsummary{ii}.durthisripple(jj) = diff(ODcsdsummary{ii}.tsthisripple(jj, :))*1000; % paso a mseg para que round fx bien
            ODcsdsummary{ii}.rangethisripple(jj) = round(ODcsdsummary{ii}.durthisripple(jj)/2);
            ODcsdsummary{ii}.csdtsthisrippleidx = ODcsd{ii}.timestamps >= -ODcsdsummary{ii}.rangethisripple(jj) & ODcsd{ii}.timestamps <= ODcsdsummary{ii}.rangethisripple(jj);
            ODcsdsummary{ii}.csdsummary(jj, :) = mean(ODcsd_events{ii}(ODcsdsummary{ii}.csdtsthisrippleidx, :, jj), 1);
        end

        % 2.3. Categorize ripples according to their sink (rad or slm)
        [val, idx] = min(ODcsdsummary{ii}.csdsummary, [], 2);
        slmpeaks = (idx >= radch + 1); fprintf("SLM ripples:  %4.2f ", sum(slmpeaks));
        ODslmsize(ii) = sum(slmpeaks);  
        radpeaks =  idx >= rippleChannel + 1 & idx <= radch; fprintf("RAD ripples:  %4.2f ", sum(radpeaks));
        ODradsize(ii) = sum(radpeaks);
        otherpeaks = ~radpeaks & ~slmpeaks; fprintf("OTHER ripples:  %4.2f ", sum(otherpeaks));
        ODothersize(ii) = sum(otherpeaks);
        if sum(radpeaks) == 0
            warning("No RAD ripples for this file")
            ODradcsd{ii} = [];
        else
            ODradcsd{ii} = buzcode_eventCSD(lfp, ODripples{ii}.peaks(radpeaks));
        end
        if sum(slmpeaks) == 0
            warning("No SLM ripples for this file")
            ODslmcsd{ii} = [];
        else
            ODslmcsd{ii} = buzcode_eventCSD(lfp, ODripples{ii}.peaks(slmpeaks));
        end
        if sum(otherpeaks) == 0
            warning("No other ripples for this file")
            ODothercsd{ii} = [];
        else
            ODothercsd{ii} = buzcode_eventCSD(lfp, ODripples{ii}.peaks(otherpeaks));
        end

        ODripples{ii}.rippleStats.data.radpeaksidx = radpeaks;
        ODripples{ii}.rippleStats.data.slmpeaksidx = slmpeaks;
        ODripples{ii}.rippleStats.data.others = otherpeaks;

            % Get measures:
        filtrip = array2table(ODripples{ii}.rippleStats.maps.ripples_filtered, 'VariableNames', string(ODripples{ii}.rippleStats.maps.timestamps));
        rawrip = array2table(ODripples{ii}.rippleStats.maps.ripples_raw, 'VariableNames', string(ODripples{ii}.rippleStats.maps.timestamps));
        thisid = repmat(string(ODlfppaths(ii)), size(filtrip, 1), 1);
        expcond = repmat("OD", size(filtrip, 1), 1);
        idtable = table(thisid, expcond,'VariableNames',["ID", "expcond"]);
        thisidtable_raw = [idtable, rawrip];
        thisidtable_filt =  [idtable, filtrip];
    
        % Put it all together:
        OD_ripplesfilter = [OD_ripplesfilter; filtrip];
        OD_ripplesraw = [OD_ripplesraw; rawrip];
    
        ODtableraw = [ODtableraw; thisidtable_raw];
        ODtablefilt = [ODtablefilt; thisidtable_filt];
    catch
        ODripples{ii} = [];
    end
    close all
end

% All have the same structure: nameoffolder.ripples.events.mat
% filename lfppaths+.ripples.events.mat
% foldername mainfolder+lfppaths
allripplesfile = strcat(ripplefolder, "\allrippledata.mat");
save(allripplesfile, "SSripples", "OSripples", "ODripples", '-mat')

%% 3. Measures table

%% SHAM-SHAM DATA:
SS_ripplemeasures = [];
for ii = 1:length(SSlfppaths)

    if ~isempty(SSripples{ii})
        ripples = SSripples{ii};
        % Get measures:
        peaknormpow = ripples.peakNormedPower;
        peakfq = ripples.rippleStats.data.peakFrequency;
        peakamp = ripples.rippleStats.data.peakAmplitude;
        duration = ripples.rippleStats.data.duration;
        intrip = ripples.rippleStats.data.interRippleFrequency;
        specEnt = ripples.rippleStats.data.spectralEntropy.';
        fastripidx = ripples.rippleStats.data.fastRippleIndex.';
        multitap = ripples.rippleStats.data.multiTapperFreq;
        % New measures
        rippleIFstat = ripples.rippleStats.data.rippleIFstat;
        rippleIF_meanFreq = ripples.rippleStats.data.rippleIF.meanFreq;
        rippleIF_sdFreq = ripples.rippleStats.data.rippleIF.stdFreq;
        rippleIF_drift = ripples.rippleStats.data.rippleIF.drift;
        % rippleIF_ridge = ripples.rippleStats.data.rippleIF.ridge.';
        % Rad vs slm sink categorization
        thisfilessink = strings(length(peaknormpow), 1);
        radidx = ripples.rippleStats.data.radpeaksidx;
        slmidx = ripples.rippleStats.data.slmpeaksidx;
        otheridx = ripples.rippleStats.data.others;
        thisfilessink(radidx) = "rad";
        thisfilessink(slmidx) = "slm";
        thisfilessink(otheridx) = "other";
        thisID = repmat(string(SSlfppaths{ii}), length(peaknormpow), 1);
        expcond = repmat("SS", length(peaknormpow), 1);
    
        % Put it all together:
        thisrectable = table(thisID,expcond, peaknormpow, peakfq, peakamp, duration, intrip, ...
            specEnt, fastripidx, multitap,...
            rippleIFstat, rippleIF_meanFreq, rippleIF_sdFreq, rippleIF_drift,...
            thisfilessink, VariableNames = ["ID","expcond","PeakNormPow", "PeakFq", "PeakAmp", ...
            "Duration", "interrippleFq", "spectralEntropy", "FastRippleIdx", "multitap", ...
            "freqstat", "meanfq_nonstat", "std_nonstat", "drift_nonstat", "sinktype"]);
        SS_ripplemeasures = [SS_ripplemeasures; thisrectable];
    else
        continue
    end
end

% Remove events that are likely noise (>50 Peak Norm Power) 
ssfilenameoriginal = strcat(ripplefolder, "\SS_rippledata_original.xlsx");
writetable(SS_ripplemeasures,ssfilenameoriginal)
SSkeepidx = SS_ripplemeasures.PeakNormPow < 50;
SS_ripples_keep = SS_ripplemeasures(SSkeepidx, :); 
% Save xlsx:

ssfilename = strcat(ripplefolder, "\SS_rippledata_keep.xlsx");
writetable(SS_ripples_keep,ssfilename)

%% ORX-SHAM DATA:

OS_ripplemeasures = [];

for ii = 1:length(OSlfppaths)
    
    if ~isempty(OSripples{ii})
        ripples = OSripples{ii};
        % Get measures:
        peaknormpow = ripples.peakNormedPower;
        peakfq = ripples.rippleStats.data.peakFrequency;
        peakamp = ripples.rippleStats.data.peakAmplitude;
        duration = ripples.rippleStats.data.duration;
        intrip = ripples.rippleStats.data.interRippleFrequency;
        specEnt = ripples.rippleStats.data.spectralEntropy.';
        fastripidx = ripples.rippleStats.data.fastRippleIndex.';
        multitap = ripples.rippleStats.data.multiTapperFreq;
        % New measures
        rippleIFstat = ripples.rippleStats.data.rippleIFstat;
        rippleIF_meanFreq = ripples.rippleStats.data.rippleIF.meanFreq;
        rippleIF_sdFreq = ripples.rippleStats.data.rippleIF.stdFreq;
        rippleIF_drift = ripples.rippleStats.data.rippleIF.drift;
        % rippleIF_ridge = ripples.rippleStats.data.rippleIF.ridge.';
        % Rad vs slm sink categorization
        thisfilessink = strings(length(peaknormpow), 1);
        radidx = ripples.rippleStats.data.radpeaksidx;
        slmidx = ripples.rippleStats.data.slmpeaksidx;
        otheridx = ripples.rippleStats.data.others;
        thisfilessink(radidx) = "rad";
        thisfilessink(slmidx) = "slm";
        thisfilessink(otheridx) = "other";
        thisID = repmat(string(OSlfppaths{ii}), length(peaknormpow), 1);
        expcond = repmat("OS", length(peaknormpow), 1);
    
        % Put it all together:
        thisrectable = table(thisID,expcond, peaknormpow, peakfq, peakamp, duration, intrip, ...
            specEnt, fastripidx, multitap,...
            rippleIFstat, rippleIF_meanFreq, rippleIF_sdFreq, rippleIF_drift,...
            thisfilessink, VariableNames = ["ID","expcond","PeakNormPow", "PeakFq", "PeakAmp", ...
            "Duration", "interrippleFq", "spectralEntropy", "FastRippleIdx", "multitap", ...
            "freqstat", "meanfq_nonstat", "std_nonstat", "drift_nonstat", "sinktype"]);
        OS_ripplemeasures = [OS_ripplemeasures; thisrectable];
    else
        continue
    end
end

osfilenameoriginal = strcat(ripplefolder, "\OS_rippledata_original.xlsx");
writetable(OS_ripplemeasures,osfilenameoriginal)
OSkeepidx = OS_ripplemeasures.PeakNormPow < 50;
OS_ripples_keep = OS_ripplemeasures(OSkeepidx, :); 
osfilename = strcat(ripplefolder, "\OS_rippledata.xlsx");
writetable(OS_ripples_keep,osfilename)


%% ORX-DHT DATA:

OD_ripplemeasures = [];

for ii = 1:length(ODlfppaths)
    if ~isempty(ODripples{ii})
        ripples = ODripples{ii};
        % Get measures:
        peaknormpow = ripples.peakNormedPower;
        peakfq = ripples.rippleStats.data.peakFrequency;
        peakamp = ripples.rippleStats.data.peakAmplitude;
        duration = ripples.rippleStats.data.duration;
        intrip = ripples.rippleStats.data.interRippleFrequency;
        specEnt = ripples.rippleStats.data.spectralEntropy.';
        fastripidx = ripples.rippleStats.data.fastRippleIndex.';
        multitap = ripples.rippleStats.data.multiTapperFreq;
        % New measures
        rippleIFstat = ripples.rippleStats.data.rippleIFstat;
        rippleIF_meanFreq = ripples.rippleStats.data.rippleIF.meanFreq;
        rippleIF_sdFreq = ripples.rippleStats.data.rippleIF.stdFreq;
        rippleIF_drift = ripples.rippleStats.data.rippleIF.drift;
        % rippleIF_ridge = ripples.rippleStats.data.rippleIF.ridge.';
        % Rad vs slm sink categorization
        thisfilessink = strings(length(peaknormpow), 1);
        radidx = ripples.rippleStats.data.radpeaksidx;
        slmidx = ripples.rippleStats.data.slmpeaksidx;
        otheridx = ripples.rippleStats.data.others;
        thisfilessink(radidx) = "rad";
        thisfilessink(slmidx) = "slm";
        thisfilessink(otheridx) = "other";
        thisID = repmat(string(ODlfppaths{ii}), length(peaknormpow), 1);
        expcond = repmat("OD", length(peaknormpow), 1);
    
        % Put it all together:
        thisrectable = table(thisID,expcond, peaknormpow, peakfq, peakamp, duration, intrip, ...
            specEnt, fastripidx, multitap,...
            rippleIFstat, rippleIF_meanFreq, rippleIF_sdFreq, rippleIF_drift,...
            thisfilessink, VariableNames = ["ID","expcond","PeakNormPow", "PeakFq", "PeakAmp", ...
            "Duration", "interrippleFq", "spectralEntropy", "FastRippleIdx", "multitap", ...
            "freqstat", "meanfq_nonstat", "std_nonstat", "drift_nonstat", "sinktype"]);
        OD_ripplemeasures = [OD_ripplemeasures; thisrectable];
    else
        continue
    end
end


ODkeepidx = OD_ripplemeasures.PeakNormPow < 50;
OD_ripples_keep = OD_ripplemeasures(ODkeepidx, :); 

odfilenameoriginal = strcat(ripplefolder, "\OD_rippledata_original.xlsx");
writetable(OD_ripplemeasures,odfilenameoriginal)
odfilename = strcat(ripplefolder, "\OD_rippledata_keep.xlsx");
writetable(OD_ripplemeasures,odfilename)


% JOIN ALL DATA AND SAVE:

all_ripples = [SS_ripplemeasures; OS_ripplemeasures; OD_ripplemeasures];
all_ripples_keep = [SS_ripples_keep; OS_ripples_keep; OD_ripples_keep];

% savedatafile = uigetdir();
allripplesfile = strcat(ripplefolder, "\SexHormones_Males_ripples_keep.xlsx");
writetable(all_ripples_keep, allripplesfile)
writetable(all_ripples, strcat(ripplefolder, "\SexHormones_Males_ripples_original.xlsx"))

SSid = repmat("SS", size(SS_ripplesraw, 1), 1);
OSid = repmat("OS", size(OS_ripplesraw, 1), 1);
ODid = repmat("OD", size(OD_ripplesraw, 1), 1);
rawripples = [SS_ripplesraw; OS_ripplesraw; OD_ripplesraw];
filteredripples = [SS_ripplesfilter; OS_ripplesfilter; OD_ripplesfilter];
id = [SSid; OSid; ODid];
waveformfile = strcat(ripplefolder, "\ripplewaveform.mat");
save(waveformfile, "rawripples", "filteredripples", '-mat')
% 


%% 4. ADD MOUSE ID TO THE TABLE

SS = ["C089", "C090", "C091", "C199", "E075"];
OS = ["C093", "C094", "C095", "D094"];
OD = ["D097", "D096", "E078", "E079"];
all_ripples = SexHormonesMalesripplesoriginal;

subject = repmat("", size(all_ripples, 1), 1);
all_ripples.ID = string(all_ripples.ID);
for ii = 1:length(SS)
    thissubject = contains(all_ripples.ID, SS(ii));
    subject(thissubject) = SS(ii);
end

for ii = 1:length(OS)
    thissubject = contains(all_ripples.ID, OS(ii));
    subject(thissubject) = OS(ii);
end

for ii = 1:length(OD)
    thissubject = contains(all_ripples.ID, OD(ii));
    subject(thissubject) = OD(ii);
end

all_ripples.subject = subject;

%% 5. Remove outliers before analysis
% [dataTable, mouseCol, recCol] = deal(all_ripples, 'subject', 'ID');
% [out, flagTable] = detect_ripple_outliers(dataTable, mouseCol, recCol);
% 
% % variables you checked
% vars = {'PeakAmp','PeakNormPow','Duration','PeakFq'};
% 
% % Per-variable indices (union flags)
% for v = 1:numel(vars)
%     var = vars{v};
%     if isfield(out.union, var)
%         idxVar = out.union.(var);           % row indices for this variable
%     else
%         idxVar = [];
%     end
%     fprintf('%s: %d flagged rows\n', var, numel(idxVar));
%     % optionally: display first 10 indices
%     if ~isempty(idxVar), fprintf('  first indices: %s\n', mat2str(idxVar(1:min(10,end)))); end
% end
% 
% % Single vector of indices to discard (union across variables)
% discardIdx = find(out.unionIdx);          % row indices to remove
% fprintf('TOTAL unique rows to discard (union across vars): %d\n', numel(discardIdx));
% 
% rowidx = sum((1:size(all_ripples, 1)) == discardIdx, 1).';
% filteredRipples = table2array(filteredripples);
% noout_filteredripples = filteredRipples(~rowidx, :);
% filteredRipples = noout_filteredripples;
% all_ripples = all_ripples(~rowidx, :);

%% 5. Hierarchical robust outlier removal (per-recording -> per-subject -> per-group -> global)
[dataTable, mouseCol, recCol] = deal(all_ripples, 'subject', 'ID');

% user parameters (adjust if needed)
vars = {'PeakAmp','PeakNormPow','Duration','PeakFq'};
logVars = {'PeakAmp','PeakNormPow','Duration'};   % variables to log-transform before MAD
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
rowidx = false(height(all_ripples),1);
rowidx(discardIdx) = true;
all_ripples = all_ripples(~rowidx, :);
filteredRipples = all_ripples;  % keep same name semantics as your code



%% 6. DESCRIPTIVE ANALYSIS PER GROUP:

% variables to plot
vars = {'PeakFq','PeakAmp','spectralEntropy','Duration'};

% desired group order and colors: SS (pink), OS (light blue), OD (dark blue)
groupLabels = {'SS','OS','OD'};
cols = [1.00 0.55 0.75; 0.65 0.85 1.00; 0.05 0.25 0.65]; % pink, light blue, dark blue
alphaVal = 0.35; % more transparent
nbins = 30;

% figure
fig = figure('Units','normalized','Position',[0.05 0.15 0.9 0.6],'Color','w');

for i = 1:numel(vars)
    ax = subplot(1,4,i);
    hold(ax,'on');
    varname = vars{i};
    
    % combined data for consistent bin edges
    combined = all_ripples.(varname);
    combined = combined(~isnan(combined) & isfinite(combined));
    if isempty(combined)
        title(ax, strrep(varname,'_',' '));
        continue;
    end
    edges = linspace(min(combined), max(combined), nbins+1);
    
    plotted = gobjects(0);
    legLabels = {};
    for g = 1:numel(groupLabels)
        grp = groupLabels{g};
        idx = categorical(all_ripples.expcond) == grp;
        grpdata = all_ripples.(varname)(idx);
        grpdata = grpdata(~isnan(grpdata) & isfinite(grpdata));
        if isempty(grpdata)
            continue;
        end
        h = histogram(ax, grpdata, edges, 'Normalization','probability', ...
            'FaceAlpha', alphaVal, 'EdgeColor','none', 'FaceColor', cols(g,:));
        uistack(h,'top');
        plotted(end+1) = h; %#ok<SAGROW>
        legLabels{end+1} = grp; %#ok<SAGROW>
    end
    
    title(ax, strrep(varname,'_',' '), 'Interpreter','none');
    xlabel(ax,'Value');
    if i == 1
        ylabel(ax,'Probability');
    end
    if ~isempty(plotted)
        legend(ax, plotted, legLabels, 'Location','best');
    end
    set(ax,'FontSize',11,'Box','off');
    hold(ax,'off');
    hold off
end

drawnow;
set(fig,'PaperPositionMode','auto');
hold off
% save as SVG (vector)
print(fig,'-dsvg','all_ripples_histograms.svg');



%% /. LINEAR-MIXED EFFECT ANALYSIS BASED ON PREVIOUSLY TESTED STATISTICAL DISTRIBUTIONS

% 6.1. One-way ANOVA (for all ripples)

data = all_ripples;

% Now analysis are performed excluding outliers:

% 1. Peak Normalized Power (average power over ripple window)
all_peaknormpow = analyze_ripple_measure_lme_figures(data, 'PeakNormPow');

% 2. Peak Ampitude (instantaneous)
all_peakamp = analyze_ripple_measure_lme_figures(data, 'PeakAmp');

% 3. Peak Frequency
all_peakfreq = analyze_ripple_measure_lme_figures(data, 'PeakFq');

% 4. Duration
all_duration = analyze_ripple_measure_glme(data, 'Duration', 'gamma');

% 5. Spectral Entropy
all_entropy = analyze_ripple_measure_lme(data, 'spectralEntropy');

% 6. Interriple Frequency
all_interriplefq = analyze_ripple_measure_glme(data, 'interrippleFq', 'weibull');

ID = unique(data.ID);

% 6.2. Two-way ANOVA (for future use to add more sink type ripples) 

% savepath = uigetdir();
% all_peaknormpow_2w = analyze_ripple_measure_lme_2way_mod(data, 'PeakNormPow', 'save_path', savepath, 'fig_name', '2way_PeakNormPow');
% all_amp_2w = analyze_ripple_measure_lme_2way_mod(data, 'PeakAmp', 'save_path', savepath, 'fig_name', '2_wayPeakAmp');
% all_peakfreq_2w = analyze_ripple_measure_lme_2way_mod(data, 'PeakFq', 'save_path', savepath, 'fig_name', '2_wayPeakFq');
% all_specEnt_2w = analyze_ripple_measure_lme_2way_mod2(data, 'spectralEntropy', 'save_path', savepath, 'fig_name', '2_wayspectralEntropy');
% all_duration_2w = analyze_ripple_measure_lme_2way_mod2(data, 'Duration', 'save_path', savepath, 'fig_name', '2_wayduration');
% all_interriplefq_2w = analyze_ripple_measure_lme_2way_mod2(data, 'interrippleFq', 'save_path', savepath, 'fig_name', '2_wayinterriplefq');

%% 8. (PART 1) General UMAP (with passband filtered LFP data)

% 8.1. Take the middle window of 127 samples (aprox 100 ms)

if size(filteredRipples, 2) == 187
    midpoint = round(size(filteredRipples, 2)/2);
    filteredRipples_centered = filteredRipples(:, midpoint-63:midpoint+63); % 127 samples
else
    filteredRipples_centered = filteredRipples;
end


% 8.2. Apply UMAP with optimal dimensionality

[reduction, umap] = run_umap(filteredRipples_centered, ...
    'n_components', 4, 'metric', 'euclidean', 'min_dist', 0.1);


% 8.3. Figure per groups

% Extract group labels from data table
if istable(data) && ismember('expcond', data.Properties.VariableNames)
    groupLabels = data.expcond;
elseif isstruct(data) && isfield(data, 'expcond')
    groupLabels = data.expcond;
else
    error('Variable "expcond" not found in data. Please ensure data contains group labels.');
end

groupLabels = categorical(groupLabels, {'SS','OS','OD'}, 'Ordinal', true);
groups = unique(groupLabels);
% groups = categorical(groups, {'SS','OS','OD'}, 'Ordinal', true);
colors = lines(numel(groups));

figure(1); hold on;
for i = 1:numel(groups)
    idx = groupLabels == groups(i);
    scatter(reduction(idx,1), reduction(idx,2), 10, colors(i,:), 'filled','MarkerFaceAlpha', 0.5);
end
xlabel('UMAP1');
ylabel('UMAP2');
title('UMAP Embedding of Sharp-Wave Ripple Measures');
legend(groups, 'Location', 'best');
% grid on;
hold off;


% 8.4. Figure per measures

UMAP_POS = [0.003472222222222, ...
            0.531797235023041, ...
            0.997685185185185, ...
            0.326267281105991];
centralPct = 95;
vars = {'PeakFq','PeakAmp','spectralEntropy','Duration'};

% Diverging colormap
nC = 256;
stops = [
    0.03 0.18 0.42;
    0.16 0.38 0.75;
    0.60 0.80 0.90;
    1.00 1.00 1.00;
    0.99 0.73 0.53;
    0.94 0.29 0.20;
    0.50 0.00 0.00
];
cmap = interp1(linspace(0,1,size(stops,1)), stops, linspace(0,1,nC));

figure('Units','normalized','Position',UMAP_POS);
for i = 1:numel(vars)
    vals = double(data.(vars{i}));
    lowerP = (100-centralPct)/2;
    clim = prctile(vals, [lowerP 100-lowerP]);
    vals_c = min(max(vals,clim(1)),clim(2));

    ax = subplot(1,4,i);
    scatter(reduction(:,1), reduction(:,2), 20, vals_c, 'filled', ...
        'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
    % axis equal tight;
    colormap(ax, cmap);
    caxis(clim);
    colorbar;
    title(vars{i}, 'Interpreter','none');
    xlabel('UMAP1'); ylabel('UMAP2');
end
sgtitle(sprintf('UMAP colored by ripple measures (central %d%%)', centralPct));


% 8.5.  Measures per group

greyColor = [0.85 0.85 0.85];

for g = 1:numel(groups)
    idx = groupLabels == groups(g);

    figure('Units','normalized','Position',UMAP_POS);
    for i = 1:numel(vars)
        vals = double(data.(vars{i}));
        clim = prctile(vals, [(100-centralPct)/2, 100-(100-centralPct)/2]);
        vals_c = min(max(vals,clim(1)),clim(2));

        ax = subplot(1,4,i);
        scatter(reduction(~idx,1), reduction(~idx,2), 20, greyColor, ...
            'filled','MarkerFaceAlpha',0.3,'MarkerEdgeColor','none');
        hold on;
        scatter(reduction(idx,1), reduction(idx,2), 20, vals_c(idx), ...
            'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
        hold off;

        % axis equal tight;
        colormap(ax,cmap);
        caxis(clim);
        colorbar;
        title(vars{i}, 'Interpreter','none');
        xlabel('UMAP1'); ylabel('UMAP2');
    end
    sgtitle(sprintf('%s (n=%d)', char(groups(g)), sum(idx)));
end


savePath = pwd;

figs = findobj('Type','figure');
for i = 1:numel(figs)
    f = figs(i);
    name = sprintf('Figure_%d', f.Number);
    saveas(f, fullfile(savePath, [name '.svg']));
end

% UMAP Visualization by Group

% Define colors for groups (SS: light pink, OS: light blue, OD: dark blue)
cols = [1.00 0.75 0.80;   % SS - light pink
        0.65 0.85 1.00;   % OS - light blue
        0.05 0.25 0.65];  % OD - dark blue

groups = {'SS', 'OS', 'OD'};
invertedgroups = { 'OD', 'OS','SS'};

figure('Position', [100 100 800 600]);
hold on;

for g = 1:3
    idx = all_ripples.expcond == invertedgroups{g};
    scatter(reduction(idx, 1), reduction(idx, 2), 10, cols(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
end

xlabel('UMAP 1');
ylabel('UMAP 2');
title('UMAP Embedding by Experimental Condition');
legend(groups, 'Location', 'best');
box on;
hold off;

% 8.6. Save figures

savePath = pwd;

figs = findobj('Type','figure');
for i = 1:numel(figs)
    f = figs(i);
    name = sprintf('Figure_%d', f.Number);
    saveas(f, fullfile(savePath, [name '.svg']));
end

% save(fullfile(savePath,'UMAP_global.mat'), ...
%     'reduction','coeff','score','explained','optimal_dim','filteredRipples_centered');
% 
% fprintf('All figures and data saved to: %s\n', savePath);

%% Repeat figures but with clim constrained by LME results:




% 8.3. Figure per groups

% Extract group labels from data table
if istable(data) && ismember('expcond', data.Properties.VariableNames)
    groupLabels = data.expcond;
elseif isstruct(data) && isfield(data, 'expcond')
    groupLabels = data.expcond;
else
    error('Variable "expcond" not found in data. Please ensure data contains group labels.');
end

groupLabels = categorical(groupLabels, {'SS','OS','OD'}, 'Ordinal', true);
groups = unique(groupLabels);
% groups = categorical(groups, {'SS','OS','OD'}, 'Ordinal', true);
colors = lines(numel(groups));

figure(1); hold on;
for i = 1:numel(groups)
    idx = groupLabels == groups(i);
    scatter(reduction(idx,1), reduction(idx,2), 10, colors(i,:), 'filled','MarkerFaceAlpha', 0.5);
end
xlabel('UMAP1');
ylabel('UMAP2');
title('UMAP Embedding of Sharp-Wave Ripple Measures');
legend(groups, 'Location', 'best');
% grid on;
hold off;


% 8.4. Figure per measures

UMAP_POS = [0.003472222222222, ...
            0.531797235023041, ...
            0.997685185185185, ...
            0.326267281105991];
centralPct = 95;
vars = {'PeakFq','PeakAmp','spectralEntropy','Duration'};

% Diverging colormap
nC = 256;
stops = [
    0.03 0.18 0.42;
    0.16 0.38 0.75;
    0.60 0.80 0.90;
    1.00 1.00 1.00;
    0.99 0.73 0.53;
    0.94 0.29 0.20;
    0.50 0.00 0.00
];
cmap = interp1(linspace(0,1,size(stops,1)), stops, linspace(0,1,nC));


figure('Units','normalized','Position',UMAP_POS);
for i = 1:numel(vars)
    vals = double(data.(vars{i}));
    lowerP = (100-centralPct)/2;
    % if  strcmp(vars{i}, 'spectralEntropy')
    %     clim = [4.9, 5.4];
    % elseif strcmp(vars{i}, 'Duration')
    %     clim = [0.035, 0.06];
    % else 
    %     clim = prctile(vals, [lowerP 100-lowerP]);
    % end
    clim = prctile(vals, [lowerP 100-lowerP]);
    vals_c = min(max(vals,clim(1)),clim(2));

    ax = subplot(1,4,i);
    scatter(reduction(:,1), reduction(:,2), 20, vals_c, 'filled', ...
        'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
    % axis equal tight;
    colormap(ax, cmap);
    caxis(clim);
    colorbar;
    title(vars{i}, 'Interpreter','none');
    xlabel('UMAP1'); ylabel('UMAP2');
end
sgtitle(sprintf('UMAP colored by ripple measures (central %d%%)', centralPct));


% 8.5.  Measures per group

greyColor = [0.85 0.85 0.85];

for g = 1:numel(groups)
    idx = groupLabels == groups(g);

    figure('Units','normalized','Position',UMAP_POS);
    for i = 1:numel(vars)
        vals = double(data.(vars{i}));
        clim = prctile(vals, [(100-centralPct)/2, 100-(100-centralPct)/2]);
        vals_c = min(max(vals,clim(1)),clim(2));

        ax = subplot(1,4,i);
        scatter(reduction(~idx,1), reduction(~idx,2), 20, greyColor, ...
            'filled','MarkerFaceAlpha',0.3,'MarkerEdgeColor','none');
        hold on;
        scatter(reduction(idx,1), reduction(idx,2), 20, vals_c(idx), ...
            'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
        hold off;

        % axis equal tight;
        colormap(ax,cmap);
        caxis(clim);
        colorbar;
        title(vars{i}, 'Interpreter','none');
        xlabel('UMAP1'); ylabel('UMAP2');
    end
    sgtitle(sprintf('%s (n=%d)', char(groups(g)), sum(idx)));
end


savePath = pwd;

figs = findobj('Type','figure');
for i = 1:numel(figs)
    f = figs(i);
    name = sprintf('Figure_%d', f.Number);
    saveas(f, fullfile(savePath, [name '.svg']));
end

% UMAP Visualization by Group

% Define colors for groups (SS: light pink, OS: light blue, OD: dark blue)
% cols = [;   % SS - light pink
%         0.65 0.85 1.00;   % OS - light blue
%         0.05 0.25 0.65];  % OD - dark blue

% I invert the onder in which they are plotted because of  visualization
% convenience
cols = [0.05 0.25 0.65;   % OD - light pink
        0.65 0.85 1.00;   % OS - light blue
        1.00 0.75 0.80];  % SS - dark blue

% groups = {'SS', 'OS', 'OD'};
groups = { 'OD', 'OS','SS'};

figure('Position', [100 100 800 600]);
hold on;

for g = 1:3
    idx = all_ripples.expcond == groups{g};
    scatter(reduction(idx, 1), reduction(idx, 2), 10, cols(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
end

xlabel('UMAP 1');
ylabel('UMAP 2');
title('UMAP Embedding by Experimental Condition');
legend(groups, 'Location', 'best');
box on;
hold off;
%% 8.7. % Compute centroids for each recording
%% --- Parameters / sanity checks ---
assert(exist('reduction','var')==1 && exist('all_ripples','var')==1, ...
    'Variables reduction and all_ripples must exist in workspace.');
nComps = size(reduction,2);
assert(height(all_ripples)==size(reduction,1), 'Row count mismatch: all_ripples vs reduction');

groups = {'SS','OS','OD'};                % desired order
cols = [0.2 0.6 0.9; 0.9 0.5 0.2; 0.3 0.7 0.4]; % colors for SS,OS,OD

% ensure expcond categorical with desired ordering
all_ripples.expcond = categorical(all_ripples.expcond, groups);

%% --- Compute recording centroids (vector centroid per recording) --- 
[recordings, ~, recIdx] = unique(all_ripples.ID, 'stable'); 
nRecs = numel(recordings);
centroids = zeros(nRecs, nComps);
for r = 1:nRecs
    centroids(r,:) = mean(reduction(recIdx==r, :), 1, 'omitnan');
end

% metadata for each recording: take first occurrence per recording
firstIdx = accumarray(recIdx, (1:length(recIdx))', [], @min);
rec_subject = all_ripples.subject(firstIdx);
rec_expcond = categorical(all_ripples.expcond(firstIdx), groups);

% build centroid table
centroid_table = table(recordings, rec_subject, rec_expcond, 'VariableNames', {'ID','subject','expcond'});
for k = 1:nComps
    centroid_table.(['UMAP' num2str(k)]) = centroids(:,k);
end
centroid_table.subject = categorical(centroid_table.subject);
centroid_table.ID = categorical(centroid_table.ID);
centroid_table.expcond = categorical(centroid_table.expcond, groups);
centroid_table.centroid_norm = sqrt(sum(centroids.^2,2));    % scalar summary



%% --- Run LME for each UMAP component + centroid_norm ---
variables = [arrayfun(@(k) sprintf('UMAP%d',k), 1:nComps, 'Uni',0), {'centroid_norm'}];
nVars = numel(variables);
pMainVals = nan(nVars,1);
pairP = nan(nVars,3);

for v = 1:nVars
    [pMainVals(v), pairP(v,:)] = lme_main_and_pairs(centroid_table, variables{v});
end

% display concise table
T = table(variables(:), pMainVals, pairP(:,1), pairP(:,2), pairP(:,3), ...
    'VariableNames', {'Variable','p_main','p_OS_vs_SS','p_OD_vs_SS','p_OD_vs_OS'});
disp('LME results (Bonferroni-corrected pairwise p-values):'); disp(T);

%% --- UMAP visualization: ripples faint + recording centroids + group centroids ---
if nComps >= 2
    figure('Position',[120 120 1000 480]);
    hold on;
    cats_rip = categories(all_ripples.expcond);
    for g = 1:numel(groups)
        thisCat = groups{g};
        idx = all_ripples.expcond == thisCat;
        scatter(reduction(idx,1), reduction(idx,2), 10, cols(g,:), 'filled',...
            'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0.06);
    end
    % recording centroids
    for g = 1:numel(groups)
        thisCat = groups{g};
        idx = centroid_table.expcond == thisCat;
        scatter(centroids(idx,1), centroids(idx,2), 90, 'filled', 'MarkerFaceColor', cols(g,:), ...
            'MarkerEdgeColor','k','LineWidth',0.8);
    end
    % group centroids (mean of recording centroids)
    for g = 1:numel(groups)
        idx = centroid_table.expcond == groups{g};
        if any(idx)
            gm = mean(centroids(idx,1:2),1,'omitnan');
            scatter(gm(1), gm(2), 300, 'p', 'MarkerFaceColor', cols(g,:), 'MarkerEdgeColor','k','LineWidth',1.2);
        end
    end
    axis equal; box on;
    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('UMAP1 vs UMAP2 — ripples (faint) + recording centroids + group centroids');
    legend(groups,'Location','best');
    hold off;
end

%% --- Boxplots (first 4 UMAP dims) with recordings->subject lines ---

%% --- LME Analysis and Boxplots for UMAP Centroids ---

% Define colors for groups (SS: light pink, OS: light blue, OD: dark blue)
cols = [1.00 0.75 0.80;   % SS - light pink
        0.65 0.85 1.00;   % OS - light blue
        0.05 0.25 0.65];  % OD - dark blue

groups = {'SS', 'OS', 'OD'};

% Run LME and collect results for each UMAP dimension
fprintf('\n========================================\n');
fprintf('LME ANALYSIS: UMAP CENTROID COMPONENTS\n');
fprintf('========================================\n');

pMainVals = zeros(nComps, 1);
pairP = zeros(nComps, 3);

for comp = 1:nComps
    varname = sprintf('UMAP%d', comp);
    centroid_table.y = centroid_table.(varname);
    
    % Fit LME model
    lme = fitlme(centroid_table, 'y ~ expcond + (1|subject)');
    
    % Main effect test
    pMainVals(comp) = coefTest(lme, [0 1 0; 0 0 1]);
    
    % Pairwise comparisons
    pairP(comp, 1) = coefTest(lme, [0 1 0]);      % OS vs SS
    pairP(comp, 2) = coefTest(lme, [0 0 1]);      % OD vs SS
    pairP(comp, 3) = coefTest(lme, [0 -1 1]);     % OD vs OS
    
    % Print results
    fprintf('\n%s:\n', varname);
    fprintf('  Main effect p-value: %.4f\n', pMainVals(comp));
    fprintf('  Pairwise comparisons:\n');
    fprintf('    OS vs SS: p = %.4f\n', pairP(comp, 1));
    fprintf('    OD vs SS: p = %.4f\n', pairP(comp, 2));
    fprintf('    OD vs OS: p = %.4f\n', pairP(comp, 3));
end
fprintf('========================================\n\n');

% Create plots
nPlot = min(4, nComps);
figure('Position', [100 100 1100 800]);

for comp = 1:nPlot
    subplot(1, 4, comp); hold on;
    varname = sprintf('UMAP%d', comp);
    y = centroid_table.(varname);
    
    x_positions = [1 2 3];
    
    for g = 1:3
        grp = groups{g};
        idx_grp = centroid_table.expcond == grp;
        subj_list = unique(centroid_table.subject(idx_grp));
        x_pos = x_positions(g);
        
        subj_means = zeros(numel(subj_list), 1);
        for s = 1:numel(subj_list)
            idx_subj = idx_grp & centroid_table.subject == subj_list(s);
            subj_means(s) = mean(y(idx_subj), 'omitnan');
            
            rec_vals = y(idx_subj);
            x_jitter = x_pos + 0.12 * (rand(numel(rec_vals), 1) - 0.5);
            
            % Plot recordings
            scatter(x_jitter, rec_vals, 30, cols(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
            
            % Connect recordings to subject mean
            for r = 1:numel(rec_vals)
                plot([x_jitter(r), x_pos], [rec_vals(r), subj_means(s)], ...
                    'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            end
        end
        
        % Plot subject means
        scatter(repmat(x_pos, numel(subj_means), 1), subj_means, 100, 'k', 'filled');
        
        % Boxplot overlay
        boxchart(x_pos * ones(numel(subj_means), 1), subj_means, ...
            'BoxFaceColor', cols(g, :), 'BoxFaceAlpha', 0.3, 'MarkerStyle', 'none', 'LineWidth', 1.2);
    end
    
    xlim([0.5 3.5]);
    xticks(x_positions);
    xticklabels(groups);
    ylabel(varname);
    title(sprintf('%s (p = %.4f)', varname, pMainVals(comp)));
    set(gca, 'FontSize', 11, 'LineWidth', 1);
    box on;
    hold off;
end

sgtitle('UMAP Component Centroids by Condition');
%% UMAP Dimension Interpretation - Ripple Characteristics

% Define colors for groups (SS: light pink, OS: light blue, OD: dark blue)
cols = [1.00 0.75 0.80;   % SS - light pink
        0.65 0.85 1.00;   % OS - light blue
        0.05 0.25 0.65];  % OD - dark blue

groups = {'SS', 'OS', 'OD'};

% Key ripple measures
measures = {'PeakAmp', 'PeakFq', 'spectralEntropy', 'Duration'};
nMeas = length(measures);

% Correlate UMAP dimensions with ripple measures
meas_corr = zeros(nComps, nMeas);
for i = 1:nMeas
    meas_corr(:, i) = corr(reduction, double(all_ripples.(measures{i})), 'rows', 'complete');
end

% Visualize
figure('Position', [100 100 800 600]);
imagesc(meas_corr');
colormap([linspace(0.05, 1, 128)', linspace(0.25, 1, 128)', linspace(0.65, 1, 128)';
          linspace(1, 0.5, 128)', linspace(1, 0, 128)', linspace(1, 0, 128)']);
colorbar;
caxis([-1 1]);
set(gca, 'XTick', 1:nComps, 'XTickLabel', arrayfun(@(x) sprintf('UMAP%d', x), 1:nComps, 'UniformOutput', false));
set(gca, 'YTick', 1:nMeas, 'YTickLabel', measures);
xlabel('UMAP Dimension');
ylabel('Ripple Measure');
title('UMAP Dimensions vs Ripple Characteristics');
set(gca, 'FontSize', 11);

% Print interpretation
fprintf('\n=== UMAP Dimension Interpretation ===\n\n');
for dim = 1:nComps
    fprintf('UMAP%d:\n', dim);
    [~, ord] = sort(abs(meas_corr(dim, :)), 'descend');
    for i = 1:nMeas
        idx = ord(i);
        fprintf('  %s: r=%.3f\n', measures{idx}, meas_corr(dim, idx));
    end
    fprintf('\n');
end

% Scatter plots: each UMAP dimension vs its strongest correlate
figure('Position', [100 100 1200 300]);
for dim = 1:nComps
    [~, best_idx] = max(abs(meas_corr(dim, :)));
    best_meas = measures{best_idx};
    
    subplot(1, nComps, dim);
    hold on;
    for g = 1:3
        idx = all_ripples.expcond == groups{g};
        scatter(reduction(idx, dim), double(all_ripples.(best_meas)(idx)), 6, cols(g, :), 'filled', 'MarkerFaceAlpha', 0.6);
    end
    xlabel(sprintf('UMAP%d', dim));
    ylabel(best_meas);
    title(sprintf('r=%.3f', meas_corr(dim, best_idx)));
    set(gca, 'FontSize', 10);
    hold off;
end
sgtitle('Each UMAP Dimension vs Its Strongest Ripple Correlate');
%% 8. (PART 2) DENSITY CONTOURS


% 8.6. Density contours per group
fig2 = figure('Position', [100 100 1200 400]);
for g = 1:numel(groups)
    subplot(1, numel(groups), g);
    idx = groupLabels == groups(g);
    
    if sum(idx) < 10
        % Skip density plot for very small groups
        text(0.5, 0.5, sprintf('Insufficient data (n=%d)', sum(idx)), ...
            'HorizontalAlignment', 'center');
        title(sprintf('%s (n=%d)', char(groups(g)), sum(idx)));
        continue;
    end
    
    % Create 2D density estimate
    pts = reduction(idx, 1:2);
    [bandwidth, density, X, Y] = kde2d(pts);
    contour(X, Y, density, 10);
    colormap(hot); colorbar;
    title(sprintf('%s (n=%d)', char(groups(g)), sum(idx)));
    xlabel('UMAP 1'); ylabel('UMAP 2');
end
sgtitle('2D Density Contours per Group');

% 8.7. Convex hulls
% fig3 = figure('Position', [100 100 600 500]); hold on;
% for i = 1:numel(groups)
%     idx = groupLabels == groups(i);
%     pts = reduction(idx, 1:2);
%     k = convhull(pts(:,1), pts(:,2));
%     plot(pts(k,1), pts(k,2), 'Color', colors(i,:), 'LineWidth', 2);
%     scatter(pts(:,1), pts(:,2), 20, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.4);
% end
% xlabel('UMAP 1'); ylabel('UMAP 2');
% title('Convex Hulls per Group');
% legend(groups, 'Location', 'best');
% grid on; axis tight;


%% 9. Load previously saved CSD data and layers if not in workspace

all_csd = all_csd_layers;
filtered_allcsd = all_csd(:, :, ~rowidx);
all_csd = filtered_allcsd;

% 9.1. Parameters
Fs = 1250;  % Sampling frequency (Hz) - ADJUST TO YOUR DATA
window_ms = 25;  % Window duration (ms)
n_samples = round(window_ms * Fs / 1000);
n_pcs = 20;  % Number of PCs to retain
n_clusters = 5;  % Number of clusters for k-means
n_ripples_total = size(all_csd, 3);

group_labels = zeros(n_ripples_total, 1);
group_labels(strcmp(string(all_ripples.expcond), "SS")) = 1;
group_labels(strcmp(string(all_ripples.expcond), "OS")) = 2;
group_labels(strcmp(string(all_ripples.expcond), "OD")) = 3;

% 9.2. Extract 25ms centered windows
fprintf('Extracting 25ms centered windows...\n');
[n_time, n_channels, ~] = size(all_csd);
center_idx = ceil(n_time / 2);
half_win = floor(n_samples / 2);
idx_start = center_idx - half_win;
idx_end = center_idx + half_win;
csd_windows = all_csd(idx_start:idx_end, :, :);

% 9.3. Vectorize
fprintf('Vectorizing CSD matrices...\n');
X = reshape(permute(csd_windows, [3 1 2]), n_ripples_total, []);

% 9.4. Z-score features
fprintf('Z-scoring features...\n');
X_zscore = zscore(X);
[coeff, score, latent, ~, explained, mu] = pca(X_zscore);
score_reduced = score(:, 1:n_pcs);
% cum_var = cumsum(explained);
% 8. K-means clustering
% 9.5 UMAP on PC scores
fprintf('Running UMAP...\n')
reduction = run_umap(score_reduced, 'n_components', 2, 'n_neighbors', 20, 'min_dist', 0.1);
% reduction = run_umap(score_reduced, 'n_components', 2, 'n_neighbors', 15, 'min_dist', 0.1); % USING PCA SCORES

fprintf('Clustering with k-means...\n');
cluster_labels = kmeans(reduction, n_clusters, 'Replicates', 20);

% 9.6.Define groups for plotting
groups = {'SS', 'OS', 'OD'};
group_colors = [0.2 0.6 0.9; 0.9 0.4 0.2; 0.3 0.8 0.3];

%% 

function [bandwidth, density, X, Y] = kde2d(data, n, MIN_XY, MAX_XY)
    % Fast 2D kernel density estimation
    % Adapted for bivariate data
    if nargin < 2 || isempty(n), n = 2^8; end
    if nargin < 3 || isempty(MIN_XY)
        MIN_XY = [min(data(:,1)), min(data(:,2))];
        MAX_XY = [max(data(:,1)), max(data(:,2))];
    end
    
    scaling = MAX_XY - MIN_XY;
    
    % Bandwidth estimation (Silverman's rule)
    sigma = median(abs(data - median(data))) / 0.6745;
    bandwidth = sigma * (size(data,1)^(-1/6));
    
    % Create grid
    x = linspace(MIN_XY(1), MAX_XY(1), n);
    y = linspace(MIN_XY(2), MAX_XY(2), n);
    [X, Y] = meshgrid(x, y);
    
    % Compute density
    density = zeros(n, n);
    for i = 1:size(data, 1)
        density = density + exp(-((X - data(i,1)).^2 + (Y - data(i,2)).^2) / (2*bandwidth(1)^2));
    end
    density = density / (size(data,1) * 2 * pi * bandwidth(1)^2);
end

% --- Utility for extracting pairwise contrasts robustly from fitlme ---
% returns pMain and three pairwise Bonferroni-adjusted p-values [OS_vs_SS, OD_vs_SS, OD_vs_OS]
function [pMain, pPairsAdj] = lme_main_and_pairs(tbl, responseVar)
    tbl.Response = tbl.(responseVar);
    lme = fitlme(tbl, 'Response ~ expcond + (1|subject) + (1|ID)'); % nested random effects
    coeffNames = lme.CoefficientNames; % e.g. {'(Intercept)','expcond_OS','expcond_OD', ...}
    % main effect: test that all expcond coefficients are zero
    expIdx = find(startsWith(coeffNames,'expcond'));
    if ~isempty(expIdx)
        C = zeros(numel(expIdx), numel(coeffNames));
        for i = 1:numel(expIdx), C(i,expIdx(i)) = 1; end
        [~, pMain] = coefTest(lme, C);
    else
        % fallback: ANOVA row lookup (less robust, kept for completeness)
        atab = anova(lme,'DFMethod','satterthwaite');
        rn = atab.Properties.RowNames;
        idx = find(contains(rn,'expcond'),1);
        if ~isempty(idx), pMain = atab.pValue(idx); else pMain = NaN; end
    end
    % pairwise contrasts (SS is reference): find indices for expcond_OS and expcond_OD
    idxOS = find(strcmp(coeffNames,'expcond_OS'),1);
    idxOD = find(strcmp(coeffNames,'expcond_OD'),1);
    p_os_ss = NaN; p_od_ss = NaN; p_od_os = NaN;
    if ~isempty(idxOS)
        c = zeros(1,numel(coeffNames)); c(idxOS) = 1;
        [~, p_os_ss] = coefTest(lme, c);
    end
    if ~isempty(idxOD)
        c = zeros(1,numel(coeffNames)); c(idxOD) = 1;
        [~, p_od_ss] = coefTest(lme, c);
    end
    if ~isempty(idxOS) && ~isempty(idxOD)
        c = zeros(1,numel(coeffNames)); c(idxOD) = 1; c(idxOS) = -1;
        [~, p_od_os] = coefTest(lme, c);
    end
    raw = [p_os_ss, p_od_ss, p_od_os];
    pPairsAdj = min(raw * 3, 1); % Bonferroni for 3 comparisons
end