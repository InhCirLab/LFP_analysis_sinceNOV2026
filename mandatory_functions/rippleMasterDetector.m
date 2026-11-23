function [ripples,SW] = rippleMasterDetector(varargin)
%   rippleMasterDetector - Wrapped function to compute different
%                           characteristics about hippocampal ripples (100
%                           ~ 200 Hz oscillations). It also computes
%                           SharpWaves based on the ripples detected.
%
% USAGE
%   [ripples] = rippleMasterDetector(<options>)
%   
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.  The estimated EMG can be used as an additional
%    exclusion criteria.
%
%   SharpWaves are detected based on the detected ripples by findSharpWaves. Radiatum lfp
%   signal is filtered (default [2 10] Hz) and szcore of the signal is
%   computed. SW peak is detected as the time when zscore rad signal
%   exceeds SWthreshold(2) during the ocurrence of a ripple (SW.peaks). If a ripple
%   does not have an associated SharpWave, nan values are in play. Onset
%   and offset of the sharpwave is computed as the time when radiatum
%   zscore signal first crosses threshold(1) both before and after the
%   peaks ( SW.timestamps)
%
% INPUTS - note these are NOT name-value pairs... just raw values
%    <options>      optional list of property-value pairs (see tables below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5]); must be integer values
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [30 100]). 
%     'minDuration' min ripple duration. Keeping this input nomenclature for backwards
%                   compatibility
%     'restrict'    interval used to compute normalization. By default,
%                       tries to remove stimulated periods. Otherwise, 'off'
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy unfiltered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%     'passband'    N x 2 matrix of frequencies to filter for ripple detection 
%                   (default = [130 200])
%     'EMGThresh'   0-1 threshold of EMG to exclude noise
%     'saveMat'     logical (default=false) to save in buzcode format
%     'plotType'   1=original version (several plots); 2=only raw lfp
%     'detector'   'filter' (default), 'cnn', 'consensus'
%    =========================================================================
%
% OUTPUT
%
%    ripples        buzcode format .event. struct with the following fields
%                   .timestamps        Nx2 matrix of start/stop times for
%                                      each ripple
%                   .detectorName      string ID for detector function used
%                   .peaks             Nx1 matrix of peak power timestamps 
%                   .stdev             standard dev used as threshold
%                   .noise             candidate ripples that were
%                                      identified as noise and removed
%                   .peakNormedPower   Nx1 matrix of peak power values
%                   .detectorParams    struct with input parameters given
%                                      to the detector
%   SW              buzcode format .event. struct with the following fields
%                   .timestamps
%                   .detectorName
%                   .peaks
% SEE ALSO
%
%    See also bz_Filter, bz_RippleStats, bz_SaveRippleEvents, bz_PlotRippleStats.
%   
%   Develop by Manu Valero and Pablo Abad 2022. Buzsaki Lab. Based on
%   bz_findRipples
warning('this function is under development and may not work... yet')

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'SWChannel',[],@isnumeric);
addParameter(p,'thresholds',[2 5],@isnumeric);
addParameter(p,'SWthresholds',[-0.5 -2], @isnumeric);
addParameter(p,'durations',[30 120],@isnumeric);
addParameter(p,'restrict',[],@isnumeric);
addParameter(p,'frequency',1250,@isnumeric);
addParameter(p,'stdev',[],@isnumeric);
addParameter(p,'show','off',@isstr);
addParameter(p,'noise',[],@ismatrix);
addParameter(p,'passband',[80 200],@isnumeric);
addParameter(p,'SWpassband',[2 10],@isnumeric);
addParameter(p,'EMGThresh',1,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'minDuration',20,@isnumeric);
addParameter(p,'plotType',2,@isnumeric);
addParameter(p,'srLfp',1250,@isnumeric);
addParameter(p,'rippleStats',true,@islogical);
addParameter(p,'debug',false,@islogical);
addParameter(p,'eventSpikeThreshold',.5);
addParameter(p,'eventSpikeThreshold_shanks','all');
addParameter(p,'force',false,@islogical);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'useCSD',false,@islogical);
addParameter(p,'stdThreshold',2,@isnumeric);
addParameter(p,'detector','filter',@ischar);
addParameter(p,'excludeIntervals',[],@isnumeric);
% -- cnn related --
addParameter(p,'cnn_channels', [], @isnumeric);
addParameter(p,'model_file', '', @ischar);
addParameter(p,'pred_every', 0.032, @isnumeric);
addParameter(p,'verbose', false, @islogical);
addParameter(p,'handle_overlap', 'mean', @ischar);
addParameter(p,'exec_env', '', @ischar);
addParameter(p,'thetaepochs', '', @isstruct);

parse(p,varargin{:})

basepath = p.Results.basepath;
rippleChannel = p.Results.rippleChannel;
SWChannel = p.Results.SWChannel;
thresholds = p.Results.thresholds;
SWthresholds = p.Results.SWthresholds;
durations = p.Results.durations;
restrict = p.Results.restrict;
frequency = p.Results.frequency;
stdev = p.Results.stdev;
show = p.Results.show;
noise = p.Results.noise;
passband = p.Results.passband;
SWpassband = p.Results.SWpassband;
EMGThresh = p.Results.EMGThresh;
saveMat = p.Results.saveMat;
minDuration = p.Results.minDuration;
plotType = p.Results.plotType;
srLfp = p.Results.srLfp;
rippleStats = p.Results.rippleStats;
debug = p.Results.debug;
eventSpikeThreshold = p.Results.eventSpikeThreshold;
eventSpikeThreshold_shanks = p.Results.eventSpikeThreshold_shanks;
force = p.Results.force;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
useCSD = p.Results.useCSD;
stdThreshold = p.Results.stdThreshold;
detector = p.Results.detector;
excludeIntervals = p.Results.excludeIntervals;
% -- cnn related --
cnn_channels = p.Results.cnn_channels;
pred_every = p.Results.pred_every;
verbose = p.Results.verbose;
model_file = p.Results.model_file;
handle_overlap = p.Results.handle_overlap;
exec_env = p.Results.exec_env;
dir_cnn = fileparts(which('detect_ripples_cnn.m'));
thetaepochs = p.Results.thetaepochs;


%% Load Session Metadata and several variables if not provided
% session = sessionTemplate(basepath,'showGUI',false);
session = loadSession(basepath);

if (exist([session.general.name '.ripples.events.mat'],'file') ...
        && ~force)
    disp(['Ripples already detected for ', session.general.name, '. Loading file.']);
    load([session.general.name '.ripples.events.mat']);
    return
end

% Ripple and SW Channel are loaded separately in case we want to provide
% only one of the
if isempty(rippleChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers();
    end
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

if isempty(SWChannel)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    else
        [hippocampalLayers] = getHippocampalLayers();
    end
    SWChannel = hippocampalLayers.bestShankLayers.radiatum;
end

if skipStimulationPeriods && ~isempty(dir('*optogeneticPulses.events.mat'))
    targetFile = dir('*optogeneticPulses.events.mat'); load(targetFile.name);
    restrict_temp = SubtractIntervals([0 Inf],optoPulses.stimulationEpochs);
    restrict =  ConsolidateIntervals([restrict; restrict_temp; restrict_temp]);
end

if useCSD
    disp('Computing CSD...');
    rippleChannel = computeCSD([],'channels',rippleChannel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing Ripples
%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(detector, 'filter')
    ripples = findRipples( rippleChannel,'thresholds',thresholds,'passband',passband,...
        'EMGThresh',EMGThresh,'durations',durations, 'saveMat',false,'restrict',restrict,'frequency',frequency,'excludeIntervals',excludeIntervals,'noise',noise);
elseif strcmp(detector, 'filtercnn')
    ripples = findRipples_wCNN(rippleChannel,'thresholds',thresholds,'passband',passband,...
        'EMGThresh',EMGThresh,'durations',durations, 'saveMat',false,'restrict',restrict,...
        'frequency',frequency,'excludeIntervals',excludeIntervals,'noise',noise,...
        'useCNN',true,'cnnChannels', rippleChannel-2:rippleChannel+5,...
        'cnnExecEnv', 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe',...
        'cnnModelFile', 'C:\Users\usuario\cnn-matlab-master\cnn',...
        'cnnThreshold', 0.7, 'cnnwinSize', 0.1);
elseif strcmp(detector, 'customcnn')
    % With basepath
    ripples = findRipples_CNN_modpeak(rippleChannel, ...
        'cnnChannels', rippleChannel-2:rippleChannel+5, ...
        'cnnExecEnv', 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe', ...
        'cnnModelFile', 'C:\Users\usuario\cnn-matlab-master\cnn', ...
        'restrict', excludeIntervals, ...
        'passband', passband,...
        'saveMat', true);
    % ripples = findRipples_CNNv1(rippleChannel, ...
    %     'cnnChannels', rippleChannel-2:rippleChannel+5, ...
    %     'cnnExecEnv', 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe', ...
    %     'cnnModelFile', 'C:\Users\usuario\cnn-matlab-master\cnn', ...
    %     'restrict', excludeIntervals, ...
    %     'saveMat', true);
elseif strcmp(detector, 'cnn')
    % Load channel configuration
    load([session.general.name,'.session.mat']);
    channelConf = session.extracellular.electrodeGroups.channels;
    % Get best shank
    % if isvarname('hippocampalLayers')
    if exist('hippocampalLayers')
        bestShank = hippocampalLayers.bestShank;
    elseif ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
        bestShank = hippocampalLayers.bestShank;
    else
        [hippocampalLayers] = getHippocampalLayers();
        bestShank = hippocampalLayers.bestShank;
    end
    % Get channels near pyr that fits in the shank
    channelsShank = channelConf{bestShank};
    rippleChPos = find(channelsShank==rippleChannel);
    if (rippleChPos < length(channelsShank)-5) && (rippleChPos > 2)
        cnn_channels = channelsShank(rippleChPos + (-2:5));
    elseif (rippleChPos < length(channelsShank)-5)
        cnn_channels = channelsShank(1:8);
    elseif (rippleChPos > 2)
        cnn_channels = channelsShank(end-7:end);
    end
    % Get LFP
    lfp = getLFP(cnn_channels,'basepath',basepath);
    % Get model directory
    if isempty(model_file)
        model_file = fullfile(dir_cnn, 'cnn');
    end
    % Predict probability
    SWRprob = detect_ripples_cnn(double(lfp.data), lfp.samplingRate, 'model_file', model_file, ...
        'pred_every', pred_every, 'verbose', verbose, 'handle_overlap', handle_overlap, ...
        'exec_env', exec_env);
    ripples.timestamps = get_intervals(SWRprob,'LFP', double(lfp.data), 'fs', lfp.samplingRate, 'win_size', 0.1, 'threshold', 0.7);
    % ripples.duration = ripples.timestamps(:, 2) - ripples.timestamps(:, 1);
    % ripples.peaks = sum(ripples.timestamps, 2)/2;
    % ripples.detectorinfo.detectorname = 'cnn';
    % % ripples.detectorinfo.detectionparams.poneraquilosparamsdelcnn = params;
    % ripples.detectorinfo.detectionchannel = rippleChPos;
    % 
    % % SOME SIMILAR PROCEDURES AS FOR CLASSICAL FILTER (FINDRIPPLES)
    % 
    % % 1. Get normalized squared signal
    % %% filter and calculate noise
    % 
    %     % Parameters
    %     frequency = lfp.samplingRate;
    %     signal = bz_Filter(lfp,'filter','butter','passband',passband,'order',3);
    %     timestamps = signal.timestamps;
    %     signal = signal.data;
    %     windowLength = frequency/frequency*11;
    % 
    %     % Square and normalize signal
    %     squaredSignal = signal.^2;
    %     % squaredSignal = abs(opsignal);
    %     window = ones(windowLength,1)/windowLength;
    %     keep = [];
    %     if ~isempty(restrict)
    %         for i=1:size(restrict,1)
    %             keep = InIntervals(timestamps,restrict);
    %         end
    %     end
    %     keep = logical(keep); 
    % 
    %     [normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),stdThreshold,keep);
    % 
    %     peakNormalizedPower = [];
    %     thirdPass = [];
    %     for i = 1:size(ripples.timestamps,1)
    %         onsetsample = find(timestamps >= ripples.timestamps(i, 1), 1, 'first');
    %         offsetsample = find(timestamps <= ripples.timestamps(i, 2), 1, 'last');
	%         [maxValue,maxIndex] = max(normalizedSquaredSignal(onsetsample:offsetsample));
	%         if maxValue > thresholds(2)
	% 	        thirdPass = [thirdPass ; ripples.timestamps(i,:)];
	% 	        peakNormalizedPower = [peakNormalizedPower ; maxValue];
	%         end
    %     end

elseif strcmp(detector, 'consensus')
    
else
    warning([detector ' is not recognised as a detector, please chose between filter/cnn/consensus'])
    return
end
%%
if skipStimulationPeriods
    try
        % Remove ripples durting stimulation artifacts
        if ~isempty(dir('*.optogeneticPulses.events.mat'))
            f = dir('*.optogeneticPulses.events.mat');
            disp('Using stimulation periods from optogeneticPulses.events.mat file');
            load(f.name);
            pulPeriods = optoPulses.stimulationEpochs;

        elseif ~isempty(dir('.pulses.events.mat'))
            f = dir('*Pulses.events.mat');
            disp('Using stimulation periods from pulses.events.mat file');
            load(f.name);
            pulPeriods = pulses.intsPeriods;
        else
            warning('No pulses epochs detected!');
        end
        nBefore = length(ripples.peaks);
        for i = 1:size(pulPeriods,1)
            a = InIntervals(ripples.peaks,pulPeriods(i,:));
            fieldsR = fields(ripples);
            for j = 1:size(fieldsR,1)
                if ~isstruct(ripples.(fieldsR{j})) && size(ripples.(fieldsR{j}),1) > 3
                    ripples.(fieldsR{j})(a,:) = [];
                end
            end
        end
        fprintf('%3.i/%3.i events during stimulation, discarted... \n',nBefore-length(ripples.peaks), nBefore); %\n
    catch
        warning('Not possible to remove ripples during stimulation epochs...');
    end
end

ripples = removeArtifactsFromEvents(ripples,'stdThreshold',stdThreshold);

if isnumeric(eventSpikeThreshold) || eventSpikeThreshold
    if islogical(eventSpikeThreshold)
        eventSpikeThreshold = 1;
    end
    ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',eventSpikeThreshold,'shanksID',eventSpikeThreshold_shanks);
end

% plotRippleChannel('rippleChannel',rippleChannel,'ripples',ripples); % to do, run this after ripple detection
% EventExplorer(pwd, ripples)

% Apply same filter as when computing CSD (if too close to the beginning or
% end of the recording, remove ripple from analysis)

events = ripples.peaks;
% Convert event times (s) to samples
samplingRate = frequency;
events_samples = round(events * samplingRate);

% Convert window [before after] (s) to samples
twin_samples = round([0.1 0.1] * samplingRate);

% Keep only events whose full time window fits inside the LFP data
lfp = getLFP(rippleChannel, 'basepath', basepath);
valid_events = (events_samples + twin_samples(2) <= size(lfp.data,1)) & ...
               (events_samples - twin_samples(1) > 0);

ripples.timestamps = ripples.timestamps(valid_events, :);
ripples.peakNormedPower = ripples.peakNormedPower(valid_events, :);
ripples.peaks = ripples.peaks(valid_events, :);

%% Ripple Stats
if rippleStats
    files = dir('*.thetaEpochs6to12.states.mat');
    fullpath = fullfile(files(1).folder, files(1).name);
    load(fullpath);
    ripples = computeRippleStats_new('ripples',ripples,'rippleChannel',rippleChannel,'bandpass',passband, 'thetaEpochs', thisthetaEpochs);
    % ripples = computeRippleStats('ripples',ripples,'rippleChannel',rippleChannel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing SharpWaves
%%%%%%%%%%%%%%%%%%%%%%%%%
try
    SW = findSharpWaves('ripples',ripples,'rippleChannel',rippleChannel,'SWChannel',SWChannel,...
        'passband',passband,'SWpassband',SWpassband);
catch
    disp('Not possible to compute Sharp Waves...');
end
%% OUTPUT
if saveMat
    disp('Saving Ripples Results...');
    save([session.general.name , '.ripples.events.mat'],'ripples');
    try
        disp('Saving SharpWaves Results...');
        save([session.general.name , '.sharpwaves.events.mat'],'SW');
    catch
    end
end


end

% function y = Filter0(b,x)
% 
% if size(x,1) == 1
% 	x = x(:);
% end
% 
% if mod(length(b),2)~=1
% 	error('filter order should be odd');
% end
% 
% shift = (length(b)-1)/2;
% 
% [y0 z] = filter(b,1,x);
% 
% y = [y0(shift+1:end,:) ; z(1:shift,:)];
% 
% function [U,stdA] = unity(A,sd,restrict)
% 
% if ~isempty(restrict),
% 	meanA = mean(A(restrict));
% 	stdA = std(A(restrict));
% else
% 	meanA = mean(A);
% 	stdA = std(A);
% end
% if ~isempty(sd),
% 	stdA = sd;
% end
% 
% U = (A - meanA)/stdA;

