function [ripples] = findRipples_CNN_modpeak(varargin)
%findRipples_CNN - Find hippocampal ripples using CNN detection
%
% USAGE
%    [ripples] = findRipples_CNN(channel, <options>)
%    [ripples] = findRipples_CNN(basepath, channel, <options>)
%
% INPUTS
%    channel          Ripple channel to use (1-indexed)
%    basepath         Path to session (default = pwd)
%
%    =========================================================================
%     Properties      Values
%    -------------------------------------------------------------------------
%     'cnnChannels'   Channels for CNN detection (required)
%     'cnnExecEnv'    Path to Python executable (required)
%     'cnnModelFile'  Path to CNN model (required)
%     'cnnThreshold'  CNN probability threshold (default = 0.7)
%     'cnnWinSize'    Window size for intervals (default = 0.1)
%     'cnnPredEvery'  Prediction interval (default = 0.01)
%     'restrict'      Intervals to EXCLUDE from detection and normalization
%     'passband'      Bandpass for normalization [Hz] (default = [130 200])
%     'frequency'     Sampling rate in Hz (default = 1250)
%     'saveMat'       Save results (default = false)
%     'basepath'      Session path (default = pwd)
%    =========================================================================
%
% OUTPUT
%    ripples.timestamps        [Nx2] start/stop times (seconds)
%    ripples.timestamps_samples [Nx2] start/stop samples
%    ripples.peaks             [Nx1] peak times (seconds)
%    ripples.peaks_samples     [Nx1] peak samples
%    ripples.peakNormedPower   [Nx1] normalized power at peaks
%    ripples.detectorinfo      Detection metadata

%% Parse inputs
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'cnnChannels',[],@isnumeric);
addParameter(p,'cnnExecEnv','',@ischar);
addParameter(p,'cnnModelFile','',@ischar);
addParameter(p,'cnnThreshold',0.7,@isnumeric);
addParameter(p,'cnnWinSize',0.1,@isnumeric);
addParameter(p,'cnnPredEvery',0.01,@isnumeric);
addParameter(p,'restrict',[],@isnumeric);
addParameter(p,'passband',[80 250],@isnumeric);
addParameter(p,'frequency',1250,@isnumeric);
addParameter(p,'saveMat',false,@islogical);

if ischar(varargin{1}) || isstring(varargin{1})
    addRequired(p,'basepath',@ischar);
    addRequired(p,'channel',@isnumeric);
    parse(p,varargin{:});
    basepath = p.Results.basepath;
    channel = p.Results.channel;
else
    addRequired(p,'channel',@isnumeric);
    parse(p,varargin{:});
    basepath = p.Results.basepath;
    channel = p.Results.channel;
end

% Extract parameters
cnnChannels = p.Results.cnnChannels;
cnnExecEnv = p.Results.cnnExecEnv;
cnnModelFile = p.Results.cnnModelFile;
cnnThreshold = p.Results.cnnThreshold;
cnnWinSize = p.Results.cnnWinSize;
cnnPredEvery = p.Results.cnnPredEvery;
restrict = p.Results.restrict;
passband = p.Results.passband;
frequency = p.Results.frequency;
saveMat = p.Results.saveMat;

% Get basename
prevPath = pwd;
cd(basepath);
basename = basenameFromBasepath(basepath);

%% CNN Detection
disp('Running CNN detection...');
LFP_all = getLFP('all','basepath',basepath);
fs = LFP_all.samplingRate;
% Para datos filtrados:
%data = bz_Filter(double(LFP_all.data),'filter','butter','passband',passband,'order',3);
% input de detector en crudo:
data = double(LFP_all.data);

% disp('Resampling LFP to 1250 Hz for CNN compatibility...');
% nChannels = size(LFP_all.data,2);
% nSamplesResampled = round(size(LFP_all.data,1) * 1250 / 20000);
% data = zeros(nSamplesResampled, nChannels);
% for ch = 1:nChannels
%     data(:,ch) = resample(double(LFP_all.data(:,ch)), 1250, 20000)';
% end

% Para datos z-scoreados:
% data = zscore(double(LFP_all.data)); % Normalizado dentro de cada canal; por columnas

y_pred = detect_ripples_cnn(data, fs, 'channels', cnnChannels, ...
    'exec_env', cnnExecEnv, 'model_file', cnnModelFile, 'pred_every', cnnPredEvery);

intervals = get_intervals(y_pred, 'LFP', data, 'fs', fs, ...
    'win_size', cnnWinSize, 'threshold', cnnThreshold);

if isempty(intervals)
    disp('No ripples detected by CNN');
    ripples = struct('timestamps', [], 'timestamps_samples', [], 'peaks', [], ...
        'peaks_samples', [], 'peakNormedPower', [], 'detectorinfo', []);
    cd(prevPath);
    return
end

% Exclude intervals in restrict periods
if ~isempty(restrict)
    keep = true(size(intervals,1),1);
    for i = 1:size(intervals,1)
        peakTime = mean(intervals(i,:));
        if any(peakTime >= restrict(:,1) & peakTime <= restrict(:,2))
            keep(i) = false;
        end
    end
    intervals = intervals(keep,:);
    disp(['After restriction exclusion: ' num2str(size(intervals,1)) ' events.']);
end

if isempty(intervals)
    disp('No ripples remain after restriction');
    ripples = struct('timestamps', [], 'timestamps_samples', [], 'peaks', [], ...
        'peaks_samples', [], 'peakNormedPower', [], 'detectorinfo', []);
    cd(prevPath);
    return
end

%% Calculate Normalized Power and Extract Peaks
disp('Computing normalized power...');
lfp = getLFP(channel,'basepath',basepath,'basename',basename);
signal = bz_Filter(double(lfp.data),'filter','butter','passband',passband,'order',3);
timestamps = lfp.timestamps;

% Square signal
squaredSignal = signal.^2;

% Smooth with moving window
windowLength = round(frequency/frequency*11);
window = ones(windowLength,1)/windowLength;
smoothedSignal = Filter0(window, squaredSignal);

% Compute normalization excluding restrict periods
keep = true(size(timestamps));
if ~isempty(restrict)
    keep = ~InIntervals(timestamps, restrict);
end

% Z-score
meanSig = mean(smoothedSignal(keep));
stdSig = std(smoothedSignal(keep));
normalizedSquaredSignal = (smoothedSignal - meanSig) / stdSig;

%% Extract peaks and normalized power
nRipples = size(intervals,1);
peaks = zeros(nRipples,1);
peaks_samples = zeros(nRipples,1);
peakNormedPower = zeros(nRipples,1);
timestamps_samples = zeros(nRipples,2);

for i = 1:nRipples
    % Find samples corresponding to interval
    startSample = find(timestamps >= intervals(i,1), 1, 'first');
    stopSample = find(timestamps <= intervals(i,2), 1, 'last');
    
    if isempty(startSample) || isempty(stopSample)
        continue;
    end
    
    % Store sample indices
    timestamps_samples(i,:) = [startSample, stopSample];
    
    % Find negative peak (trough) of filtered signal - like original findRipples
    [~, minIndex] = min(signal(startSample:stopSample));
    peakSample = startSample + minIndex - 1;
    
    peaks_samples(i) = peakSample;
    peaks(i) = timestamps(peakSample);
    
    % Get normalized power at that peak
    peakNormedPower(i) = normalizedSquaredSignal(peakSample);
end

%% Build output structure
ripples.timestamps = intervals;
ripples.timestamps_samples = timestamps_samples;
ripples.peaks = peaks;
ripples.peaks_samples = peaks_samples;
ripples.peakNormedPower = peakNormedPower;
ripples.stdev = stdSig;

% Detector info
detectorinfo.detectorname = 'findRipples_CNN';
detectorinfo.detectiondate = date;
detectorinfo.detectionintervals = restrict;
detectorinfo.detectionchannel = channel;
detectorinfo.detectionparms.cnnChannels = cnnChannels;
detectorinfo.detectionparms.cnnThreshold = cnnThreshold;
detectorinfo.detectionparms.cnnWinSize = cnnWinSize;
detectorinfo.detectionparms.cnnPredEvery = cnnPredEvery;
detectorinfo.detectionparms.passband = passband;
detectorinfo.detectionparms.frequency = frequency;

ripples.detectorinfo = detectorinfo;

disp(['CNN detection complete: ' num2str(nRipples) ' ripples detected.']);

%% Save
if saveMat
    save(fullfile(basepath, [basename '.ripples.events.mat']), 'ripples');
end

cd(prevPath);

end

%% Helper functions
function y = Filter0(b,x)
if size(x,1) == 1
    x = x(:);
end
if mod(length(b),2)~=1
    error('filter order should be odd');
end
shift = (length(b)-1)/2;
[y0, z] = filter(b,1,x);
y = [y0(shift+1:end,:) ; z(1:shift,:)];
end