function [ripples, maps] = computeRippleStats_new(varargin)

%RippleStats - Compute descriptive stats for ripples (100~200Hz oscillations).
%
%  USAGE
%
%    [maps,data,stats] = bz_RippleStats(filtered,timestamps,ripples,<options>)
%
%    filtered       ripple-band filtered samples (one channel)
%    timestamps	    timestamps (in seconds) to match filtered vector
%    ripples        ripple timing information (STRUCT VERSION) (obtained using <a href="matlab:help bz_FindRipples">bz_FindRipples</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bandpass'    Default, [80 200]
%     'ripples'     By default tries to load ripple.event.mat structure or
%                   runs findRipples. Otherwise, ripple timing information structure(see
%                       findRipples)
%     'durations'   durations before and after ripple peak (in s)
%                   (default = [-0.075 0.075])
%     'basepath'    Default, pwd
%     'saveSummary' Default, true
%     'saveMat'     Save or updates ripples.event.mat structure with status info.
%                       Default, true
%     'plotOpt'     Default, true
%     'rippleChannel'
%                   By default loads getHippocampalLayers, and uses
%                   hippocampalLayer.bestShankLayers.pyramidal (1-index)
%    =========================================================================
%
%  OUTPUT
%
%    maps.ripples               instantaneous amplitude (one ripple per row)
%    maps.frequency             instantaneous frequency
%    maps.phase                 instantaneous phase
%    maps.amplitude             enveloppe amplitude
%    data.peakFrequency         frequency at peak
%    data.peakAmplitude         amplitude at peak
%    data.duration              durations
%    stats.durationAmplitude    correlation between duration and amplitude (rho, p)
%    stats.durationFrequency    correlation between duration and frequency (rho, p)
%    stats.amplitudeFrequency   correlation between amplitude and frequency (rho, p)
%    stats.acg.data             autocorrelogram data
%    stats.acg.t                autocorrelogram time bins
%    and more...

% edited by David Tingley to fit buzcode formatting standards, 2017
% edited by Manuel Valero 2022 (fix inputs, include FRindex and entropy, as in Foffani et al, 2007; Valero et al, 2017)
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'bandpass',[120 200],@isnumeric);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'durations',[-0.075 0.075], @isnumeric);
addParameter(p,'ripples',[], @isstruct);
addParameter(p,'rippleChannel',[], @isnumeric);
addParameter(p,'plotOpt',[], @islogical);
addParameter(p, 'thetaEpochs', [], @isstruct);

parse(p,varargin{:})
basepath = p.Results.basepath;
bandpass = p.Results.bandpass;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
durations = p.Results.durations;
ripples = p.Results.ripples;
rippleChannel = p.Results.rippleChannel;
plotOpt = p.Results.plotOpt;
thetaEpochs = p.Results.thetaEpochs;

% Dealing with inputs
prevBasepath = pwd;
cd(basepath);

% session = sessionTemplate(basepath,'showGUI',false);
if isempty(rippleChannel)
    hippocampalLayers = getHippocampalLayers;
    rippleChannel = hippocampalLayers.bestShankLayers.pyramidal;
end

if isempty(ripples)
    ripples = findRipples(rippleChannel,'thresholds',[1.5 3.5],'passband',bandpass,...
        'EMGThresh',1,'durations',[20 250], 'saveMat', saveMat);
end

% some fixed parameters
corrBinSize = 0.01;

% % get non-loc lfp (including ripples)
% lfp = getLFP(rippleChannel);
% ts = lfp.timestamps;
% % Load theta epochs
% % 
% % thetaEpochs = thisthetaEpochs;
% 
% % Get non-loc lfp
% nonloclfp = lfp;
% nonlocint = thetaEpochs.intervals_nonTheta;
% nonlocidx = events_to_logical(ts, nonlocint);
% nonloclfp.data = lfp.data(nonlocidx);
% nonloclfp.timestamps = ts(nonlocidx);
% lfp = nonloclfp;

    lfp = getLFP(rippleChannel); 
    if ~isempty(thetaEpochs.thetaRun.newintervals)
        [a,b] = InIntervals(lfp.timestamps,thetaEpochs.thetaRun.newintervals);
        lfp.timestamps(a) = [];
        lfp.data(a,:) = [];
    end
    


% get filtered signal
filteredLFP = bz_Filter(lfp,'channels',rippleChannel,'filter','butter','passband',bandpass,'order',3);
samplingRate = filteredLFP.samplingRate;
timestamps = filteredLFP.timestamps;
filtered = filteredLFP.data;
unfiltered = lfp.data;

nBins = floor(samplingRate*diff(durations)/2)*2+1; % must be odd
nHalfCenterBins = 3;
centerBin = ceil(nBins/2);
centerBins = centerBin-nHalfCenterBins:centerBin+nHalfCenterBins;

%% SNIPPET FROM CHATGPT TO COMPUTE INST FREQUENCY ASSUMING RIPPLES ARE NON-STATIONARY

filteredLFP = bz_Filter(lfp,'channels',rippleChannel,'filter','butter','passband',bandpass,'order',3);
filtdata  = filteredLFP.data;

nRip = length(ripples.peaks);
ridgewin = [-0.05 0.05];

for k = 1:nRip

    % --- extract ripple waveform ---
    idx = timestamps >= ripples.peaks(k)+ridgewin(1) & timestamps <= ripples.peaks(k)+ridgewin(2);
    sig = filtdata(idx);

    % ── FIX: sig vacio o demasiado corto para cwt (ripple en borde del ──
    % recording o adyacente a intervalo theta eliminado)
    if numel(sig) < 10
        rippleIF.meanFreq(k,1) = NaN;
        rippleIF.stdFreq(k,1)  = NaN;
        rippleIF.drift(k,1)    = NaN;
        rippleIF.ridge{k}      = [];
        continue;
    end

    % --- wavelet transform ---
    try
        [wt,f] = cwt(sig, samplingRate);
    catch
        rippleIF.meanFreq(k,1) = NaN;
        rippleIF.stdFreq(k,1)  = NaN;
        rippleIF.drift(k,1)    = NaN;
        rippleIF.ridge{k}      = [];
        continue;
    end

    % keep ripple band only
    keep = f>=bandpass(1) & f<=bandpass(2);

    % ── FIX: keep vacio — ninguna frecuencia del CWT cae en la banda ────
    if ~any(keep)
        rippleIF.meanFreq(k,1) = NaN;
        rippleIF.stdFreq(k,1)  = NaN;
        rippleIF.drift(k,1)    = NaN;
        rippleIF.ridge{k}      = [];
        continue;
    end

    wt = abs(wt(keep,:));
    f  = f(keep);

    % --- ridge extraction ---
    [~,I] = max(wt,[],1);
    ridge = f(I);

    % ── FIX: ridge vacio por alguna razon residual ───────────────────────
    if isempty(ridge)
        rippleIF.meanFreq(k,1) = NaN;
        rippleIF.stdFreq(k,1)  = NaN;
        rippleIF.drift(k,1)    = NaN;
        rippleIF.ridge{k}      = [];
        continue;
    end

    % --- metrics ---
    rippleIF.meanFreq(k,1) = mean(ridge);
    rippleIF.stdFreq(k,1)  = std(ridge);
    rippleIF.drift(k,1)    = ridge(end)-ridge(1);
    rippleIF.ridge{k}      = ridge;
end

%%

% Get instantaneous phase and amplitude
phase = filteredLFP.phase;
amplitude = zscore(abs(filteredLFP.hilb));
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = bz_Diff(medfilt1(unwrapped,12),timestamps,'smooth',0);
frequency = frequency/(2*pi);

%% MULTITAPER/SPECTROGRAM BASED INTRINSIC FQ COMPUTATION (STATIONARY) 
lfpdata  = double(lfp.data);

% 1-s ripple-centered windows (required for multitaper)
[r,i] = Sync([timestamps lfpdata], ripples.peaks, 'durations', [-0.5 0.5]);
ripMat = SyncMap(r,i,'nbins',1250,'smooth',0);

% Multitaper spectra
spec = ripSpectrogram(ripMat', samplingRate, 0);

freqs = spec.Freq_range;
power = spec.Suma;         % freq × ripple

% Ripple band
bandIdx = freqs >= bandpass(1) & freqs <= bandpass(2);

% Intrinsic frequency = dominant spectral peak
[~,idx] = max(power(bandIdx,:),[],1);
f_band = freqs(bandIdx);

f_intrinsic = f_band(idx);

%% 
% idx = ceil((ripples.peaks-ripples.timestamps(:,1))*ripples.detectorinfo.detectionparms.frequency);

% Compute ripple map
[r,i] = Sync([timestamps filtered],ripples.peaks,'durations',durations);
maps.ripples_filtered = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);
[r,i] = Sync([timestamps double(unfiltered)],ripples.peaks,'durations',durations);
maps.ripples_raw = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute frequency Map
[f,i] = Sync([timestamps frequency],ripples.peaks,'durations',durations);
maps.frequency = SyncMap(f,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute phase map
[p,i] = Sync([timestamps phase],ripples.peaks,'durations',durations);
maps.phase = SyncMap(p,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute amplitude map
[a,i] = Sync([timestamps amplitude],ripples.peaks,'durations',durations);
maps.amplitude = SyncMap(a,i,'durations',durations,'nbins',nBins,'smooth',0);

maps.timestamps = linspace(durations(1),durations(2),nBins);

data.peakFrequency = maps.frequency(:,centerBin);
data.peakAmplitude = maps.amplitude(:,centerBin);

% Ripple durations
data.duration = abs(diff(ripples.timestamps'))';

% InterRipple Frequency
data.interRippleFrequency = [NaN; 1./diff(ripples.peaks)];

% stats
[stats.acg.data,stats.acg.t] = CCG(ripples.peaks,ones(length(ripples.peaks),1),'binSize',corrBinSize);

[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);

% fastRipple index and entropy (as computed in Foffani et al, 2007; Valero et al, 2017)
[r,i] = Sync([timestamps double(unfiltered)],ripples.peaks,'durations',[-0.5 0.5]);
matDouble = SyncMap(r,i,'durations',[-0.5 0.5],'smooth',0,'nbins',1250);

disp('Computing ripple features... ');
ripSpec = ripSpectrogram(double(matDouble'), samplingRate,0);

% collect results
data.spectralEntropy = ripSpec.entropyDada;
data.fastRippleIndex = ripSpec.frippindex;
data.multiTapperFreq = ripSpec.freqData;
data.rippleIF = rippleIF;
data.rippleIFstat = f_intrinsic;

maps.multitaperSpecs = ripSpec;

rippleStats.stats = stats;
rippleStats.data = data;
rippleStats.maps = maps;
ripples.rippleStats = rippleStats;

if plotOpt
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,['SummaryFigures\ripplesDetection.png']);
end
    
    h = figure;
    subplot(3,2,1)
    computeWavelet(mean(maps.ripples_raw), maps.timestamps, [80 300]);
    title(['Ch: ' num2str(ripples.detectorinfo.detectionchannel)],'FontWeight','normal');

    subplot(3,2,2)
    scatter(data.spectralEntropy, data.fastRippleIndex,10,[.7 .7 .7],"filled");
    xlabel('Spectral Entropy'); ylabel('Fast ripple index');

    subplot(3,2,3)
    [~, idx] = sort(data.fastRippleIndex);
    imagesc(maps.multitaperSpecs.Freq_range,[1 length(idx)],maps.multitaperSpecs.Suma(:,idx)');
    xlabel('Hz'); ylabel('#');

    subplot(3,2,4)
    imagesc(maps.timestamps,[1 length(idx)],maps.ripples_filtered(idx,:));
    xlabel('s');

    subplot(3,2,5)
    hold on
    plot(maps.timestamps, maps.ripples_raw/1E3,'color', [.8 .8 .8 .1]);
    plot(maps.timestamps, median(maps.ripples_raw)/1E3,'color', [0 0 0],'LineWidth',1.5);
    axis tight
    xlabel('s'); ylabel('Amp');

    subplot(3,2,6)
    hold on
    scatter(data.peakAmplitude, data.peakFrequency,10,data.duration,"filled");
    ylabel('Peak amplitude'); xlabel('Peak Freq');
    cb = colorbar;
    cb.Label.String = 's';
    
    if saveSummary
        mkdir('SummaryFigures'); % create folder
        saveas(gcf,'SummaryFigures\ripplesDetection.png');
    end

    if ~plotOpt
        close(h);
    end
    
    if saveMat
        disp('Saving Ripples Results...');
        save([basenameFromBasepath(pwd) , '.ripples.events.mat'],'ripples');
    end

    cd(prevBasepath);
end

