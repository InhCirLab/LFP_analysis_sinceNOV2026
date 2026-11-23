function [csd, lfpAvg, csd_events, lfp_events] = buzcode_eacheventCSD(lfp, events, varargin)

% [ csd, lfpAvg, csd_events, lfp_events ] = buzcode_eventCSD(lfp, events, varargin)
% Calculates event-triggered CSD map from linear array LFPs
% Now includes individual event data for later classification

% INPUT
%    lfp            buzcode structure with fields: data, timestamps, samplingRate
%   events          event timestamps (in sec)
%    <options>      
%       channels    vector with channels. Default: all
%       twin        time window [before after] in sec. Default [0.1 0.1]
%       spat_sm     spatial smoothing degree. Default = 11
%       temp_sm     temporal smoothing degree. Default = 11
%       doDetrend   logical. Default false
%       plotCSD     logical. Default true
%       plotLFP     logical. Default true
%       samplingRate  if lfp is numeric. Default 1250

% OUTPUT:
%    csd            average CSD structure
%    lfpAvg         average LFP structure  
%    csd_events     [time x channels x events] CSD for each event
%    lfp_events     [time x channels x events] LFP for each event

% Antonio FR, 7/18; Manu Valero, 2022; Modified for per-event storage, 2024

%% Parse inputs
p = inputParser;
addParameter(p,'channels',1:size(lfp.data,2),@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',11,@isnumeric);
addParameter(p,'temp_sm',11,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
parse(p,varargin{:});

channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;

% Handle lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp)
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = (1:length(lfp))' / samplingRate;
end

twin_samples = round(p.Results.twin * samplingRate);
events_samples = round(events * samplingRate);

%% Extract event-triggered LFP segments
valid_events = (events_samples + twin_samples(2) <= size(data,1)) & (events_samples - twin_samples(1) > 0);
events_samples = events_samples(valid_events);
n_events = length(events_samples);
n_samples = twin_samples(1) + twin_samples(2) + 1;
n_channels = length(channels);

% Preallocate
lfp_temp = nan(n_samples, n_channels, n_events);

% Extract all events
for e = 1:n_events
    lfp_temp(:,:,e) = data(events_samples(e)-twin_samples(1):events_samples(e)+twin_samples(2), channels) * -1;
end

%% Process each event individually
lfp_events = nan(n_samples, n_channels, n_events);
csd_events = nan(n_samples, n_channels-2, n_events);

for e = 1:n_events
    lfp_evt = lfp_temp(:,:,e);
    
    % Detrend
    if doDetrend
        lfp_evt = detrend(lfp_evt')';
    end
    
    % Temporal smoothing
    if temp_sm > 0
        for ch = 1:n_channels
            lfp_evt(:,ch) = smooth(lfp_evt(:,ch), temp_sm, 'sgolay');
        end
    end
    
    % Spatial smoothing
    if spat_sm > 0
        for t = 1:n_samples
            lfp_evt(t,:) = smooth(lfp_evt(t,:), spat_sm, 'lowess');
        end
    end
    
    lfp_events(:,:,e) = lfp_evt;
    csd_events(:,:,e) = diff(lfp_evt, 2, 2);
end

%% Compute averages
lfp_avg = nanmean(lfp_events, 3);
CSD = nanmean(csd_events, 3);

%% Generate output structures
taxis = -twin_samples(1):twin_samples(2);

csd.data = CSD;
csd.timestamps = taxis;
csd.samplingRate = samplingRate;
csd.channels = channels;
csd.params.spat_sm = spat_sm;
csd.params.temp_sm = temp_sm;
csd.params.detrend = doDetrend;

lfpAvg.data = lfp_avg;
lfpAvg.timestamps = taxis;
lfpAvg.samplingRate = samplingRate;
lfpAvg.channels = channels;
lfpAvg.params.spat_sm = spat_sm;
lfpAvg.params.temp_sm = temp_sm;
lfpAvg.params.detrend = doDetrend;

%% Plot
if plotLFP || plotCSD
    taxis_ms = (taxis / samplingRate) * 1e3;
    cmax = max(abs(CSD(:)));
    
    figure;
    
    if plotLFP
        subplot(1,2,1);
    end
    
    % CSD heatmap (channel 1 on top, higher channels at bottom)
    contourf(taxis_ms, 1:size(CSD,2), CSD', 40, 'LineColor', 'none');
    colormap jet; caxis([-cmax cmax]);
    set(gca, 'YDir', 'reverse');
    xlabel('time (ms)'); ylabel('channel'); title('CSD');
    hold on; plot([0 0], [1 size(CSD,2)], '--k');
    
    if plotLFP
        subplot(1,2,2);
        offset = 0;
        for ch = 1:size(lfp_avg,2)
            plot(taxis_ms, 1e0*lfp_avg(:,ch) + offset, 'k', 'LineWidth', 1.5); hold on;
            offset = offset + 400;
        end
        set(gca, 'YDir', 'reverse', 'YTickLabel', []);
        ylim([-1000 offset+1000]); xlim([taxis_ms(1) taxis_ms(end)]);
        xlabel('time (ms)'); ylabel('channel'); title('LFP');
        plot([0 0], ylim, '--r');
    end
end

end