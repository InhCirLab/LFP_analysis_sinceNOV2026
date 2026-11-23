function [groupeddata] = ps1to120_Paria(thisgrouplabel, mainFolder)
% PS1TO120  Compute 1-120 Hz LFP power spectra per hippocampal layer.
%
%  - Accepts mainFolder as an optional argument (otherwise opens uigetdir).
%  - delta_hires_indb / delta_hires_inuv fields per layer (pwelch, 1-3.5 Hz).
%  - If a layer (or, pyr, rad, slm) is missing, it is skipped and the rest
%    are processed (that session gets NaN for the missing layer).
%  - SERIAL: each session is processed one at a time (plain for loop). Slower
%    but easier to debug. (To restore parallel: change the 'for' loop marked
%    "SERIAL (was parfor)" below back to parfor and open a pool before calling.)
%  - Reads <session>.thetaEpochs6to12.states.mat (the one refreshed with new ripples).

% ── Main folder ───────────────────────────────────────────────────────────
if nargin < 2 || isempty(mainFolder)
    mainFolder = uigetdir();
end
if mainFolder == 0
    error('ps1to120: no folder selected. Aborted.');
end

folderInfo = dir(mainFolder);
subfolders = folderInfo([folderInfo.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
lfppaths   = {subfolders.name};
nSess      = numel(lfppaths);

DELTA_RANGE = [1 3.5];
ALL_LAYERS  = {'or','pyr','rad','slm'};
layerFieldMap = struct('or','oriens', 'pyr','pyramidal', 'rad','radiatum', 'slm','slm');

% ── Process sessions in parallel ──────────────────────────────────────────
res = cell(nSess, 1);
for ii = 1:nSess   % SERIAL (was parfor)
    res{ii} = processSession(fullfile(mainFolder, lfppaths{ii}), lfppaths{ii}, ...
                             ALL_LAYERS, layerFieldMap, DELTA_RANGE);
end

% ── Messages (workers do not print on the client) ─────────────────────────
for ii = 1:nSess
    if isempty(res{ii}), continue; end
    for m = 1:numel(res{ii}.msg), fprintf('  [%s] %s\n', lfppaths{ii}, res{ii}.msg{m}); end
end

% ── Assemble session x frequency accumulators ─────────────────────────────
okSess = find(~cellfun(@(r) isempty(r) || ~isfield(r,'f') || isempty(r.f), res));
if isempty(okSess)
    error('ps1to120: no session produced data.');
end
freqs          = res{okSess(1)}.f;
freqs_delta_hr = res{okSess(1)}.f_delta;

varList = {'allrec_db','allrec_uv','loc_db','loc_uv','nlnr_db','nlnr_uv','nl_db','nl_uv','delta_db','delta_uv'};
acc = struct();
for li = 1:numel(ALL_LAYERS)
    lname = ALL_LAYERS{li};
    for vi = 1:numel(varList)
        vname = varList{vi};
        % vector length (first session that has this layer)
        len = [];
        for k = okSess(:)'
            if isfield(res{k}.layers, lname)
                len = numel(res{k}.layers.(lname).(vname));
                break
            end
        end
        if isempty(len), continue; end   % layer never available
        M = nan(nSess, len);
        for k = okSess(:)'
            if isfield(res{k}.layers, lname)
                M(k, :) = res{k}.layers.(lname).(vname);
            end
        end
        acc.(lname).(vname) = M;
    end
end

% ── ID TABLE ──────────────────────────────────────────────────────────────
ssid      = lfppaths.';
ssexpcond = repmat(thisgrouplabel, nSess, 1);
SSidtable = table(ssid, ssexpcond, 'VariableNames', ["ID", "ExpCond"]);

% ── PACKAGE OUTPUT ─────────────────────────────────────────────────────────
groupeddata.idtable           = SSidtable;
groupeddata.freqs             = freqs;
groupeddata.freqs_delta_hires = freqs_delta_hr;
groupeddata.label             = thisgrouplabel;

for li = 1:length(ALL_LAYERS)
    lname = ALL_LAYERS{li};
    if ~isfield(acc, lname)
        groupeddata.(lname) = struct();
        continue
    end
    a = acc.(lname);
    groupeddata.(lname).allrec              = a.allrec_db;
    groupeddata.(lname).allrecinuv          = a.allrec_uv;
    groupeddata.(lname).locindb             = a.loc_db;
    groupeddata.(lname).locinuv             = a.loc_uv;
    groupeddata.(lname).nonlocnonrippleindb = a.nlnr_db;
    groupeddata.(lname).nonlocnonrippleinuv = a.nlnr_uv;
    groupeddata.(lname).nonlocindb          = a.nl_db;
    groupeddata.(lname).nonlocinuv          = a.nl_uv;
    groupeddata.(lname).delta_hires_indb    = a.delta_db;
    groupeddata.(lname).delta_hires_inuv    = a.delta_uv;
end
end


function r = processSession(sessPath, sessName, ALL_LAYERS, layerFieldMap, DELTA_RANGE)
% Process ONE session.
r = struct('f',[], 'f_delta',[], 'layers',struct());
r.msg = {};

cd(sessPath);

S = load(fullfile(sessPath, [sessName '.thetaEpochs6to12.states.mat']), 'thisthetaEpochs');
thetaEpochs = S.thisthetaEpochs;

% ── Available channels ─────────────────────────────────────────────────────
hippocampalLayers = getHippocampalLayers;
layerDef = hippocampalLayers.layers{1,1};

availableLayers = {};
chanMap = struct();
for li = 1:length(ALL_LAYERS)
    lname = ALL_LAYERS{li};
    fname = layerFieldMap.(lname);

    if ~isfield(layerDef, fname)
        r.msg{end+1} = sprintf('WARN: field "%s" does not exist, layer "%s" skipped.', fname, lname);
        continue
    end
    chNum = layerDef.(fname);
    if isempty(chNum) || ~isnumeric(chNum) || any(isnan(chNum)) || any(chNum <= 0) || any(chNum ~= floor(chNum))
        r.msg{end+1} = sprintf('WARN: invalid channel (%s = %s) for layer "%s", skipped.', fname, mat2str(chNum), lname);
        continue
    end
    availableLayers{end+1} = lname; %#ok<AGROW>
    chanMap.(lname) = chNum;
end

if isempty(availableLayers)
    r.msg{end+1} = 'WARN: no layer available, session skipped.';
    return
end

% ── LFP of available layers ────────────────────────────────────────────────
lfpData = struct();
for li = 1:length(availableLayers)
    lname = availableLayers{li};
    lfpData.(lname) = getLFP(chanMap.(lname), 'noPrompts', true);
    lfpData.(lname).data = removeLineNoise(lfpData.(lname).data, ...
                                           lfpData.(lname).samplingRate);  % 50/100/150 Hz
end
Fs = lfpData.(availableLayers{1}).samplingRate;

params.Fs     = Fs;
params.fpass  = [1 120];
params.tapers = [3 5];
params.pad    = 1;

thisIDints = thetaEpochs.thetaRun.newintervals;

% ── NON-LOCOMOTION NON-RIPPLES ─────────────────────────────────────────────
nlnr_int      = thetaEpochs.intervals_nonTheta_nonRipples;
dif_nlnr      = nlnr_int(:,2) - nlnr_int(:,1);
ntnr_int_keep = round(nlnr_int(dif_nlnr >= 1, :));

% ── NON-LOCOMOTION (INCLUDING RIPPLES) ─────────────────────────────────────
nl_int      = thetaEpochs.intervals_nonTheta;
dif_nl      = nl_int(:,2) - nl_int(:,1);
nl_int_keep = round(nl_int(dif_nl >= 1, :));

% ── High-resolution delta pwelch parameters ────────────────────────────────
nfft_hires = 2^nextpow2(Fs / 0.05);
noverlap   = floor(nfft_hires * 0.5);
win_hires  = hann(nfft_hires);
didx_hr    = [];

for li = 1:length(availableLayers)
    lname = availableLayers{li};
    lfp   = lfpData.(lname);

    % Multitaper spectrogram
    [Sx, ~, f] = mtspecgramc_fast(single(lfp.data), [2 1], params);
    Sx(Sx == 0) = NaN;
    Sdb = pow2db(Sx);

    if isempty(r.f), r.f = f; end

    % Clamp bounds for non-loc intervals
    checkend = size(Sx, 1);
    if ntnr_int_keep(end) > checkend
        ntnr_int_keep = ntnr_int_keep(1:end-1, :);
    end
    if ntnr_int_keep(1,1) == 0
        ntnr_int_keep(1,1) = 1;
    end
    if nl_int_keep(end) > checkend
        nl_int_keep = nl_int_keep(1:end-1, :);
    end
    if nl_int_keep(1,1) == 0
        nl_int_keep(1,1) = 1;
    end

    L = struct();
    L.allrec_db = median(Sdb, 1, 'omitnan');
    L.allrec_uv = median(Sx,  1, 'omitnan');
    L.loc_db    = computeLocalizedMedian(Sdb, thisIDints);
    L.loc_uv    = computeLocalizedMedian(Sx,  thisIDints);
    L.nlnr_db   = computeLocalizedMedian(Sdb, ntnr_int_keep);
    L.nlnr_uv   = computeLocalizedMedian(Sx,  ntnr_int_keep);
    L.nl_db     = computeLocalizedMedian(Sdb, nl_int_keep);
    L.nl_uv     = computeLocalizedMedian(Sx,  nl_int_keep);

    % High-resolution delta (pwelch)
    sig = double(lfp.data(:,1));
    [pxx, f_hr] = pwelch(sig, win_hires, noverlap, nfft_hires, Fs);
    if isempty(didx_hr)
        didx_hr = f_hr >= DELTA_RANGE(1) & f_hr <= DELTA_RANGE(2);
        if isempty(r.f_delta), r.f_delta = f_hr(didx_hr); end
    end
    pxx_delta  = pxx(didx_hr);
    L.delta_db = pow2db(pxx_delta(:))';
    L.delta_uv = pxx_delta(:)';

    % Force row vectors (1 x nFreq)
    fn = fieldnames(L);
    for k = 1:numel(fn), L.(fn{k}) = L.(fn{k})(:)'; end

    r.layers.(lname) = L;
end
r.msg{end+1} = 'OK';
end
