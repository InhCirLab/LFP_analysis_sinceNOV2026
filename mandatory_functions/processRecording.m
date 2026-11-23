function o = processRecording(recPath, recName, P, runOnly)
% PROCESSRECORDING  Process ONE recording. Returns struct 'o'.
%
% Standalone function (extracted from LFPpipeline.m so sections can be run
% one at a time). Depends on: detectThetaEpochs, getLFP, getWavelet,
% extractEpochPower, storeRunPower, storeNonRunPower, buildRow,
% robust_zscore_matrix, saveStates.

o = struct('zAll',nan(1,3), 'zRun',nan(1,15), 'zNonRun',nan(1,15), ...
           'nR',NaN, 'totDur',NaN, 'meanDur',NaN, 'recLen',NaN, 'pct',NaN, ...
           'runPow',[], 'nonrunPow',[]);
o.msg = {};

cd(recPath);
statesFile = fullfile(recPath, [recName '.thetaEpochs6to12.states.mat']);

% Old version of nonTheta_nonRipples (to compare after the refresh)
oldNLNR = [];
if ~runOnly && isfile(statesFile)
    tmpOld = load(statesFile, 'thisthetaEpochs');
    if isfield(tmpOld.thisthetaEpochs, 'intervals_nonTheta_nonRipples')
        oldNLNR = tmpOld.thisthetaEpochs.intervals_nonTheta_nonRipples;
    end
end

te         = detectThetaEpochs(P.DETECT_OPTS{:});
powts      = te.timestamps(:);
run_int    = te.thetaRun.intervals;
nonrun_int = te.intervals_nonTheta_nonRipples;

if ~runOnly
    if isempty(oldNLNR)
        o.msg{end+1} = sprintf('nonTheta_nonRipples: no previous version to compare (%d intervals)', size(nonrun_int,1));
    elseif isequal(size(oldNLNR), size(nonrun_int)) && all(abs(oldNLNR(:) - nonrun_int(:)) < 1e-9)
        o.msg{end+1} = sprintf('nonTheta_nonRipples: SAME as previous version (%d intervals)', size(nonrun_int,1));
    else
        o.msg{end+1} = sprintf('nonTheta_nonRipples: UPDATED (%d -> %d intervals)', size(oldNLNR,1), size(nonrun_int,1));
    end
end

S = load(fullfile(recPath, [recName '.hippocampalLayers.channelinfo.mat']), 'hippocampalLayers');
hippocampalLayers = S.hippocampalLayers;
nCh = hippocampalLayers.channels(end);

% Get SLM, fallback if missing or NaN
if isfield(hippocampalLayers.layers{1,1}, 'slm') && ~isempty(hippocampalLayers.layers{1,1}.slm)
    slmCh = hippocampalLayers.layers{1,1}.slm;
    if isnan(slmCh)
        o.msg{end+1} = sprintf('WARN: SLM = NaN in %s. Using last channel (%d).', recName, nCh);
        slmCh = nCh;
    end
else
    o.msg{end+1} = sprintf('WARN: SLM not found in %s. Using last channel (%d).', recName, nCh);
    slmCh = nCh;
end
lfpT         = getLFP(slmCh, 'noPrompts', true);
samplingRate = lfpT.samplingRate;
lfpData      = double(lfpT.data(:,1));
lfpData      = removeLineNoise(lfpData, samplingRate);   % remove 50/100/150 Hz line noise

run_int = run_int((run_int(:,2) - run_int(:,1)) >= P.MIN_RUN_DUR, :);

[tW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, P.THETA_BAND(1),  P.THETA_BAND(2),  P.N_CYCLES, 0);
[lW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, P.LGAMMA_BAND(1), P.LGAMMA_BAND(2), P.N_CYCLES, 0);
[hW,~,~,~,~,~,~,~,~,~] = getWavelet(lfpData, samplingRate, P.HGAMMA_BAND(1), P.HGAMMA_BAND(2), P.N_CYCLES, 0);

rawT = median(tW,1);  rawL = median(lW,1);  rawH = median(hW,1);
tDB = pow2db(tW);  tDB_med = median(tDB,1);
lDB = pow2db(lW);  lDB_med = median(lDB,1);
hDB = pow2db(hW);  hDB_med = median(hDB,1);
tZ  = mean(zscore(tDB),1);  lZ = mean(zscore(lDB),1);  hZ = mean(zscore(hDB),1);
tZm = median(robust_zscore_matrix(tDB),1);
lZm = median(robust_zscore_matrix(lDB),1);
hZm = median(robust_zscore_matrix(hDB),1);

runPow = extractEpochPower(run_int, powts, rawT,rawL,rawH, tDB_med,lDB_med,hDB_med, tZ,lZ,hZ, tZm,lZm,hZm);
te     = storeRunPower(te, run_int, runPow);

nR   = size(run_int,1);
rDur = run_int(:,2) - run_int(:,1);

if ~runOnly
    nonrunPow = extractEpochPower(nonrun_int, powts, rawT,rawL,rawH, tDB_med,lDB_med,hDB_med, tZ,lZ,hZ, tZm,lZm,hZm);
    te        = storeNonRunPower(te, nonrun_int, nonrunPow);

    nNR  = size(nonrun_int,1);
    nDur = nonrun_int(:,2) - nonrun_int(:,1);

    o.zAll    = [median(tDB_med) median(lDB_med) median(hDB_med)];
    o.zRun    = buildRow(nR,  rDur, runPow);
    o.zNonRun = buildRow(nNR, nDur, nonrunPow);
    o.nR      = nR;
    o.totDur  = sum(rDur);
    o.meanDur = mean(rDur);
    o.recLen  = te.thetaRun.timestamps(end);
    o.pct     = o.totDur / o.recLen;
    o.runPow    = runPow;
    o.nonrunPow = nonrunPow;
end

saveStates(statesFile, te);
o.msg{end+1} = 'OK';
end
