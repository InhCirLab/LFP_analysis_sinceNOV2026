function [groupeddata] = ps1to120(thisgrouplabel)

mainFolder = uigetdir();
folderInfo = dir(mainFolder);

subfolders = folderInfo([folderInfo.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));

lfppaths = {subfolders.name};

for ii = 1:length(lfppaths)

    thisfoldername = strcat(mainFolder, "\", lfppaths(ii));
    cd(thisfoldername)
    thetafilename = strcat(thisfoldername, "\", lfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(thetafilename);
    thetaEpochs   = thisthetaEpochs;
    hippocampalLayers = getHippocampalLayers;
    slmchannel = hippocampalLayers.layers{1,1}.slm;
    orchannel  = hippocampalLayers.layers{1,1}.oriens;
    pyrchannel = hippocampalLayers.layers{1,1}.pyramidal;
    radchannel = hippocampalLayers.layers{1,1}.radiatum;
    thisIDints = thetaEpochs.thetaRun.newintervals;

    slmlfpT = getLFP(slmchannel, 'noPrompts', true);
    orlfpT  = getLFP(orchannel,  'noPrompts', true);
    pyrlfpT = getLFP(pyrchannel, 'noPrompts', true);
    radlfpT = getLFP(radchannel, 'noPrompts', true);

    slmparams.Fs    = slmlfpT.samplingRate;
    slmparams.fpass = [1 120];
    slmparams.tapers = [3 5];
    slmparams.pad   = 1;
    orparams  = slmparams;
    pyrparams = slmparams;
    radparams = slmparams;

    % ── OR ────────────────────────────────────────────────────────────────
    [orS, ort, orf] = mtspecgramc_fast(single(orlfpT.data), [2 1], orparams);
    orS(orS == 0)   = NaN;
    orSdb           = pow2db(orS);

    orallrecthisIDps_inuv(ii, :) = median(orS,   1, 'omitnan');   % uV  (was here, now saved)
    orallrecthisIDps_indb(ii, :) = median(orSdb, 1, 'omitnan');   % dB  (new — median of dB)

    % ── PYR ───────────────────────────────────────────────────────────────
    [pyrS, pyrt, pyrf] = mtspecgramc_fast(single(pyrlfpT.data), [2 1], pyrparams);
    pyrS(pyrS == 0)    = NaN;
    pyrSdb             = pow2db(pyrS);

    pyrallrecthisIDps_inuv(ii, :) = median(pyrS,   1, 'omitnan');
    pyrallrecthisIDps_indb(ii, :) = median(pyrSdb, 1, 'omitnan');

    % ── RAD ───────────────────────────────────────────────────────────────
    [radS, radt, radf] = mtspecgramc_fast(single(radlfpT.data), [2 1], radparams);
    radS(radS == 0)    = NaN;
    radSdb             = pow2db(radS);

    radallrecthisIDps_inuv(ii, :) = median(radS,   1, 'omitnan');
    radallrecthisIDps_indb(ii, :) = median(radSdb, 1, 'omitnan');

    % ── SLM ───────────────────────────────────────────────────────────────
    [slmS, slmt, slmf] = mtspecgramc_fast(single(slmlfpT.data), [2 1], slmparams);
    slmS(slmS == 0)    = NaN;
    slmSdb             = pow2db(slmS);

    slmallrecthisIDps_inuv(ii, :) = median(slmS,   1, 'omitnan');
    slmallrecthisIDps_indb(ii, :) = median(slmSdb, 1, 'omitnan');

    % ── LOCOMOTION ────────────────────────────────────────────────────────
    orlocS_median_db(ii, :)  = computeLocalizedMedian(orSdb,  thisIDints);
    pyrlocS_median_db(ii, :) = computeLocalizedMedian(pyrSdb, thisIDints);
    radlocS_median_db(ii, :) = computeLocalizedMedian(radSdb, thisIDints);
    slmlocS_median_db(ii, :) = computeLocalizedMedian(slmSdb, thisIDints);

    orlocS_median_uv(ii, :)  = computeLocalizedMedian(orS,  thisIDints);
    pyrlocS_median_uv(ii, :) = computeLocalizedMedian(pyrS, thisIDints);
    radlocS_median_uv(ii, :) = computeLocalizedMedian(radS, thisIDints);
    slmlocS_median_uv(ii, :) = computeLocalizedMedian(slmS, thisIDints);

    % ── NON-LOCOMOTION NON-RIPPLES ────────────────────────────────────────
    nlnr_int       = thetaEpochs.intervals_nonTheta_nonRipples;
    dif_nlnr_int   = nlnr_int(:,2) - nlnr_int(:,1);
    ntnr_int_keep  = nlnr_int(dif_nlnr_int >= 1, :);
    ntnr_int_keep  = round(ntnr_int_keep);
    checkend       = size(slmS, 1);
    if ntnr_int_keep(end) > checkend
        ntnr_int_keep = ntnr_int_keep(1:end-1, :);
    end
    if ntnr_int_keep(1,1) == 0
        ntnr_int_keep(1,1) = 1;
    end

    ornlnr_median_db(ii, :)  = computeLocalizedMedian(orSdb,  ntnr_int_keep);
    pyrnlnr_median_db(ii, :) = computeLocalizedMedian(pyrSdb, ntnr_int_keep);
    radnlnr_median_db(ii, :) = computeLocalizedMedian(radSdb, ntnr_int_keep);
    slmnlnr_median_db(ii, :) = computeLocalizedMedian(slmSdb, ntnr_int_keep);

    ornlnr_median_uv(ii, :)  = computeLocalizedMedian(orS,  ntnr_int_keep);
    pyrnlnr_median_uv(ii, :) = computeLocalizedMedian(pyrS, ntnr_int_keep);
    radnlnr_median_uv(ii, :) = computeLocalizedMedian(radS, ntnr_int_keep);
    slmnlnr_median_uv(ii, :) = computeLocalizedMedian(slmS, ntnr_int_keep);

    % ── NON-LOCOMOTION (INCLUDING RIPPLES) ────────────────────────────────
    nl_int       = thisthetaEpochs.intervals_nonTheta;
    dif_nl_int   = nl_int(:,2) - nl_int(:,1);
    nl_int_keep  = nl_int(dif_nl_int >= 1, :);
    nl_int_keep  = round(nl_int_keep);
    if nl_int_keep(end) > checkend
        nl_int_keep = nl_int_keep(1:end-1, :);
    end
    if nl_int_keep(1,1) == 0
        nl_int_keep(1,1) = 1;
    end

    ornl_median_db(ii, :)  = computeLocalizedMedian(orSdb,  nl_int_keep);
    pyrnl_median_db(ii, :) = computeLocalizedMedian(pyrSdb, nl_int_keep);
    radnl_median_db(ii, :) = computeLocalizedMedian(radSdb, nl_int_keep);
    slmnl_median_db(ii, :) = computeLocalizedMedian(slmSdb, nl_int_keep);

    ornl_median_uv(ii, :)  = computeLocalizedMedian(orS,  nl_int_keep);
    pyrnl_median_uv(ii, :) = computeLocalizedMedian(pyrS, nl_int_keep);
    radnl_median_uv(ii, :) = computeLocalizedMedian(radS, nl_int_keep);
    slmnl_median_uv(ii, :) = computeLocalizedMedian(slmS, nl_int_keep);

end

freqs = orf;

% ── ID TABLE ──────────────────────────────────────────────────────────────
ssid       = lfppaths.';
ssexpcond  = repmat(thisgrouplabel, length(lfppaths), 1);
SSidtable  = table(ssid, ssexpcond, 'VariableNames', ["ID", "ExpCond"]);

% ── EMPAQUETAR OUTPUT ─────────────────────────────────────────────────────

groupeddata.idtable = SSidtable;
groupeddata.freqs   = freqs;
groupeddata.label   = thisgrouplabel;

% OR
groupeddata.or.allrec              = orallrecthisIDps_indb;    % dB
groupeddata.or.allrecinuv          = orallrecthisIDps_inuv;    % uV  ← NUEVO
groupeddata.or.locindb             = orlocS_median_db;
groupeddata.or.locinuv             = orlocS_median_uv;
groupeddata.or.nonlocnonrippleindb = ornlnr_median_db;
groupeddata.or.nonlocnonrippleinuv = ornlnr_median_uv;
groupeddata.or.nonlocindb          = ornl_median_db;
groupeddata.or.nonlocinuv          = ornl_median_uv;

% PYR
groupeddata.pyr.allrec              = pyrallrecthisIDps_indb;
groupeddata.pyr.allrecinuv          = pyrallrecthisIDps_inuv;  % ← NUEVO
groupeddata.pyr.locindb             = pyrlocS_median_db;
groupeddata.pyr.locinuv             = pyrlocS_median_uv;
groupeddata.pyr.nonlocnonrippleindb = pyrnlnr_median_db;
groupeddata.pyr.nonlocnonrippleinuv = pyrnlnr_median_uv;
groupeddata.pyr.nonlocindb          = pyrnl_median_db;
groupeddata.pyr.nonlocinuv          = pyrnl_median_uv;

% RAD
groupeddata.rad.allrec              = radallrecthisIDps_indb;
groupeddata.rad.allrecinuv          = radallrecthisIDps_inuv;  % ← NUEVO
groupeddata.rad.locindb             = radlocS_median_db;
groupeddata.rad.locinuv             = radlocS_median_uv;
groupeddata.rad.nonlocnonrippleindb = radnlnr_median_db;
groupeddata.rad.nonlocnonrippleinuv = radnlnr_median_uv;
groupeddata.rad.nonlocindb          = radnl_median_db;
groupeddata.rad.nonlocinuv          = radnl_median_uv;

% SLM
groupeddata.slm.allrec              = slmallrecthisIDps_indb;
groupeddata.slm.allrecinuv          = slmallrecthisIDps_inuv;  % ← NUEVO
groupeddata.slm.locindb             = slmlocS_median_db;
groupeddata.slm.locinuv             = slmlocS_median_uv;
groupeddata.slm.nonlocnonrippleindb = slmnlnr_median_db;
groupeddata.slm.nonlocnonrippleinuv = slmnlnr_median_uv;
groupeddata.slm.nonlocindb          = slmnl_median_db;
groupeddata.slm.nonlocinuv          = slmnl_median_uv;

clearvars -except groupeddata
end