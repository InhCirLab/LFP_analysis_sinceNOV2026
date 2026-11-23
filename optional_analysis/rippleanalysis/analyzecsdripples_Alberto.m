%% ========================================================================
%  PARTE VI — CSD DE RIPPLES (PARALLEL VERSION)
%  + VI.5: WAVEFORMS MEDIAS POR REGISTRO (todos los ripples juntos)
%  ========================================================================

clearvars; close all; clc;

%% ── VI.0  PARAMETROS ────────────────────────────────────────────────────
CSD_WIN_MS          = 100;
EXCLUDE_THETA       = true;
SMOOTH_TEMPORAL_WIN = 9;
SMOOTH_SPATIAL_WIN  = 5;
WAV_WIN_MS          = 150;
WAV_FREQS           = 100:1:300;
WAV_CYCLES          = 7;
WAV_SR_TARGET       = 1250; %#ok<NASGU>
WM_SAMPLE_RIP       = 200;   % nº de ripples individuales dibujados por registro en VI.5

%% ── Iniciar parallel pool ────────────────────────────────────────────────
if isempty(gcp('nocreate'))
    parpool('local');
end

%% ── VI.1  CONFIGURACION ──────────────────────────────────────────────────
nGroupsStr = inputdlg('Cuantos grupos experimentales?','CSD Ripples',1,{'2'});
if isempty(nGroupsStr), error('Abortado.'); end
nGroups = str2double(nGroupsStr{1});

groupFolders = cell(nGroups,1); groupLabels = cell(nGroups,1);
for g = 1:nGroups
    lbl = inputdlg(sprintf('Etiqueta grupo %d:',g),'',1,{sprintf('G%d',g)});
    if isempty(lbl), error('Abortado.'); end
    groupLabels{g} = strtrim(lbl{1});
    folder = uigetdir(pwd, sprintf('Carpeta grupo: %s', groupLabels{g}));
    if folder==0, error('Abortado.'); end
    groupFolders{g} = folder;
end
outFolder = uigetdir(groupFolders{1},'OUTPUT — Parte VI');
if outFolder==0, error('Abortado.'); end

figFolder      = fullfile(outFolder,'figures','csd');
figFolderAllCh = fullfile(outFolder,'figures','csd_allch');
figFolderWav   = fullfile(outFolder,'figures','wavelet');
for fd = {figFolder, figFolderAllCh, figFolderWav}
    if ~exist(fd{1},'dir'), mkdir(fd{1}); end
end

%% ── VI.2  BUCLE POR GRUPO ────────────────────────────────────────────────
allMeasures = table(); allWavMeas = table();

for g = 1:nGroups
    grpLabel   = groupLabels{g};
    mainFolder = groupFolders{g};

    entries  = dir(mainFolder); entries = entries([entries.isdir]);
    entries  = entries(~ismember({entries.name},{'.','..'}));
    sesNames = {entries.name}; nSess = numel(sesNames);
    assert(nSess>0,'No se encontraron subcarpetas en: %s',mainFolder);
    fprintf('\n==== Grupo "%s": %d sesion(es) ====\n',grpLabel,nSess);

    grpCSD_all = cell(nSess,1); grpCSD_rad = cell(nSess,1);
    grpCSD_slm = cell(nSess,1); grpCSD_oth = cell(nSess,1);
    grpWav_all = cell(nSess,1); grpWav_rad = cell(nSess,1);
    grpWav_slm = cell(nSess,1); grpMeanProf = nan(nSess,1);

    measCell   = cell(nSess,1);
    wavMeasCell= cell(nSess,1);

    wav_freqs   = WAV_FREQS;
    wav_cycles  = WAV_CYCLES;
    wav_win_ms  = WAV_WIN_MS;
    csd_win_ms  = CSD_WIN_MS;
    excl_theta  = EXCLUDE_THETA;

    parfor ii = 1:nSess
        sesName = sesNames{ii};
        sesPath = fullfile(mainFolder, sesName);
        fprintf('  [%d/%d] %s ...\n', ii, nSess, sesName);

        layersFile = fullfile(sesPath,[sesName '.hippocampalLayers.channelinfo.mat']);
        if ~isfile(layersFile)
            warning('Sin capas: %s. Saltando.',sesName); continue
        end
        S = load(layersFile,'hippocampalLayers');
        hippocampalLayers = S.hippocampalLayers;
        orCh  = getLayerCh(hippocampalLayers,'oriens');
        pyrCh = getLayerCh(hippocampalLayers,'pyramidal');
        radCh = getLayerCh(hippocampalLayers,'radiatum');
        slmCh = getLayerCh(hippocampalLayers,'slm');
        layerInfo = buildLayerInfo(orCh,pyrCh,radCh,slmCh);

        ripEvFile = fullfile(sesPath,[sesName '.ripples.events.mat']);
        if ~isfile(ripEvFile)
            warning('Sin ripples: %s. Saltando.',sesName); continue
        end
        try
            R = load(ripEvFile,'ripples');
            peaks      = R.ripples.peaks;
            timestamps = R.ripples.timestamps;
        catch ME
            warning('ripplesEvents no valido en %s: %s. Saltando.',sesName,ME.message);
            continue
        end
        if isempty(peaks), warning('Sin peaks en %s.',sesName); continue; end

        origDir = pwd;
        cd(sesPath);
        try; lfp = getLFP('all','noPrompts',true); catch; lfp = getLFP('all'); end
        cd(origDir);
        srLFP = lfp.samplingRate;

        if excl_theta
            thetaFile = fullfile(sesPath,[sesName '.thetaEpochs6to12.states.mat']);
            if isfile(thetaFile)
                T2 = load(thetaFile,'thisthetaEpochs');
                thetaInt = T2.thisthetaEpochs.thetaRun.newintervals;
                inTheta  = false(size(peaks));
                for tt = 1:size(thetaInt,1)
                    inTheta = inTheta | (peaks>=thetaInt(tt,1) & peaks<=thetaInt(tt,2));
                end
                peaks      = peaks(~inTheta);
                timestamps = timestamps(~inTheta,:);
                fprintf('    (theta excluido: %d ripples)\n',numel(peaks));
            end
        end
        if isempty(peaks), warning('Sin ripples tras theta: %s.',sesName); continue; end

        try
            [csd,~,csdEvents,~] = buzcode_eacheventCSD(lfp,peaks,...
                'plotCSD',false,'plotLFP',false);
        catch ME
            warning('CSD fallo en %s: %s',sesName,ME.message); continue
        end
        srCSD   = csd.samplingRate;
        nCSDchan= size(csdEvents,2);
        nRip    = size(csdEvents,3);
        if nRip < numel(peaks)
            peaks      = peaks(1:nRip);
            timestamps = timestamps(1:nRip,:);
        end

        halfWin = round((csd_win_ms/2)/1000*srCSD);
        winIdx  = csd.timestamps >= -halfWin & csd.timestamps <= halfWin;
        csdWin  = csdEvents(winIdx,:,:);

        csdSummary = nan(nRip,nCSDchan);
        for jj = 1:nRip
            durMs = diff(timestamps(jj,:))*1000;
            hw    = max(round((durMs/2)/1000*srCSD),1);
            idx   = csd.timestamps >= -hw & csd.timestamps <= hw;
            if ~any(idx), idx = true(size(csd.timestamps)); end
            csdSummary(jj,:) = mean(csdEvents(idx,:,jj),1);
        end
        [sinkAmp,sinkIdx] = min(csdSummary,[],2);

        radMask = false(nRip,1); slmMask = false(nRip,1);
        if ~isnan(pyrCh)&&~isnan(radCh)
            radMask = sinkIdx>=pyrCh+1 & sinkIdx<=radCh;
        elseif ~isnan(radCh)
            radMask = sinkIdx<=radCh;
        end
        if ~isnan(radCh), slmMask = sinkIdx>=radCh+1; end
        otherMask = ~radMask & ~slmMask;
        sinkType  = repmat("other",nRip,1);
        sinkType(radMask) = "rad"; sinkType(slmMask) = "slm";
        fprintf('    %s — RAD:%d SLM:%d Other:%d\n',...
            sesName,sum(radMask),sum(slmMask),sum(otherMask));

        grpCSD_all{ii} = makeCSDstruct(csdWin,winIdx,csd,nRip,layerInfo,srCSD);
        grpCSD_rad{ii} = makeCSDstruct(csdWin,winIdx,csd,sum(radMask),layerInfo,srCSD,radMask);
        grpCSD_slm{ii} = makeCSDstruct(csdWin,winIdx,csd,sum(slmMask),layerInfo,srCSD,slmMask);
        grpCSD_oth{ii} = makeCSDstruct(csdWin,winIdx,csd,sum(otherMask),layerInfo,srCSD,otherMask);
        grpMeanProf(ii) = mean(sinkIdx);

        thisID  = repmat(string(sesName),nRip,1);
        thisGrp = repmat(string(grpLabel),nRip,1);
        measCell{ii} = table(thisID,thisGrp,peaks(:),diff(timestamps,1,2)*1000,...
            sinkAmp,double(sinkIdx),sinkType,...
            'VariableNames',["ID","Group","PeakTime_s","Duration_ms",...
            "SinkAmplitude","SinkChannel","SinkType"]);

        wavMeasRows = nan(nRip,3);
        if ~isnan(pyrCh) && pyrCh>=1 && pyrCh<=size(lfp.data,2)
            lfpPyr      = double(lfp.data(:,pyrCh));
            halfWavSamp = round((wav_win_ms/2)/1000*srLFP);
            nSampWin    = 2*halfWavSamp+1;
            tsWav_ms    = linspace(-wav_win_ms/2, wav_win_ms/2, nSampWin);
            nFreqs      = numel(wav_freqs);
            nTotSmp     = size(lfpPyr,1);
            wavPow      = nan(nFreqs,nSampWin,nRip);
            waveforms   = nan(nSampWin,nRip);

            morletKernels = buildMorletKernels(wav_freqs, wav_cycles, srLFP);

            validRip = false(nRip,1);
            snippets = cell(nRip,1);
            for jj = 1:nRip
                peakSamp = round(peaks(jj)*srLFP);
                i1 = peakSamp-halfWavSamp; i2 = peakSamp+halfWavSamp;
                if i1<1||i2>nTotSmp, continue; end
                snippets{jj}   = lfpPyr(i1:i2);
                waveforms(:,jj)= lfpPyr(i1:i2);
                validRip(jj)   = true;
            end

            convResults = cell(nFreqs,1);
            for fi = 1:nFreqs
                convResults{fi} = conv(lfpPyr, morletKernels{fi}, 'same');
            end

            for jj = 1:nRip
                if ~validRip(jj), continue; end
                peakSamp = round(peaks(jj)*srLFP);
                i1 = peakSamp-halfWavSamp; i2 = peakSamp+halfWavSamp;
                for fi = 1:nFreqs
                    wavPow(fi,:,jj) = abs(convResults{fi}(i1:i2)).^2;
                end
            end

            grpWav_all{ii} = buildWavStruct(wavPow,waveforms,true(nRip,1),...
                wav_freqs,tsWav_ms,sesName,grpLabel,pyrCh);
            grpWav_rad{ii} = buildWavStruct(wavPow,waveforms,radMask,...
                wav_freqs,tsWav_ms,sesName,grpLabel,pyrCh);
            grpWav_slm{ii} = buildWavStruct(wavPow,waveforms,slmMask,...
                wav_freqs,tsWav_ms,sesName,grpLabel,pyrCh);

            centMask = abs(tsWav_ms) <= 20;
            for jj = 1:nRip
                pw = wavPow(:,:,jj); if all(isnan(pw(:))), continue; end
                pwC = pw(:,centMask); pwC(isnan(pwC))=0;
                [~,pkIdx] = max(pwC(:)); [fIdx,~] = ind2sub(size(pwC),pkIdx);
                wavMeasRows(jj,1) = wav_freqs(fIdx);
                wavMeasRows(jj,2) = max(pwC(:));
                tp = sum(pwC,2);
                if sum(tp)>0
                    wavMeasRows(jj,3) = sum(wav_freqs(:).*tp)/sum(tp);
                end
            end
            wavMeasCell{ii} = table(thisID,thisGrp,peaks(:),sinkType,...
                wavMeasRows(:,1),wavMeasRows(:,2),wavMeasRows(:,3),...
                'VariableNames',["ID","Group","PeakTime_s","SinkType",...
                "PeakFreq_Hz","PeakPower","WeightedFreq_Hz"]);
            fprintf('    Wavelet OK (%d ripples)\n',nRip);
        else
            warning('pyrCh no valido en %s.',sesName);
        end
    end  % fin parfor sesiones

    for ii = 1:nSess
        if ~isempty(measCell{ii}),    allMeasures = [allMeasures; measCell{ii}]; end    %#ok<AGROW>
        if ~isempty(wavMeasCell{ii}), allWavMeas  = [allWavMeas;  wavMeasCell{ii}]; end %#ok<AGROW>
    end

    set(0,'DefaultFigureVisible','off');
    for ii = 1:nSess
        sesName = sesNames{ii};
        if isempty(grpCSD_all{ii}), continue; end
        saveCSDFig(grpCSD_all{ii},grpCSD_rad{ii},grpCSD_slm{ii},...
            sesName,grpLabel,figFolder);
        saveCSDFigAllCh(grpCSD_all{ii},grpCSD_rad{ii},grpCSD_slm{ii},...
            sesName,grpLabel,figFolderAllCh,SMOOTH_TEMPORAL_WIN,SMOOTH_SPATIAL_WIN);
        if ~isempty(grpWav_all{ii})
            saveWavFig(grpWav_all{ii},grpWav_rad{ii},grpWav_slm{ii},...
                sesName,grpLabel,figFolderWav);
        end
    end
    close all; drawnow;
    set(0,'DefaultFigureVisible','on');

    csdGrpData = struct('label',grpLabel,'sessions',{sesNames},...
        'all',{grpCSD_all},'rad',{grpCSD_rad},'slm',{grpCSD_slm},...
        'other',{grpCSD_oth},'meanSinkChannel',grpMeanProf);
    save(fullfile(outFolder,sprintf('csd_ripples_%s.mat',grpLabel)),'csdGrpData','-mat');

    wavGrpData = struct('label',grpLabel,'sessions',{sesNames},...
        'all',{grpWav_all},'rad',{grpWav_rad},'slm',{grpWav_slm},...
        'freqs',WAV_FREQS,'winMs',WAV_WIN_MS,'cycles',WAV_CYCLES);
    save(fullfile(outFolder,sprintf('wavelet_ripples_%s.mat',grpLabel)),'wavGrpData','-mat');

    saveGroupCSDFig(grpCSD_all,grpCSD_rad,grpCSD_slm,grpLabel,sesNames,figFolder);
    saveGroupCSDFigAllCh(grpCSD_all,grpCSD_rad,grpCSD_slm,grpLabel,sesNames,...
        figFolderAllCh,SMOOTH_TEMPORAL_WIN,SMOOTH_SPATIAL_WIN);
    saveGroupWavFig(grpWav_all,grpWav_rad,grpWav_slm,grpLabel,sesNames,figFolderWav);
    fprintf('  Grupo "%s" completado.\n',grpLabel);
end

%% ── VI.3  EXCEL CSD ──────────────────────────────────────────────────────
fprintf('\n[VI.3] Exportando Excel CSD...\n');
writetable(allMeasures,fullfile(outFolder,'csd_ripples_measures.xlsx'));
for g = 1:nGroups
    msk = strcmp(allMeasures.Group,groupLabels{g});
    writetable(allMeasures(msk,:),fullfile(outFolder,...
        sprintf('csd_ripples_measures_%s.xlsx',groupLabels{g})));
end

sinkTypes  = {'all','rad','slm'};
sinkLabels = {'ALL','RAD','SLM'};

for g = 1:nGroups
    grpLabel = groupLabels{g};
    matFile  = fullfile(outFolder,sprintf('csd_ripples_%s.mat',grpLabel));
    if ~isfile(matFile), continue; end
    load(matFile,'csdGrpData');
    csdCells  = {csdGrpData.all, csdGrpData.rad, csdGrpData.slm};
    sesNamesG = csdGrpData.sessions;

    for st = 1:3
        fileB = fullfile(outFolder,sprintf('csd_bysession_%s.xlsx',grpLabel));
        tOut  = buildCSDExcelBySess(csdCells{st},sesNamesG,grpLabel,sinkTypes{st});
        if ~isempty(tOut), writetable(tOut,fileB,'Sheet',sinkLabels{st}); end

        fileC = fullfile(outFolder,sprintf('csd_group_%s.xlsx',grpLabel));
        tOut  = buildCSDExcelGroup(csdCells{st},grpLabel,sinkTypes{st});
        if ~isempty(tOut), writetable(tOut,fileC,'Sheet',sinkLabels{st}); end

        fileH = fullfile(outFolder,sprintf('csd_heatmap_%s.xlsx',grpLabel));
        tOut  = buildCSDHeatmapExcelBySess(csdCells{st},sesNamesG,grpLabel,sinkTypes{st},...
            SMOOTH_TEMPORAL_WIN,SMOOTH_SPATIAL_WIN);
        if ~isempty(tOut), writetable(tOut,fileH,'Sheet',sinkLabels{st}); end

        fileHG = fullfile(outFolder,sprintf('csd_heatmap_group_%s.xlsx',grpLabel));
        tOut   = buildCSDHeatmapExcelGroup(csdCells{st},grpLabel,sinkTypes{st},...
            SMOOTH_TEMPORAL_WIN,SMOOTH_SPATIAL_WIN);
        if ~isempty(tOut), writetable(tOut,fileHG,'Sheet',sinkLabels{st}); end

        % ── NUEVO: waveforms CSD allch por sesion ─────────────────────────
        fileWF = fullfile(outFolder,sprintf('csd_allch_waveforms_%s.xlsx',grpLabel));
        tOut   = buildCSDAllChWaveformsBySess(csdCells{st},sesNamesG,grpLabel,sinkTypes{st});
        if ~isempty(tOut), writetable(tOut,fileWF,'Sheet',sinkLabels{st}); end

        % ── NUEVO: waveforms CSD allch media grupal ───────────────────────
        fileWFG = fullfile(outFolder,sprintf('csd_allch_waveforms_group_%s.xlsx',grpLabel));
        tOut    = buildCSDAllChWaveformsGroup(csdCells{st},grpLabel,sinkTypes{st});
        if ~isempty(tOut), writetable(tOut,fileWFG,'Sheet',sinkLabels{st}); end
    end
    fprintf('  [%s] Excel CSD completados.\n',grpLabel);
end

fileFig = fullfile(outFolder,'csd_plot_data.xlsx');
if isfile(fileFig), delete(fileFig); end
for st = 1:3
    tOut = buildCSDExcelAllGroups(outFolder,groupLabels,sinkTypes{st});
    if ~isempty(tOut), writetable(tOut,fileFig,'Sheet',sinkLabels{st}); end
end

%% ── VI.4  EXCEL WAVELET ──────────────────────────────────────────────────
fprintf('\n[VI.4] Exportando Excel Wavelet...\n');
if ~isempty(allWavMeas)
    writetable(allWavMeas,fullfile(outFolder,'wavelet_measures.xlsx'));
end

for g = 1:nGroups
    grpLabel = groupLabels{g};
    matFile  = fullfile(outFolder,sprintf('wavelet_ripples_%s.mat',grpLabel));
    if ~isfile(matFile), continue; end
    load(matFile,'wavGrpData');
    wavCells  = {wavGrpData.all, wavGrpData.rad, wavGrpData.slm};
    sesNamesG = wavGrpData.sessions;
    freqs     = wavGrpData.freqs;

    for st = 1:3
        fileF = fullfile(outFolder,sprintf('wavelet_bysession_%s.xlsx',grpLabel));
        tOut  = buildWavExcelBySess(wavCells{st},sesNamesG,grpLabel,sinkTypes{st},freqs);
        if ~isempty(tOut), writetable(tOut,fileF,'Sheet',sinkLabels{st}); end

        fileG = fullfile(outFolder,sprintf('wavelet_group_%s.xlsx',grpLabel));
        tOut  = buildWavExcelGroup(wavCells{st},grpLabel,sinkTypes{st},freqs);
        if ~isempty(tOut), writetable(tOut,fileG,'Sheet',sinkLabels{st}); end

        fileW = fullfile(outFolder,sprintf('waveforms_%s.xlsx',grpLabel));
        tOut  = buildWaveformsExcel(wavCells{st},sesNamesG,grpLabel,sinkTypes{st});
        if ~isempty(tOut), writetable(tOut,fileW,'Sheet',sinkLabels{st}); end

        fileWH = fullfile(outFolder,sprintf('wavelet_heatmap_%s.xlsx',grpLabel));
        tOut   = buildWavHeatmapExcelBySess(wavCells{st},sesNamesG,grpLabel,sinkTypes{st},...
            freqs,SMOOTH_TEMPORAL_WIN);
        if ~isempty(tOut), writetable(tOut,fileWH,'Sheet',sinkLabels{st}); end

        fileWHG = fullfile(outFolder,sprintf('wavelet_heatmap_group_%s.xlsx',grpLabel));
        tOut    = buildWavHeatmapExcelGroup(wavCells{st},grpLabel,sinkTypes{st},...
            freqs,SMOOTH_TEMPORAL_WIN);
        if ~isempty(tOut), writetable(tOut,fileWHG,'Sheet',sinkLabels{st}); end

        fileWS = fullfile(outFolder,sprintf('waveforms_smoothed_%s.xlsx',grpLabel));
        tOut   = buildWaveformsSmoothedExcel(wavCells{st},sesNamesG,grpLabel,sinkTypes{st},...
            SMOOTH_TEMPORAL_WIN);
        if ~isempty(tOut), writetable(tOut,fileWS,'Sheet',sinkLabels{st}); end
    end
    fprintf('  [%s] Excel Wavelet completados.\n',grpLabel);
end

fileWavFig = fullfile(outFolder,'wavelet_plot_data.xlsx');
if isfile(fileWavFig), delete(fileWavFig); end
for st = 1:3
    tOut = buildWavExcelAllGroups(outFolder,groupLabels,sinkTypes{st});
    if ~isempty(tOut), writetable(tOut,fileWavFig,'Sheet',sinkLabels{st}); end
end

%% ── VI.5  WAVEFORMS MEDIAS POR REGISTRO (TODOS LOS RIPPLES JUNTOS) ───────
%  LFP del canal piramidal alineado al pico del ripple, promediado por
%  registro con TODOS los ripples (sin separar por sinkType).
%  Salidas:
%    waveforms_mean_allripples.xlsx   (una hoja por tabla)
%    csv_waveforms/*.csv              (mismas tablas en CSV)
%    waveforms_mean_allripples.mat    (struct WM)
%    figures/waveform_mean/*          (.png, .fig, .pdf vectorial)
fprintf('\n[VI.5] Waveforms medias por registro (todos los ripples)...\n');
figFolderWM = fullfile(outFolder,'figures','waveform_mean');
csvFolderWM = fullfile(outFolder,'csv_waveforms');
for fd = {figFolderWM, csvFolderWM}
    if ~exist(fd{1},'dir'), mkdir(fd{1}); end
end

WM = collectMeanWaveforms(outFolder, groupLabels, SMOOTH_TEMPORAL_WIN, WM_SAMPLE_RIP);

if isempty(WM.rec)
    warning('No hay waveforms para exportar en VI.5.');
else
    TW = buildWMTables(WM, groupLabels);

    xlsWM = fullfile(outFolder,'waveforms_mean_allripples.xlsx');
    if isfile(xlsWM), delete(xlsWM); end
    tblNames = fieldnames(TW);
    for k = 1:numel(tblNames)
        writetable(TW.(tblNames{k}), xlsWM, 'Sheet', tblNames{k});
        writetable(TW.(tblNames{k}), fullfile(csvFolderWM, ...
            ['waveforms_mean_allripples_' tblNames{k} '.csv']));
    end
    save(fullfile(outFolder,'waveforms_mean_allripples.mat'),'WM','TW','-mat');

    set(0,'DefaultFigureVisible','off');
    for r = 1:numel(WM.rec)
        saveWMSessFig(WM.rec(r), WM.tsMs, figFolderWM);
    end
    saveWMGroupFig(WM, groupLabels, figFolderWM);
    saveWMCompareFig(TW.Group_MeanSEM, groupLabels, figFolderWM);
    close all; drawnow;
    set(0,'DefaultFigureVisible','on');

    fprintf('  %d registros exportados (xlsx + csv + figuras).\n', numel(WM.rec));
end

fprintf('\n[Parte VI] Completada. Resultados en: %s\n',outFolder);

%% =========================================================================
%  FUNCIONES LOCALES
%% =========================================================================

function ch = getLayerCh(hippocampalLayers, fieldName)
    try; ch=hippocampalLayers.layers{1,1}.(fieldName);
         if isempty(ch), ch=NaN; end; ch=ch(1);
    catch; ch=NaN; end
end

function info = buildLayerInfo(orCh, pyrCh, radCh, slmCh)
    rawCh=[orCh,pyrCh,radCh,slmCh]; rawLbl={'SO','SP','SR','SLM'};
    valid=~isnan(rawCh);
    info.channels=rawCh(valid); info.labels=rawLbl(valid); info.nLayers=sum(valid);
end

function s = makeCSDstruct(csdWin, winIdx, csd, nRip, layerInfo, srCSD, mask)
    if nargin<7, mask=true(nRip,1); end
    if sum(mask)==0, s=[]; return; end
    sub=csdWin(:,:,mask);
    s.mean=mean(sub,3); s.sem=std(sub,0,3)/sqrt(sum(mask));
    s.timestamps=csd.timestamps(winIdx); s.nRip=sum(mask);
    s.layerInfo=layerInfo; s.samplingRate=srCSD;
end

function out = smoothCSDdata(data, tempWin, spatWin)
    out=data;
    if tempWin>1
        gw=gausswin(tempWin); gw=gw/sum(gw);
        for ch=1:size(out,2)
            col=out(:,ch); col(isnan(col))=0;
            out(:,ch)=filtfilt(gw,1,col);
        end
    end
    if spatWin>1, out=movmean(out,spatWin,2,'omitnan'); end
end

function out = smoothWaveform(data, tempWin)
    if tempWin<=1, out=data; return; end
    gw=gausswin(tempWin); gw=gw/sum(gw);
    wasNaN=isnan(data); data(wasNaN)=0;
    out=filtfilt(gw,1,data(:)); out(wasNaN)=NaN;
end

function out = smoothWavPow(meanPow, tempWin)
    out=meanPow;
    if tempWin<=1, return; end
    gw=gausswin(tempWin); gw=gw/sum(gw);
    for fi=1:size(out,1)
        row=out(fi,:); row(isnan(row))=0;
        out(fi,:)=filtfilt(gw,1,row);
    end
end

function kernels = buildMorletKernels(freqs, nCycles, sr)
    nFreqs=numel(freqs); kernels=cell(nFreqs,1);
    for fi=1:nFreqs
        f=freqs(fi); sigma=nCycles/(2*pi*f);
        tWin=-3*sigma:1/sr:3*sigma;
        kern=(exp(-tWin.^2/(2*sigma^2)).*exp(2i*pi*f*tWin));
        kernels{fi}=kern/sum(abs(kern));
    end
end

function s = buildWavStruct(wavPow, waveforms, mask, freqs, tsMs, sesName, grpLabel, pyrCh)
    nRipSub=sum(mask);
    if nRipSub==0, s=[]; return; end
    subPow=wavPow(:,:,mask); subWav=waveforms(:,mask);
    s.meanPow=mean(subPow,3,'omitnan'); s.semPow=std(subPow,0,3,'omitnan')/sqrt(nRipSub);
    s.meanWaveform=mean(subWav,2,'omitnan'); s.semWaveform=std(subWav,0,2,'omitnan')/sqrt(nRipSub);
    s.allWaveforms=subWav; s.freqs=freqs(:); s.tsMs=tsMs(:);
    s.nRip=nRipSub; s.session=sesName; s.group=grpLabel; s.pyrCh=pyrCh;
end

function out = stackCSDCells(csdCell, validMask)
    idxV=find(validMask); nS=numel(idxV);
    layMats=cellfun(@(c) extractLayerData(c),csdCell(idxV),'UniformOutput',false);
    minR=min(cellfun(@(m)size(m,1),layMats)); minC=min(cellfun(@(m)size(m,2),layMats));
    stk=nan(minR,minC,nS);
    for k=1:nS; m=layMats{k}; stk(:,:,k)=m(1:minR,1:minC); end
    ts=csdCell{idxV(1)}.timestamps; ts=ts(1:minR);
    info=[];
    if isfield(csdCell{idxV(1)},'layerInfo')
        info=csdCell{idxV(1)}.layerInfo;
        if info.nLayers>minC, info.channels=info.channels(1:minC); info.labels=info.labels(1:minC); info.nLayers=minC; end
    end
    out.mean=mean(stk,3,'omitnan'); out.sem=std(stk,0,3)/sqrt(nS);
    out.timestamps=ts; out.nRip=sum(cellfun(@(c)c.nRip,csdCell(idxV)));
    out.layerInfo=info;
    if isfield(csdCell{idxV(1)},'samplingRate'), out.samplingRate=csdCell{idxV(1)}.samplingRate; end
end

function out = stackCSDCellsAllCh(csdCell, validMask)
    idxV=find(validMask); nS=numel(idxV);
    minR=min(cellfun(@(c)size(c.mean,1),csdCell(idxV)));
    minC=min(cellfun(@(c)size(c.mean,2),csdCell(idxV)));
    stk=nan(minR,minC,nS);
    for k=1:nS; m=csdCell{idxV(k)}.mean; stk(:,:,k)=m(1:minR,1:minC); end
    ts=csdCell{idxV(1)}.timestamps; ts=ts(1:minR);
    info=[]; if isfield(csdCell{idxV(1)},'layerInfo'), info=csdCell{idxV(1)}.layerInfo; end
    out.mean=mean(stk,3,'omitnan'); out.sem=std(stk,0,3)/sqrt(nS);
    out.timestamps=ts; out.nRip=sum(cellfun(@(c)c.nRip,csdCell(idxV)));
    out.layerInfo=info;
    if isfield(csdCell{idxV(1)},'samplingRate'), out.samplingRate=csdCell{idxV(1)}.samplingRate; end
end

function out = stackWavCells(wavCell, validMask)
    idxV=find(validMask); nS=numel(idxV);
    ref=wavCell{idxV(1)}; nF=numel(ref.freqs); nT=numel(ref.tsMs);
    stkP=nan(nF,nT,nS); stkW=nan(nT,nS); allWav=[];
    for k=1:nS
        c=wavCell{idxV(k)};
        stkP(1:min(nF,size(c.meanPow,1)),1:min(nT,size(c.meanPow,2)),k)=...
            c.meanPow(1:min(nF,size(c.meanPow,1)),1:min(nT,size(c.meanPow,2)));
        stkW(1:min(nT,numel(c.meanWaveform)),k)=c.meanWaveform(1:min(nT,numel(c.meanWaveform)));
        allWav=[allWav,c.allWaveforms]; %#ok<AGROW>
    end
    out.meanPow=mean(stkP,3,'omitnan'); out.semPow=std(stkP,0,3,'omitnan')/sqrt(nS);
    out.meanWaveform=mean(stkW,2,'omitnan'); out.semWaveform=std(stkW,0,2,'omitnan')/sqrt(nS);
    out.allWaveforms=allWav; out.freqs=ref.freqs; out.tsMs=ref.tsMs;
    out.nRip=sum(cellfun(@(c)c.nRip,wavCell(idxV)));
    out.group=ref.group; out.pyrCh=ref.pyrCh;
end

function layerData = extractLayerData(csdStruct)
    data=csdStruct.mean;
    if ~isfield(csdStruct,'layerInfo')||isempty(csdStruct.layerInfo)||csdStruct.layerInfo.nLayers==0
        layerData=data; return; end
    info=csdStruct.layerInfo; nChTot=size(data,2); nLay=info.nLayers; T=size(data,1);
    layerData=nan(T,nLay);
    for li=1:nLay; ch=info.channels(li); if ch>=1&&ch<=nChTot, layerData(:,li)=data(:,ch); end; end
end

% ══════════════════════════════════════════════════════════════════════════
%  VI.5 — WAVEFORMS MEDIAS POR REGISTRO (todos los ripples juntos)
% ══════════════════════════════════════════════════════════════════════════

% Lee wavelet_ripples_<grupo>.mat y calcula, por registro, la waveform media
% (y SEM entre ripples) de TODOS los ripples validos, en una rejilla temporal comun.
function WM = collectMeanWaveforms(outFolder, groupLabels, stW, nSample)
    WM.tsMs = [];
    WM.rec  = [];
    for g = 1:numel(groupLabels)
        mf = fullfile(outFolder,sprintf('wavelet_ripples_%s.mat',groupLabels{g}));
        if ~isfile(mf), continue; end
        S  = load(mf,'wavGrpData');
        wc = S.wavGrpData.all; sn = S.wavGrpData.sessions;
        for ii = 1:numel(wc)
            w = wc{ii};
            if isempty(w) || ~isfield(w,'allWaveforms'), continue; end
            t  = w.tsMs(:);
            if isempty(WM.tsMs), WM.tsMs = t; end
            wf = w.allWaveforms;                       % [T x nRip]
            wf = wf(:, any(~isnan(wf),1));             % solo ripples validos
            if isempty(wf), continue; end
            wf = interpToGrid(t, wf, WM.tsMs);         % robusto si cambia srLFP
            nR = size(wf,2);
            m  = mean(wf,2,'omitnan');
            s  = std(wf,0,2,'omitnan')/sqrt(nR);
            idx = round(linspace(1,nR,min(nR,nSample)));
            r = struct('Group',string(groupLabels{g}),'Session',string(sn{ii}),...
                'nRip',nR,'mean',m,'sem',s,'meanSm',smoothWaveform(m,stW),...
                'sample',wf(:,idx));
            if isempty(WM.rec), WM.rec = r; else, WM.rec(end+1) = r; end %#ok<AGROW>
        end
    end
end

function wfOut = interpToGrid(t, wf, tRef)
    if numel(t)==numel(tRef) && all(abs(t(:)-tRef(:))<1e-6)
        wfOut = wf; return;
    end
    wfOut = interp1(t(:), wf, tRef(:), 'linear', NaN);
end

% Tablas listas para copiar/pegar. Nombre de campo = nombre de hoja Excel.
function TW = buildWMTables(WM, groupLabels)
    rec = WM.rec; nS = numel(rec); ts = WM.tsMs(:); T = numel(ts);
    grp = string({rec.Group})'; ses = string({rec.Session})'; nR = [rec.nRip]';
    M = nan(nS,T); S = nan(nS,T); Ms = nan(nS,T);
    for k = 1:nS
        M(k,:)  = rec(k).mean(:)';
        S(k,:)  = rec(k).sem(:)';
        Ms(k,:) = rec(k).meanSm(:)';
    end
    colN = string(arrayfun(@(t)sprintf('t_%.3fms',t),ts','UniformOutput',false));
    meta = table(grp,ses,nR,'VariableNames',{'Group','Session','nRipples'});

    % 1) Filas = registros, columnas = tiempo
    TW.Mean_bySession       = [meta, array2table(M, 'VariableNames',colN)];
    TW.SEM_bySession        = [meta, array2table(S, 'VariableNames',colN)];
    TW.MeanSmooth_bySession = [meta, array2table(Ms,'VariableNames',colN)];

    % 2) Filas = tiempo, columnas = registros (directo a Prism/Origin/Excel)
    vn = cellstr(matlab.lang.makeValidName(grp + "_" + ses));
    vn = matlab.lang.makeUniqueStrings(vn(:), {'Time_ms'});
    vn = [{'Time_ms'}; vn(:)];
    TW.Mean_TimeXSession       = array2table([ts, M'],  'VariableNames',vn);
    TW.MeanSmooth_TimeXSession = array2table([ts, Ms'], 'VariableNames',vn);

    % 3) Formato largo (R / Python / GraphPad)
    TW.Long = table(repelem(grp,T,1), repelem(ses,T,1), repelem(nR,T,1), ...
        repmat(ts,nS,1), reshape(M',[],1), reshape(S',[],1), reshape(Ms',[],1), ...
        'VariableNames',{'Group','Session','nRipples','Time_ms','Mean','SEM','MeanSmoothed'});

    % 4) Media +- SEM por grupo (n = registros)
    gT = table(ts,'VariableNames',{'Time_ms'});
    gl = strings(0,1); gn = zeros(0,1); gr = zeros(0,1);
    for g = 1:numel(groupLabels)
        msk = grp==string(groupLabels{g}); if ~any(msk), continue; end
        n  = sum(msk); nm = matlab.lang.makeValidName(groupLabels{g});
        gT.([nm '_mean'])   = mean(M(msk,:),1,'omitnan')';
        gT.([nm '_sem'])    = (std(M(msk,:),0,1,'omitnan')/sqrt(n))';
        gT.([nm '_meanSm']) = mean(Ms(msk,:),1,'omitnan')';
        gT.([nm '_semSm'])  = (std(Ms(msk,:),0,1,'omitnan')/sqrt(n))';
        gl(end+1,1) = string(groupLabels{g}); %#ok<AGROW>
        gn(end+1,1) = n;                      %#ok<AGROW>
        gr(end+1,1) = sum(nR(msk));           %#ok<AGROW>
    end
    TW.Group_MeanSEM = gT;
    TW.Group_Info    = table(gl,gn,gr,'VariableNames',{'Group','nSessions','nRipplesTotal'});
end

function saveWMSessFig(rec, ts, figFolder)
    t = ts(:); m = rec.mean(:); s = rec.sem(:);
    fig = figure('Color','w','Position',[100 100 900 380]);
    ax1 = subplot(1,2,1); hold(ax1,'on');
    plot(ax1,t,rec.sample,'Color',[0.8 0.8 0.8],'LineWidth',0.3);
    plot(ax1,t,m,'r-','LineWidth',1.8);
    xline(ax1,0,'k--');
    xlim(ax1,[t(1) t(end)]); box(ax1,'off'); set(ax1,'TickDir','out','FontSize',9);
    xlabel(ax1,'Time from ripple peak (ms)'); ylabel(ax1,'LFP pyrCh (raw units)');
    title(ax1,sprintf('Ripples individuales (max %d) + media',size(rec.sample,2)),'FontSize',10);
    ax2 = subplot(1,2,2); hold(ax2,'on');
    fill(ax2,[t;flipud(t)],[m-s;flipud(m+s)],[1 0.7 0.7],'EdgeColor','none','FaceAlpha',0.6);
    h1 = plot(ax2,t,m,'r-','LineWidth',1.8);
    h2 = plot(ax2,t,rec.meanSm,'k-','LineWidth',1.0);
    xline(ax2,0,'k--');
    xlim(ax2,[t(1) t(end)]); box(ax2,'off'); set(ax2,'TickDir','out','FontSize',9);
    xlabel(ax2,'Time from ripple peak (ms)'); ylabel(ax2,'LFP pyrCh (raw units)');
    legend(ax2,[h1 h2],{'Media \pm SEM','Suavizada'},'Box','off','Location','best');
    title(ax2,sprintf('Media (n=%d ripples)',rec.nRip),'FontSize',10);
    sgtitle(sprintf('%s — %s',rec.Group,rec.Session),'Interpreter','none','FontSize',11,'FontWeight','bold');
    fname = fullfile(figFolder,sprintf('waveform_mean_%s_%s',rec.Group,rec.Session));
    exportWMFig(fig,fname); close(fig);
end

function saveWMGroupFig(WM, groupLabels, figFolder)
    rec = WM.rec; grp = string({rec.Group}); ts = WM.tsMs(:); nG = numel(groupLabels);
    cols = lines(max(nG,1));
    fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 min(0.3*nG,0.95) 0.8]);
    sgtitle('Waveform media por registro — todos los ripples','FontSize',12,'FontWeight','bold');
    for g = 1:nG
        msk = find(grp==string(groupLabels{g})); if isempty(msk), continue; end
        Mg = cell2mat(arrayfun(@(k)rec(k).mean(:)',msk(:),'UniformOutput',false)); % [nS x T]
        ax1 = subplot(2,nG,g); hold(ax1,'on');
        plot(ax1,ts,Mg','Color',cols(g,:)*0.5+0.5,'LineWidth',0.6);
        plot(ax1,ts,mean(Mg,1,'omitnan'),'Color',cols(g,:),'LineWidth',2);
        xline(ax1,0,'k--');
        xlim(ax1,[ts(1) ts(end)]); box(ax1,'off'); set(ax1,'TickDir','out','FontSize',9);
        xlabel(ax1,'Time from ripple peak (ms)'); ylabel(ax1,'LFP pyrCh (raw units)');
        title(ax1,sprintf('%s (n=%d registros)',groupLabels{g},numel(msk)),'Interpreter','none','FontSize',10);
        ax2 = subplot(2,nG,nG+g);
        imagesc(ax2,ts,1:numel(msk),Mg); colormap(ax2,divergingCmap());
        cl = max(abs(Mg(:)),[],'omitnan'); if cl>0, set(ax2,'CLim',cl*[-1 1]); end
        cb = colorbar(ax2); cb.Label.String = 'LFP'; cb.FontSize = 8;
        xline(ax2,0,'k--');
        set(ax2,'YTick',1:numel(msk),'YTickLabel',cellstr(string({rec(msk).Session})),...
            'TickLabelInterpreter','none','FontSize',7,'TickDir','out');
        xlabel(ax2,'Time from ripple peak (ms)'); ylabel(ax2,'Registro');
    end
    exportWMFig(fig,fullfile(figFolder,'waveform_mean_groups_overview')); close(fig);
end

function saveWMCompareFig(gT, groupLabels, figFolder)
    ts = gT.Time_ms; cols = lines(max(numel(groupLabels),1));
    sfx = {'','Sm'}; ttl = {'Media cruda \pm SEM (entre registros)','Suavizada \pm SEM'};
    fig = figure('Color','w','Position',[100 100 900 380]);
    for p = 1:2
        ax = subplot(1,2,p); hold(ax,'on'); h = []; lg = {};
        for g = 1:numel(groupLabels)
            nm = matlab.lang.makeValidName(groupLabels{g});
            if ~ismember([nm '_mean'],gT.Properties.VariableNames), continue; end
            m = gT.([nm '_mean' sfx{p}]); s = gT.([nm '_sem' sfx{p}]);
            fill(ax,[ts;flipud(ts)],[m-s;flipud(m+s)],cols(g,:),'FaceAlpha',0.25,...
                'EdgeColor','none','HandleVisibility','off');
            h(end+1) = plot(ax,ts,m,'Color',cols(g,:),'LineWidth',1.8); %#ok<AGROW>
            lg{end+1} = groupLabels{g};                                   %#ok<AGROW>
        end
        xline(ax,0,'k--','HandleVisibility','off');
        if ~isempty(h), legend(ax,h,lg,'Interpreter','none','Box','off'); end
        xlim(ax,[ts(1) ts(end)]); box(ax,'off'); set(ax,'TickDir','out','FontSize',9);
        xlabel(ax,'Time from ripple peak (ms)'); ylabel(ax,'LFP pyrCh (raw units)');
        title(ax,ttl{p},'FontSize',10);
    end
    exportWMFig(fig,fullfile(figFolder,'waveform_mean_groups_compare')); close(fig);
end

function exportWMFig(fig, fname)
    saveas(fig,[fname '.png']); savefig(fig,[fname '.fig']);
    try; exportgraphics(fig,[fname '.pdf'],'ContentType','vector'); catch; end
end

function cmap = divergingCmap()
    n=256; h=n/2;
    cmap=[[linspace(0,1,h)',linspace(0,1,h)',ones(h,1)];[ones(h,1),linspace(1,0,h)',linspace(1,0,h)']];
end

% ── NUEVA: waveform CSD allch por sesion ──────────────────────────────────
% Formato: filas = canales, columnas = time points
% Metadatos: Group, Session, SinkType, Stat (mean/sem), Channel
function tOut = buildCSDAllChWaveformsBySess(csdCell, sesNames, grpLabel, stype)
    grpA=[]; sesA=[]; stA=[]; statA=[]; chA=[]; numA=[]; colN=[];
    for ii = 1:numel(csdCell)
        if isempty(csdCell{ii}), continue; end
        c   = csdCell{ii};
        sr  = if_else(isfield(c,'samplingRate') && ~isempty(c.samplingRate), c.samplingRate, 1250);
        ts  = c.timestamps(:)' / sr * 1000;   % ms
        nCh = size(c.mean, 2);
        if isempty(colN)
            colN = string(arrayfun(@(t) sprintf('t_%.2fms', t), ts, 'UniformOutput', false));
        end
        for statStr = {"mean","sem"}
            mat = if_else(strcmp(statStr{1},"mean"), c.mean, c.sem);  % [T x nCh]
            grpA  = [grpA;  repmat(string(grpLabel),   nCh, 1)];      %#ok<AGROW>
            sesA  = [sesA;  repmat(string(sesNames{ii}),nCh, 1)];     %#ok<AGROW>
            stA   = [stA;   repmat(string(stype),       nCh, 1)];     %#ok<AGROW>
            statA = [statA; repmat(string(statStr{1}),  nCh, 1)];     %#ok<AGROW>
            chA   = [chA;   (1:nCh)'];                                 %#ok<AGROW>
            numA  = [numA;  mat'];                                     %#ok<AGROW>  [nCh x T]
        end
    end
    if isempty(numA), tOut = []; return; end
    tOut = array2table(numA, 'VariableNames', colN);
    tOut.Group    = grpA;
    tOut.Session  = sesA;
    tOut.SinkType = stA;
    tOut.Stat     = statA;
    tOut.Channel  = chA;
    tOut = movevars(tOut, {'Group','Session','SinkType','Stat','Channel'}, 'Before', 1);
end

% ── NUEVA: waveform CSD allch media grupal ────────────────────────────────
function tOut = buildCSDAllChWaveformsGroup(csdCell, grpLabel, stype)
    vm = ~cellfun(@isempty, csdCell);
    if ~any(vm), tOut = []; return; end
    gs  = stackCSDCellsAllCh(csdCell, vm);
    sr  = if_else(isfield(gs,'samplingRate') && ~isempty(gs.samplingRate), gs.samplingRate, 1250);
    ts  = gs.timestamps(:)' / sr * 1000;
    nCh = size(gs.mean, 2);
    nSess = sum(vm);
    colN = string(arrayfun(@(t) sprintf('t_%.2fms', t), ts, 'UniformOutput', false));
    grpA=[]; stA=[]; statA=[]; chA=[]; nSA=[]; numA=[];
    for statStr = {"mean","sem"}
        mat = if_else(strcmp(statStr{1},"mean"), gs.mean, gs.sem);  % [T x nCh]
        grpA  = [grpA;  repmat(string(grpLabel),        nCh, 1)];   %#ok<AGROW>
        stA   = [stA;   repmat(string(stype),            nCh, 1)];   %#ok<AGROW>
        statA = [statA; repmat(string(statStr{1}),       nCh, 1)];   %#ok<AGROW>
        chA   = [chA;   (1:nCh)'];                                    %#ok<AGROW>
        nSA   = [nSA;   repmat(nSess, nCh, 1)];                      %#ok<AGROW>
        numA  = [numA;  mat'];                                        %#ok<AGROW>
    end
    tOut = array2table(numA, 'VariableNames', colN);
    tOut.Group     = grpA;
    tOut.SinkType  = stA;
    tOut.Stat      = statA;
    tOut.Channel   = chA;
    tOut.nSessions = nSA;
    tOut = movevars(tOut, {'Group','SinkType','Stat','Channel','nSessions'}, 'Before', 1);
end

function saveCSDFig(csdAll, csdRad, csdSlm, sesName, grpLabel, figFolder)
    structs={csdAll,csdRad,csdSlm}; labels={'All','Rad^{sink}','LM^{sink}'};
    nP=sum(~cellfun(@isempty,structs)); if nP==0, return; end
    allD=[]; for k=1:3; if ~isempty(structs{k}); tmp=extractLayerData(structs{k}); allD=[allD;tmp(:)]; end; end %#ok<AGROW>
    cl=max(abs(allD))*[-1 1]; if diff(cl)==0, cl=[-1 1]; end
    fig=figure('Units','normalized','Position',[0.05 0.1 min(0.32*nP,0.92) 0.62],'Color','w');
    sgtitle(sprintf('%s — %s',grpLabel,strrep(sesName,'_','\_')),'FontSize',11,'FontWeight','bold');
    pi=1;
    for k=1:3
        if isempty(structs{k}), continue; end
        ax=subplot(1,nP,pi); pi=pi+1;
        nStr=if_else(k>1,sprintf(' (n=%d)',structs{k}.nRip),'');
        plotCSDPanel(ax,structs{k},[labels{k} nStr],cl,pi-1==nP);
    end
    fname=fullfile(figFolder,sprintf('csd_%s_%s',grpLabel,sesName));
    saveas(fig,[fname '.png']); saveas(fig,[fname '.fig']); close(fig);
end

function saveGroupCSDFig(csdAll, csdRad, csdSlm, grpLabel, ~, figFolder)
    vA=~cellfun(@isempty,csdAll); vR=~cellfun(@isempty,csdRad); vS=~cellfun(@isempty,csdSlm);
    if ~any(vA), return; end
    mA=stackCSDCells(csdAll,vA); mR=[]; mS=[];
    if any(vR), mR=stackCSDCells(csdRad,vR); end
    if any(vS), mS=stackCSDCells(csdSlm,vS); end
    structs={mA,mR,mS}; labels={'All avg','Rad^{sink} avg','LM^{sink} avg'};
    nP=sum(~cellfun(@isempty,structs));
    allD=[]; for k=1:3; if ~isempty(structs{k}); tmp=extractLayerData(structs{k}); allD=[allD;tmp(:)]; end; end %#ok<AGROW>
    cl=max(abs(allD))*[-1 1]; if diff(cl)==0, cl=[-1 1]; end
    fig=figure('Units','normalized','Position',[0.05 0.1 min(0.32*nP,0.92) 0.62],'Color','w');
    sgtitle(sprintf('GROUP: %s — %d sesiones',grpLabel,sum(vA)),'FontSize',12,'FontWeight','bold');
    pi=1;
    for k=1:3
        if isempty(structs{k}), continue; end
        ax=subplot(1,nP,pi); pi=pi+1;
        plotCSDPanel(ax,structs{k},[labels{k} sprintf(' (n=%d)',structs{k}.nRip)],cl,pi-1==nP);
    end
    fname=fullfile(figFolder,sprintf('csd_group_%s',grpLabel));
    saveas(fig,[fname '.png']); savefig(fig,[fname '.fig']); close(fig);
    fprintf('  Figura CSD grupal: csd_group_%s.png\n',grpLabel);
end

function saveCSDFigAllCh(csdAll, csdRad, csdSlm, sesName, grpLabel, figFolder, stW, spW)
    structs={csdAll,csdRad,csdSlm}; labels={'All','Rad^{sink}','LM^{sink}'};
    nP=sum(~cellfun(@isempty,structs)); if nP==0, return; end
    allD=[]; for k=1:3; if ~isempty(structs{k}); allD=[allD;structs{k}.mean(:)]; end; end %#ok<AGROW>
    cl=max(abs(allD))*[-1 1]; if diff(cl)==0, cl=[-1 1]; end
    fig=figure('Units','normalized','Position',[0.05 0.1 min(0.32*nP,0.92) 0.72],'Color','w');
    sgtitle(sprintf('%s — %s [all ch]',grpLabel,strrep(sesName,'_','\_')),'FontSize',11,'FontWeight','bold');
    pi=1;
    for k=1:3
        if isempty(structs{k}), continue; end
        ax=subplot(1,nP,pi); pi=pi+1;
        nStr=if_else(k>1,sprintf(' (n=%d)',structs{k}.nRip),'');
        plotCSDPanelAllCh(ax,structs{k},[labels{k} nStr],cl,pi-1==nP,stW,spW);
    end
    fname=fullfile(figFolder,sprintf('csd_allch_%s_%s',grpLabel,sesName));
    saveas(fig,[fname '.png']); saveas(fig,[fname '.fig']); close(fig);
end

function saveGroupCSDFigAllCh(csdAll, csdRad, csdSlm, grpLabel, ~, figFolder, stW, spW)
    vA=~cellfun(@isempty,csdAll); vR=~cellfun(@isempty,csdRad); vS=~cellfun(@isempty,csdSlm);
    if ~any(vA), return; end
    mA=stackCSDCellsAllCh(csdAll,vA); mR=[]; mS=[];
    if any(vR), mR=stackCSDCellsAllCh(csdRad,vR); end
    if any(vS), mS=stackCSDCellsAllCh(csdSlm,vS); end
    structs={mA,mR,mS}; labels={'All avg','Rad^{sink} avg','LM^{sink} avg'};
    nP=sum(~cellfun(@isempty,structs));
    allD=[]; for k=1:3; if ~isempty(structs{k}); allD=[allD;structs{k}.mean(:)]; end; end %#ok<AGROW>
    cl=max(abs(allD))*[-1 1]; if diff(cl)==0, cl=[-1 1]; end
    fig=figure('Units','normalized','Position',[0.05 0.1 min(0.32*nP,0.92) 0.72],'Color','w');
    sgtitle(sprintf('GROUP: %s — %d sesiones [all ch]',grpLabel,sum(vA)),'FontSize',12,'FontWeight','bold');
    pi=1;
    for k=1:3
        if isempty(structs{k}), continue; end
        ax=subplot(1,nP,pi); pi=pi+1;
        plotCSDPanelAllCh(ax,structs{k},[labels{k} sprintf(' (n=%d)',structs{k}.nRip)],cl,pi-1==nP,stW,spW);
    end
    fname=fullfile(figFolder,sprintf('csd_allch_group_%s',grpLabel));
    saveas(fig,[fname '.png']); savefig(fig,[fname '.fig']); close(fig);
    fprintf('  Figura CSD allch grupal: csd_allch_group_%s.png\n',grpLabel);
end

function plotCSDPanel(ax, csdStruct, titleStr, clim, showCbar)
    if isempty(csdStruct), return; end
    data=extractLayerData(csdStruct); ts=csdStruct.timestamps(:); nLay=size(data,2);
    if isfield(csdStruct,'layerInfo')&&csdStruct.layerInfo.nLayers==nLay
        yL=csdStruct.layerInfo.labels;
    else; yL=arrayfun(@(i)sprintf('Ch%d',i),1:nLay,'UniformOutput',false); end
    sr=if_else(isfield(csdStruct,'samplingRate')&&~isempty(csdStruct.samplingRate),csdStruct.samplingRate,1250);
    tsMs=ts/sr*1000;
    n=256; cmap=[linspace(0,1,n/2)',linspace(0,1,n/2)',ones(n/2,1);ones(n/2,1),linspace(1,0,n/2)',linspace(1,0,n/2)'];
    imagesc(ax,tsMs,1:nLay,data'); set(ax,'YDir','reverse'); colormap(ax,flipud(cmap));
    if ~isempty(clim)&&diff(clim)>0, set(ax,'CLim',clim); end
    if showCbar, cb=colorbar(ax,'Location','eastoutside'); cb.Label.String='CSD'; cb.FontSize=8;
        lo=cb.Limits(1); hi=cb.Limits(2);
        if lo<0&&hi>0, cb.Ticks=[lo,0,hi]; cb.TickLabels={'Sink','0','Source'};
        else; cb.Ticks=[lo,hi]; cb.TickLabels={'Min','Max'}; end; end
    hold(ax,'on');
    for li=1:nLay; tr=data(:,li); mx=max(abs(tr)); if mx>0, tr=tr/mx; end
        plot(ax,tsMs,li+tr*0.42,'k-','LineWidth',0.5); end
    xline(ax,0,'k--','LineWidth',1.2,'Alpha',0.5); hold(ax,'off');
    xlim(ax,[tsMs(1) tsMs(end)]); ylim(ax,[0.5 nLay+0.5]);
    xTk=unique([-100,-50,0,50,100]); xTk=xTk(xTk>=tsMs(1)&xTk<=tsMs(end));
    set(ax,'XTick',xTk,'YTick',1:nLay,'YTickLabel',yL,'FontSize',9,'Box','off','TickDir','out');
    xlabel(ax,'Time from ripple (ms)','FontSize',9); ylabel(ax,'Layer','FontSize',9);
    title(ax,titleStr,'Interpreter','tex','FontSize',10,'FontWeight','bold');
end

function plotCSDPanelAllCh(ax, csdStruct, titleStr, clim, showCbar, stW, spW)
    if isempty(csdStruct), return; end
    data=csdStruct.mean; ts=csdStruct.timestamps(:); nCh=size(data,2);
    sr=if_else(isfield(csdStruct,'samplingRate')&&~isempty(csdStruct.samplingRate),csdStruct.samplingRate,1250);
    tsMs=ts/sr*1000;
    sd=smoothCSDdata(data,stW,spW);
    n=256; h=n/2;
    cmap=[[linspace(0,1,h)',linspace(0,1,h)',ones(h,1)];[ones(h,1),linspace(1,0,h)',linspace(1,0,h)']];
    imagesc(ax,tsMs,1:nCh,sd'); set(ax,'YDir','reverse'); colormap(ax,cmap);
    if ~isempty(clim)&&diff(clim)>0, set(ax,'CLim',clim); end
    if showCbar, cb=colorbar(ax,'Location','eastoutside'); cb.Label.String='CSD (a.u.)'; cb.FontSize=8;
        lo=cb.Limits(1); hi=cb.Limits(2);
        if lo<0&&hi>0, cb.Ticks=[lo,0,hi]; cb.TickLabels={'Sink','0','Source'};
        else; cb.Ticks=[lo,hi]; cb.TickLabels={'Min','Max'}; end; end
    hold(ax,'on');
    for ch=1:nCh; tr=sd(:,ch); mx=max(abs(tr)); if mx>0, tr=tr/mx; end
        plot(ax,tsMs,ch+tr*0.45,'k-','LineWidth',0.35); end
    xline(ax,0,'w-','LineWidth',1.5,'Alpha',0.8);
    if isfield(csdStruct,'layerInfo')&&~isempty(csdStruct.layerInfo)&&csdStruct.layerInfo.nLayers>0
        info=csdStruct.layerInfo; xL=tsMs(1)+0.03*(tsMs(end)-tsMs(1));
        for li=1:info.nLayers; ch=info.channels(li);
            if ch>=1&&ch<=nCh
                yline(ax,ch-0.5,'w--','LineWidth',0.9,'Alpha',0.75);
                text(ax,xL,ch,info.labels{li},'Color','w','FontSize',8,'FontWeight','bold','VerticalAlignment','middle');
            end; end
    end
    hold(ax,'off');
    xlim(ax,[tsMs(1) tsMs(end)]); ylim(ax,[0.5 nCh+0.5]);
    xTk=unique([-100,-50,0,50,100]); xTk=xTk(xTk>=tsMs(1)&xTk<=tsMs(end));
    set(ax,'XTick',xTk,'Box','off','TickDir','out','FontSize',9);
    if isfield(csdStruct,'layerInfo')&&~isempty(csdStruct.layerInfo)&&csdStruct.layerInfo.nLayers>0
        info=csdStruct.layerInfo; vch=info.channels(info.channels>=1&info.channels<=nCh);
        set(ax,'YTick',vch,'YTickLabel',info.labels(info.channels>=1&info.channels<=nCh));
    else; set(ax,'YTick',[]); end
    xlabel(ax,'Time from ripple (ms)','FontSize',9); ylabel(ax,'Channel','FontSize',9);
    title(ax,titleStr,'Interpreter','tex','FontSize',10,'FontWeight','bold');
end

function saveWavFig(wavAll, wavRad, wavSlm, sesName, grpLabel, figFolder)
    structs={wavAll,wavRad,wavSlm}; labels={'All','Rad^{sink}','LM^{sink}'};
    nP=sum(~cellfun(@isempty,structs)); if nP==0, return; end
    allPow=[]; for k=1:3; if ~isempty(structs{k}); allPow=[allPow;structs{k}.meanPow(:)]; end; end %#ok<AGROW>
    cl=[0,prctile(allPow,99)];
    fig=figure('Units','normalized','Position',[0.05 0.1 min(0.35*nP,0.95) 0.55],'Color','w');
    sgtitle(sprintf('%s — %s [Wavelet pyrCh]',grpLabel,strrep(sesName,'_','\_')),'FontSize',11,'FontWeight','bold');
    pi=1;
    for k=1:3
        if isempty(structs{k}), continue; end
        ax=subplot(1,nP,pi); pi=pi+1;
        plotWavPanel(ax,structs{k},[labels{k} sprintf(' (n=%d)',structs{k}.nRip)],cl,pi-1==nP);
    end
    fname=fullfile(figFolder,sprintf('wav_%s_%s',grpLabel,sesName));
    saveas(fig,[fname '.png']); saveas(fig,[fname '.fig']); close(fig);
end

function saveGroupWavFig(wavAll, wavRad, wavSlm, grpLabel, ~, figFolder)
    vA=~cellfun(@isempty,wavAll); vR=~cellfun(@isempty,wavRad); vS=~cellfun(@isempty,wavSlm);
    if ~any(vA), return; end
    mA=stackWavCells(wavAll,vA); mR=[]; mS=[];
    if any(vR), mR=stackWavCells(wavRad,vR); end
    if any(vS), mS=stackWavCells(wavSlm,vS); end
    structs={mA,mR,mS}; labels={'All avg','Rad^{sink} avg','LM^{sink} avg'};
    nP=sum(~cellfun(@isempty,structs));
    allPow=[]; for k=1:3; if ~isempty(structs{k}); allPow=[allPow;structs{k}.meanPow(:)]; end; end %#ok<AGROW>
    cl=[0,prctile(allPow,99)];
    fig=figure('Units','normalized','Position',[0.05 0.1 min(0.35*nP,0.95) 0.55],'Color','w');
    sgtitle(sprintf('GROUP: %s — %d sesiones [Wavelet]',grpLabel,sum(vA)),'FontSize',12,'FontWeight','bold');
    pi=1;
    for k=1:3
        if isempty(structs{k}), continue; end
        ax=subplot(1,nP,pi); pi=pi+1;
        plotWavPanel(ax,structs{k},[labels{k} sprintf(' (n=%d)',structs{k}.nRip)],cl,pi-1==nP);
    end
    fname=fullfile(figFolder,sprintf('wav_group_%s',grpLabel));
    saveas(fig,[fname '.png']); savefig(fig,[fname '.fig']); close(fig);
    fprintf('  Figura wavelet grupal: wav_group_%s.png\n',grpLabel);
end

function plotWavPanel(ax, wavStruct, titleStr, climWav, showCbar)
    if isempty(wavStruct), return; end
    pow=wavStruct.meanPow; tsMs=wavStruct.tsMs(:)'; freq=wavStruct.freqs(:);
    wav=wavStruct.meanWaveform(:)';
    imagesc(ax,tsMs,freq,pow); set(ax,'YDir','normal'); colormap(ax,jet(256));
    if ~isempty(climWav)&&diff(climWav)>0, set(ax,'CLim',climWav); end
    set(ax,'Color',[0 0 0]);
    if showCbar, cb=colorbar(ax,'Location','eastoutside'); cb.Label.String='Power (a.u.)'; cb.FontSize=8; end
    hold(ax,'on');
    fMin=min(freq); fMax=max(freq); fR=fMax-fMin;
    wavN=wav-mean(wav,'omitnan'); mx=max(abs(wavN));
    if mx>0, wavN=wavN/mx*(fR*0.38); end
    plot(ax,tsMs,(fMin+fMax)/2+wavN,'w-','LineWidth',1.5);
    xline(ax,0,'w--','LineWidth',1.0,'Alpha',0.7); hold(ax,'off');
    xlim(ax,[tsMs(1) tsMs(end)]); ylim(ax,[fMin-5 fMax+5]);
    xTk=unique([-100,-50,0,50,100]); xTk=xTk(xTk>=tsMs(1)&xTk<=tsMs(end));
    set(ax,'XTick',xTk,'Box','off','TickDir','out','FontSize',9,'XColor','k','YColor','k');
    xlabel(ax,'Time from ripple (ms)','FontSize',9,'Color','k');
    ylabel(ax,'Freq (Hz)','FontSize',9,'Color','k');
    title(ax,titleStr,'Interpreter','tex','FontSize',10,'FontWeight','bold','Color','k');
end

function tOut = buildCSDExcelBySess(csdCell, sesNames, grpLabel, stype)
    grpA=[]; sesA=[]; stA=[]; layA=[]; numA=[]; colN=[];
    for ii=1:numel(csdCell)
        if isempty(csdCell{ii}), continue; end
        ld=extractLayerData(csdCell{ii}); ts=csdCell{ii}.timestamps(:)'/1250*1000; nL=size(ld,2);
        if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false)); end
        ln=if_else(isfield(csdCell{ii},'layerInfo')&&csdCell{ii}.layerInfo.nLayers==nL,...
            string(csdCell{ii}.layerInfo.labels),string(arrayfun(@(i)sprintf('L%d',i),1:nL,'UniformOutput',false)));
        grpA=[grpA;repmat(string(grpLabel),nL,1)]; sesA=[sesA;repmat(string(sesNames{ii}),nL,1)]; %#ok<AGROW>
        stA=[stA;repmat(string(stype),nL,1)]; layA=[layA;ln(:)]; numA=[numA;ld']; %#ok<AGROW>
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.Layer=layA;
    tOut=movevars(tOut,{'Group','Session','SinkType','Layer'},'Before',1);
end

function tOut = buildCSDExcelGroup(csdCell, grpLabel, stype)
    vm=~cellfun(@isempty,csdCell); if ~any(vm), tOut=[]; return; end
    gs=stackCSDCells(csdCell,vm); lM=gs.mean; lS=gs.sem;
    ts=gs.timestamps(:)'/1250*1000; nL=size(lM,2); nSess=sum(vm);
    colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false));
    ln=if_else(isfield(gs,'layerInfo')&&gs.layerInfo.nLayers==nL,...
        string(gs.layerInfo.labels),string(arrayfun(@(i)sprintf('L%d',i),1:nL,'UniformOutput',false)));
    grpA=repmat(string(grpLabel),nL*2,1); stA=repmat(string(stype),nL*2,1);
    layA=[ln(:);ln(:)]; statA=[repmat("mean",nL,1);repmat("sem",nL,1)];
    tOut=array2table([lM';lS'],'VariableNames',colN);
    tOut.Group=grpA; tOut.SinkType=stA; tOut.Layer=layA; tOut.Stat=statA; tOut.nSessions=repmat(nSess,nL*2,1);
    tOut=movevars(tOut,{'Group','SinkType','Layer','Stat','nSessions'},'Before',1);
end

function tOut = buildCSDExcelAllGroups(outFolder, groupLabels, stype)
    grpA=[]; sesA=[]; stA=[]; layA=[]; statA=[]; nSA=[]; numA=[]; colN=[];
    for g=1:numel(groupLabels)
        mf=fullfile(outFolder,sprintf('csd_ripples_%s.mat',groupLabels{g}));
        if ~isfile(mf), continue; end
        load(mf,'csdGrpData'); cc=csdGrpData.(stype); sn=csdGrpData.sessions;
        for ii=1:numel(cc)
            if isempty(cc{ii}), continue; end
            ld=extractLayerData(cc{ii}); ts=cc{ii}.timestamps(:)'/1250*1000; nL=size(ld,2);
            if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false)); end
            ln=if_else(isfield(cc{ii},'layerInfo')&&cc{ii}.layerInfo.nLayers==nL,...
                string(cc{ii}.layerInfo.labels),string(arrayfun(@(i)sprintf('L%d',i),1:nL,'UniformOutput',false)));
            grpA=[grpA;repmat(string(groupLabels{g}),nL,1)]; sesA=[sesA;repmat(string(sn{ii}),nL,1)]; %#ok<AGROW>
            stA=[stA;repmat(string(stype),nL,1)]; layA=[layA;ln(:)]; %#ok<AGROW>
            statA=[statA;repmat("session_mean",nL,1)]; nSA=[nSA;ones(nL,1)]; numA=[numA;ld']; %#ok<AGROW>
        end
        vm=~cellfun(@isempty,cc); if ~any(vm), continue; end
        gs=stackCSDCells(cc,vm); nL=size(gs.mean,2); nSess=sum(vm);
        ln=if_else(isfield(gs,'layerInfo')&&gs.layerInfo.nLayers==nL,...
            string(gs.layerInfo.labels),string(arrayfun(@(i)sprintf('L%d',i),1:nL,'UniformOutput',false)));
        for ss={"group_mean","group_sem"}
            mat=if_else(strcmp(ss{1},"group_mean"),gs.mean',gs.sem');
            grpA=[grpA;repmat(string(groupLabels{g}),nL,1)]; sesA=[sesA;repmat("ALL",nL,1)]; %#ok<AGROW>
            stA=[stA;repmat(string(stype),nL,1)]; layA=[layA;ln(:)]; %#ok<AGROW>
            statA=[statA;repmat(string(ss{1}),nL,1)]; nSA=[nSA;repmat(nSess,nL,1)]; numA=[numA;mat]; %#ok<AGROW>
        end
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.Layer=layA; tOut.Stat=statA; tOut.nSessions=nSA;
    tOut=movevars(tOut,{'Group','Session','SinkType','Layer','Stat','nSessions'},'Before',1);
end

function tOut = buildCSDHeatmapExcelBySess(csdCell, sesNames, grpLabel, stype, stW, spW)
    grpA=[]; sesA=[]; stA=[]; statA=[]; chA=[]; numA=[]; colN=[];
    for ii=1:numel(csdCell)
        if isempty(csdCell{ii}), continue; end
        sr=if_else(isfield(csdCell{ii},'samplingRate'),csdCell{ii}.samplingRate,1250);
        ts=csdCell{ii}.timestamps(:)'/sr*1000; nCh=size(csdCell{ii}.mean,2);
        if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false)); end
        for statStr={"mean","sem"}
            if strcmp(statStr{1},"mean"), raw=csdCell{ii}.mean; else; raw=csdCell{ii}.sem; end
            sd=smoothCSDdata(raw,stW,spW);
            grpA=[grpA;repmat(string(grpLabel),nCh,1)]; %#ok<AGROW>
            sesA=[sesA;repmat(string(sesNames{ii}),nCh,1)]; %#ok<AGROW>
            stA=[stA;repmat(string(stype),nCh,1)]; %#ok<AGROW>
            statA=[statA;repmat(string(statStr{1}),nCh,1)]; %#ok<AGROW>
            chA=[chA;(1:nCh)']; numA=[numA;sd']; %#ok<AGROW>
        end
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.Stat=statA; tOut.Channel=chA;
    tOut=movevars(tOut,{'Group','Session','SinkType','Stat','Channel'},'Before',1);
end

function tOut = buildCSDHeatmapExcelGroup(csdCell, grpLabel, stype, stW, spW)
    vm=~cellfun(@isempty,csdCell); if ~any(vm), tOut=[]; return; end
    gs=stackCSDCellsAllCh(csdCell,vm);
    sr=if_else(isfield(gs,'samplingRate'),gs.samplingRate,1250);
    ts=gs.timestamps(:)'/sr*1000; nCh=size(gs.mean,2); nSess=sum(vm);
    colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false));
    grpA=[]; stA=[]; statA=[]; chA=[]; nSA=[]; numA=[];
    for statStr={"mean","sem"}
        raw=if_else(strcmp(statStr{1},"mean"),gs.mean,gs.sem);
        sd=smoothCSDdata(raw,stW,spW);
        grpA=[grpA;repmat(string(grpLabel),nCh,1)]; %#ok<AGROW>
        stA=[stA;repmat(string(stype),nCh,1)]; %#ok<AGROW>
        statA=[statA;repmat(string(statStr{1}),nCh,1)]; %#ok<AGROW>
        chA=[chA;(1:nCh)']; nSA=[nSA;repmat(nSess,nCh,1)]; numA=[numA;sd']; %#ok<AGROW>
    end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.SinkType=stA; tOut.Stat=statA; tOut.Channel=chA; tOut.nSessions=nSA;
    tOut=movevars(tOut,{'Group','SinkType','Stat','Channel','nSessions'},'Before',1);
end

function tOut = buildWavExcelBySess(wavCell, sesNames, grpLabel, stype, freqs)
    grpA=[]; sesA=[]; stA=[]; frA=[]; statA=[]; numA=[]; colN=[];
    for ii=1:numel(wavCell)
        if isempty(wavCell{ii}), continue; end
        w=wavCell{ii}; ts=w.tsMs(:)'; nF=numel(freqs);
        if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false)); end
        for ss={"mean","sem"}
            mat=if_else(strcmp(ss{1},"mean"),w.meanPow,w.semPow);
            grpA=[grpA;repmat(string(grpLabel),nF,1)]; sesA=[sesA;repmat(string(sesNames{ii}),nF,1)]; %#ok<AGROW>
            stA=[stA;repmat(string(stype),nF,1)]; frA=[frA;freqs(:)]; %#ok<AGROW>
            statA=[statA;repmat(string(ss{1}),nF,1)]; numA=[numA;mat]; %#ok<AGROW>
        end
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.Freq_Hz=frA; tOut.Stat=statA;
    tOut=movevars(tOut,{'Group','Session','SinkType','Freq_Hz','Stat'},'Before',1);
end

function tOut = buildWavExcelGroup(wavCell, grpLabel, stype, freqs)
    vm=~cellfun(@isempty,wavCell); if ~any(vm), tOut=[]; return; end
    gs=stackWavCells(wavCell,vm); ts=gs.tsMs(:)'; nF=numel(freqs); nSess=sum(vm);
    colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false));
    grpA=[]; stA=[]; frA=[]; statA=[]; nSA=[]; numA=[];
    for ss={"mean","sem"}
        mat=if_else(strcmp(ss{1},"mean"),gs.meanPow,gs.semPow);
        grpA=[grpA;repmat(string(grpLabel),nF,1)]; stA=[stA;repmat(string(stype),nF,1)]; %#ok<AGROW>
        frA=[frA;freqs(:)]; statA=[statA;repmat(string(ss{1}),nF,1)]; %#ok<AGROW>
        nSA=[nSA;repmat(nSess,nF,1)]; numA=[numA;mat]; %#ok<AGROW>
    end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.SinkType=stA; tOut.Freq_Hz=frA; tOut.Stat=statA; tOut.nSessions=nSA;
    tOut=movevars(tOut,{'Group','SinkType','Freq_Hz','Stat','nSessions'},'Before',1);
end

function tOut = buildWaveformsExcel(wavCell, sesNames, grpLabel, stype)
    grpA=[]; sesA=[]; stA=[]; ridA=[]; numA=[]; colN=[];
    for ii=1:numel(wavCell)
        if isempty(wavCell{ii}), continue; end
        w=wavCell{ii}; wf=w.allWaveforms; nR=size(wf,2); ts=w.tsMs(:)';
        if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.3fms',t),ts,'UniformOutput',false)); end
        grpA=[grpA;repmat(string(grpLabel),nR,1)]; sesA=[sesA;repmat(string(sesNames{ii}),nR,1)]; %#ok<AGROW>
        stA=[stA;repmat(string(stype),nR,1)]; ridA=[ridA;(1:nR)']; numA=[numA;wf']; %#ok<AGROW>
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.RippleIdx=ridA;
    tOut=movevars(tOut,{'Group','Session','SinkType','RippleIdx'},'Before',1);
end

function tOut = buildWavExcelAllGroups(outFolder, groupLabels, stype)
    grpA=[]; sesA=[]; stA=[]; frA=[]; statA=[]; nSA=[]; numA=[]; colN=[];
    for g=1:numel(groupLabels)
        mf=fullfile(outFolder,sprintf('wavelet_ripples_%s.mat',groupLabels{g}));
        if ~isfile(mf), continue; end
        load(mf,'wavGrpData'); wc=wavGrpData.(stype); sn=wavGrpData.sessions; freqs=wavGrpData.freqs; nF=numel(freqs);
        for ii=1:numel(wc)
            if isempty(wc{ii}), continue; end
            w=wc{ii}; ts=w.tsMs(:)';
            if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false)); end
            for ss={"mean","sem"}
                mat=if_else(strcmp(ss{1},"mean"),w.meanPow,w.semPow);
                grpA=[grpA;repmat(string(groupLabels{g}),nF,1)]; sesA=[sesA;repmat(string(sn{ii}),nF,1)]; %#ok<AGROW>
                stA=[stA;repmat(string(stype),nF,1)]; frA=[frA;freqs(:)]; %#ok<AGROW>
                statA=[statA;repmat(string(ss{1}),nF,1)]; nSA=[nSA;ones(nF,1)]; numA=[numA;mat]; %#ok<AGROW>
            end
        end
        vm=~cellfun(@isempty,wc); if ~any(vm), continue; end
        gs=stackWavCells(wc,vm); nSess=sum(vm);
        for ss={"group_mean","group_sem"}
            mat=if_else(strcmp(ss{1},"group_mean"),gs.meanPow,gs.semPow);
            grpA=[grpA;repmat(string(groupLabels{g}),nF,1)]; sesA=[sesA;repmat("ALL",nF,1)]; %#ok<AGROW>
            stA=[stA;repmat(string(stype),nF,1)]; frA=[frA;freqs(:)]; %#ok<AGROW>
            statA=[statA;repmat(string(ss{1}),nF,1)]; nSA=[nSA;repmat(nSess,nF,1)]; numA=[numA;mat]; %#ok<AGROW>
        end
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.Freq_Hz=frA; tOut.Stat=statA; tOut.nSessions=nSA;
    tOut=movevars(tOut,{'Group','Session','SinkType','Freq_Hz','Stat','nSessions'},'Before',1);
end

function tOut = buildWavHeatmapExcelBySess(wavCell, sesNames, grpLabel, stype, freqs, stW)
    grpA=[]; sesA=[]; stA=[]; frA=[]; statA=[]; numA=[]; colN=[];
    for ii=1:numel(wavCell)
        if isempty(wavCell{ii}), continue; end
        w=wavCell{ii}; ts=w.tsMs(:)'; nF=numel(freqs);
        if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false)); end
        for ss={"mean","sem"}
            raw=if_else(strcmp(ss{1},"mean"),w.meanPow,w.semPow);
            sm=smoothWavPow(raw,stW);
            grpA=[grpA;repmat(string(grpLabel),nF,1)]; sesA=[sesA;repmat(string(sesNames{ii}),nF,1)]; %#ok<AGROW>
            stA=[stA;repmat(string(stype),nF,1)]; frA=[frA;freqs(:)]; %#ok<AGROW>
            statA=[statA;repmat(string(ss{1}),nF,1)]; numA=[numA;sm]; %#ok<AGROW>
        end
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.Freq_Hz=frA; tOut.Stat=statA;
    tOut=movevars(tOut,{'Group','Session','SinkType','Freq_Hz','Stat'},'Before',1);
end

function tOut = buildWavHeatmapExcelGroup(wavCell, grpLabel, stype, freqs, stW)
    vm=~cellfun(@isempty,wavCell); if ~any(vm), tOut=[]; return; end
    gs=stackWavCells(wavCell,vm); ts=gs.tsMs(:)'; nF=numel(freqs); nSess=sum(vm);
    colN=string(arrayfun(@(t)sprintf('t_%.2fms',t),ts,'UniformOutput',false));
    grpA=[]; stA=[]; frA=[]; statA=[]; nSA=[]; numA=[];
    for ss={"mean","sem"}
        raw=if_else(strcmp(ss{1},"mean"),gs.meanPow,gs.semPow);
        sm=smoothWavPow(raw,stW);
        grpA=[grpA;repmat(string(grpLabel),nF,1)]; stA=[stA;repmat(string(stype),nF,1)]; %#ok<AGROW>
        frA=[frA;freqs(:)]; statA=[statA;repmat(string(ss{1}),nF,1)]; %#ok<AGROW>
        nSA=[nSA;repmat(nSess,nF,1)]; numA=[numA;sm]; %#ok<AGROW>
    end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.SinkType=stA; tOut.Freq_Hz=frA; tOut.Stat=statA; tOut.nSessions=nSA;
    tOut=movevars(tOut,{'Group','SinkType','Freq_Hz','Stat','nSessions'},'Before',1);
end

function tOut = buildWaveformsSmoothedExcel(wavCell, sesNames, grpLabel, stype, stW)
    grpA=[]; sesA=[]; stA=[]; ridA=[]; numA=[]; colN=[];
    for ii=1:numel(wavCell)
        if isempty(wavCell{ii}), continue; end
        w=wavCell{ii}; wf=w.allWaveforms; nR=size(wf,2); ts=w.tsMs(:)';
        if isempty(colN), colN=string(arrayfun(@(t)sprintf('t_%.3fms',t),ts,'UniformOutput',false)); end
        wfSm=nan(size(wf));
        for jj=1:nR; wfSm(:,jj)=smoothWaveform(wf(:,jj),stW); end
        grpA=[grpA;repmat(string(grpLabel),nR,1)]; sesA=[sesA;repmat(string(sesNames{ii}),nR,1)]; %#ok<AGROW>
        stA=[stA;repmat(string(stype),nR,1)]; ridA=[ridA;(1:nR)']; numA=[numA;wfSm']; %#ok<AGROW>
    end
    if isempty(numA), tOut=[]; return; end
    tOut=array2table(numA,'VariableNames',colN);
    tOut.Group=grpA; tOut.Session=sesA; tOut.SinkType=stA; tOut.RippleIdx=ridA;
    tOut=movevars(tOut,{'Group','Session','SinkType','RippleIdx'},'Before',1);
end

function out = if_else(cond, a, b)
    if cond, out=a; else; out=b; end
end