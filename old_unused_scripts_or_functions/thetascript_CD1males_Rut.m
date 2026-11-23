
%%%%%%% GET THETA MEASURES: USE AFTER HIPPOCOOKBOOK PRE-PROCESSING %%%%%%%

%% 1. Get paths:

mainFolder = uigetdir(); % Change this to your main folder path
folderInfo = dir(mainFolder);

% Filter only directories and exclude '.' and '..'
subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

% Main folders:
lfppaths = {subfolders.name};
% lfppaths = lfppaths(1:end-1);

%% 2. Get theta and gamma measures:

for ii = 1:length(lfppaths)

    % 1. Load theta file:

    thisfoldername = strcat(mainFolder, "\", lfppaths(ii));
    % Go to folder:
    cd(thisfoldername)
    % thetaEpochs{ii} = detectThetaEpochs('saveMat', false, 'force', true, 'bandpass', [4 12], 'delta_bandpass', [1 3.5], 'powerThreshold', 1);
    thetaEpochs{ii} = detectThetaEpochs('saveMat', false, 'force', true, 'bandpass', [6 12], 'delta_bandpass', [1 3.5], 'powerThreshold', 1);
    % thetaEpochs{ii} = detectThetaEpochs('saveMat', false, 'force', true, 'bandpass', [6 9], 'delta_bandpass', [0.5 3.5], 'powerThreshold', 1);

    % 
    % thisfilename = strcat(thisfoldername, "\", lfppaths(ii), ".thetaEpochs.states.mat");
    % if ~exist(thisfilename)
    %     disp(strcat("This File: '", thisfoldername, "' has no theta file"))
    %     continue
    % else
    %     load(thisfilename)
    % end
    % 
    % thetapower = thetaEpochs{ii}.power;
    powts = thetaEpochs{ii}.timestamps.';
    run_intervals = thetaEpochs{ii}.thetaRun.intervals;
    nonrun_intervals = thetaEpochs{ii}.intervals_nonTheta_nonRipples;

    % 2. Load LFP to compute low and high gamma power:

    % cd(thisfoldername)
    layersfile = strcat(thisfoldername, "\", lfppaths(ii), ".hippocampalLayers.channelinfo.mat");    
    load(layersfile)
    slmch = hippocampalLayers.layers{1, 1}.slm;

    lfpT = getLFP(slmch,'noPrompts',true);
    samplingRate = lfpT.samplingRate;

    % 3. Control for events shorter than 1 seg:

    shortidx = (run_intervals(:, 2) - run_intervals(:, 1)) < 1;
    run_intervals = run_intervals(~shortidx, :);
    % 
    % run_thetapower = zeros(size(run_intervals, 1), 1);
    % nonrun_thetapower = zeros(size(nonrun_intervals, 1), 1);
    % run_lgammapower = run_thetapower;
    % nonrun_lgammapower = nonrun_thetapower;
    % run_hgammapower = run_thetapower;
    % nonrun_hgammapower = nonrun_thetapower;


    % 4. Get mean theta, lgamma and hgamma power (run and non-run) 

    % theta 
    % theta = [6 9];
    theta = [4 12];
    [theta_wave,~,~,~,~,~,~,~,~,~]=getWavelet(double(lfpT.data(:,1)),samplingRate,theta(1),theta(2),8,0);
    % theta_wave = rms(abs(theta_wave))';  
    rawthetapowerinuv = median(theta_wave, 1);
    thetapowerindb = pow2db(theta_wave);
    thetapowerindbmedian = median(thetapowerindb, 1);
    zscorethetapower = zscore(thetapowerindb);
    zscorethetamean = mean(zscorethetapower, 1);
    modzthetapower = robust_zscore_matrix(thetapowerindb);
    modztheta = median(modzthetapower, 1);

    % lgamma
    lgamma = [20 60];
    [lgamma_wave,~,~,~,~,~,~,~,~,~]=getWavelet(double(lfpT.data(:,1)),samplingRate,lgamma(1),lgamma(2),8,0);
    % lgamma_wave = rms(abs(lgamma_wave))';
    rawlgammapowerinuv = median(lgamma_wave, 1);
    lgammapowerindb = pow2db(lgamma_wave);
    lgammapowerindbmedian = median(lgammapowerindb,1);
    zscorelgammapower = zscore(lgammapowerindb);
    zscorelgammamean = mean(zscorelgammapower, 1);
    modzlgammapower = robust_zscore_matrix(lgammapowerindb);
    modzlgamma = median(modzlgammapower, 1);

    % hgamma
    hgamma = [60 100]; 
    [hgamma_wave,~,~,~,~,~,~,~,~,~]=getWavelet(double(lfpT.data(:,1)),samplingRate,hgamma(1),hgamma(2),8,0);
    % hgamma_wave = rms(abs(hgamma_wave))';
    % in uv:
    rawhgammapowerinuv = median(hgamma_wave, 1);
    % in db
    hgammapowerindb = pow2db(hgamma_wave);
    hgammapowerindbmedian = median(hgammapowerindb ,1);
    % in zscore mean
    zscorehgammapower = zscore(hgammapowerindb);
    zscorehgammamean = mean(zscorehgammapower, 1);
    % in zscore median
    modzhgammapower = robust_zscore_matrix(hgammapowerindb);
    modzhgamma = median(modzhgammapower, 1);

    % 4.1 Other measures:

    runevents = size(run_intervals, 1);
    totalruntime = sum(run_intervals(:, 2) - run_intervals(:, 1));
    meanruntime = mean(run_intervals(:, 2) - run_intervals(:, 1));
    nonrunevents = size(nonrun_intervals, 1);
    totalnonruntime = sum(nonrun_intervals(:, 2) - nonrun_intervals(:, 1));
    meannonruntime = mean(nonrun_intervals(:, 2) - nonrun_intervals(:, 1));
    if runevents == 0
        run_rawthetapowermedian = nan;run_rawlgammapowermedian = nan;
        run_rawhgammapowermedian = nan;run_thetapowerindb_median = nan;
        run_lgammapowerindb_median = nan; run_hgammapowerindb_median = nan;
        run_thetapowerindb_zscore = nan; run_lgammapowerindb_zscore = nan;
        run_hgammapowerindb_zscore = nan; run_thetapowerindb_zmod = nan;
        run_lgammapowerindb_zmod = nan; run_hgammapowerindb_zmod = nan;
    else 
        for jj = 1:runevents
            temponset = find(powts >= run_intervals(jj, 1), 1, 'first');
            tempoffset = find(powts <= run_intervals(jj, 2), 1, 'last'); % si da problemas cambiar el menor o igual por solo menor
    
            % 1. Raw in uv/s^2:
            run_rawthetapowermedian(jj) = median(rawthetapowerinuv(temponset:tempoffset));
            run_rawlgammapowermedian(jj) = median(rawlgammapowerinuv(temponset:tempoffset));
            run_rawhgammapowermedian(jj) = median(rawhgammapowerinuv(temponset:tempoffset));
    
            % 2. In power (db) >> only from 6 to 9 Hz:
    
            % Mean
            % run_thetapower(jj) = mean(thetapower(temponset:tempoffset));
            % run_lgammapower(jj) = mean(lgamma_wave(temponset:tempoffset));
            % run_hgammapower(jj) = mean(hgamma_wave(temponset:tempoffset));
    
            % Median
            run_thetapowerindb_median(jj) = median(thetapowerindbmedian(temponset:tempoffset));
            run_lgammapowerindb_median(jj) = median(lgammapowerindbmedian(temponset:tempoffset));
            run_hgammapowerindb_median(jj) = median(hgammapowerindbmedian(temponset:tempoffset));
            
            % 3. Zscore (mean)
    
            run_thetapowerindb_zscore(jj) = median(zscorethetamean(temponset:tempoffset));
            run_lgammapowerindb_zscore(jj) = median(zscorelgammamean(temponset:tempoffset));
            run_hgammapowerindb_zscore(jj) = median(zscorehgammamean(temponset:tempoffset));
    
            % 3. Zscore (robust - median + MAD)
    
            run_thetapowerindb_zmod(jj) = median(modztheta(temponset:tempoffset));
            run_lgammapowerindb_zmod(jj) = median(modzlgamma(temponset:tempoffset));
            run_hgammapowerindb_zmod(jj) = median(modzhgamma(temponset:tempoffset));
        end
    end

    if nonrunevents == 0
        
        nonrun_rawthetapowermedian = nan;nonrun_rawlgammapowermedian = nan;
        nonrun_rawhgammapowermedian = nan;nonrun_thetapowerindb_median = nan;
        nonrun_lgammapowerindb_median = nan; nonrun_hgammapowerindb_median = nan;
        nonrun_thetapowerindb_zscore = nan; nonrun_lgammapowerindb_zscore = nan;
        nonrun_hgammapowerindb_zscore = nan; nonrun_thetapowerindb_zmod = nan;
        nonrun_lgammapowerindb_zmod = nan; nonrun_hgammapowerindb_zmod = nan;
    else

        for kk = 1:nonrunevents
            nontemponset = find(powts >= nonrun_intervals(kk, 1), 1, 'first');
            nontempoffset = find(powts <= nonrun_intervals(kk, 2), 1, 'last'); % si da problemas cambiar el menor o igual por solo menor
    
            % 1. Raw in uv/s^2:
            nonrun_rawthetapowermedian(kk) = median(rawthetapowerinuv(nontemponset:nontempoffset));
            nonrun_rawlgammapowermedian(kk) = median(rawlgammapowerinuv(nontemponset:nontempoffset));
            nonrun_rawhgammapowermedian(kk) = median(rawhgammapowerinuv(nontemponset:nontempoffset));
    
            % 2. In power (db) >> only from 6 to 9 Hz:
    
            % Mean
            % run_thetapower(kk) = mean(thetapower(temponset:tempoffset));
            % run_lgammapower(kk) = mean(lgamma_wave(temponset:tempoffset));
            % run_hgammapower(kk) = mean(hgamma_wave(temponset:tempoffset));
    
            % Median
            nonrun_thetapowerindb_median(kk) = median(thetapowerindbmedian(nontemponset:nontempoffset));
            nonrun_lgammapowerindb_median(kk) = median(lgammapowerindbmedian(nontemponset:nontempoffset));
            nonrun_hgammapowerindb_median(kk) = median(hgammapowerindbmedian(nontemponset:nontempoffset));
            
            % 3. Zscore (mean)
    
            nonrun_thetapowerindb_zscore(kk) = median(zscorethetamean(nontemponset:nontempoffset));
            nonrun_lgammapowerindb_zscore(kk) = median(zscorelgammamean(nontemponset:nontempoffset));
            nonrun_hgammapowerindb_zscore(kk) = median(zscorehgammamean(nontemponset:nontempoffset));
    
            % 3. Zscore (robust - median + MAD)
    
            nonrun_thetapowerindb_zmod(kk) = median(modztheta(nontemponset:nontempoffset));
            nonrun_lgammapowerindb_zmod(kk) = median(modzlgamma(nontemponset:nontempoffset));
            nonrun_hgammapowerindb_zmod(kk) = median(modzhgamma(nontemponset:nontempoffset));
    
        end

    end

    % 5. Save new intervals after duration criteria:

    % RUN DATA
    thetaEpochs{ii}.thetaRun.newintervals = run_intervals;
    thetaEpochs{ii}.thetaRun.thetapowerinuv= run_rawthetapowermedian;
    thetaEpochs{ii}.thetaRun.lgammapowerinuv = run_rawlgammapowermedian;
    thetaEpochs{ii}.thetaRun.hgammapowerinuv = run_rawhgammapowermedian;
    thetaEpochs{ii}.thetaRun.thetapowerindb = run_thetapowerindb_median;
    thetaEpochs{ii}.thetaRun.lgammapowerindb  = run_lgammapowerindb_median;
    thetaEpochs{ii}.thetaRun.hgammapowerindb  = run_hgammapowerindb_median;
    thetaEpochs{ii}.thetaRun.thetapowerz = run_thetapowerindb_zscore;
    thetaEpochs{ii}.thetaRun.lgammapowerz  = run_lgammapowerindb_zscore;
    thetaEpochs{ii}.thetaRun.hgammapowerz  = run_hgammapowerindb_zscore;
    thetaEpochs{ii}.thetaRun.thetapowerzmod= run_thetapowerindb_zmod;
    thetaEpochs{ii}.thetaRun.lgammapowerzmod  = run_lgammapowerindb_zmod;
    thetaEpochs{ii}.thetaRun.hgammapowerzmod  = run_hgammapowerindb_zmod;

    % NONRUN DATA
    thetaEpochs{ii}.thetanonrun.thetapower_nonTheta_nonRipples_inuv = nonrun_rawthetapowermedian;
    thetaEpochs{ii}.thetanonrun.lgammapower_nonTheta_nonRipples_inuv = nonrun_rawlgammapowermedian;
    thetaEpochs{ii}.thetanonrun.hgammapower_nonTheta_nonRipples_inuv = nonrun_rawhgammapowermedian;
    thetaEpochs{ii}.thetanonrun.thetapower_nonTheta_nonRipples_indb = nonrun_thetapowerindb_median;
    thetaEpochs{ii}.thetanonrun.lgammapower_nonTheta_nonRipples_indb = nonrun_lgammapowerindb_median;
    thetaEpochs{ii}.thetanonrun.hgammapower_nonTheta_nonRipples_indb = nonrun_hgammapowerindb_median;
    thetaEpochs{ii}.thetanonrun.thetapower_nonTheta_nonRipples_z = nonrun_thetapowerindb_zscore;
    thetaEpochs{ii}.thetanonrun.lgammapower_nonTheta_nonRipples_z = nonrun_lgammapowerindb_zscore;
    thetaEpochs{ii}.thetanonrun.hgammapower_nonTheta_nonRipples_z = nonrun_hgammapowerindb_zscore;
    thetaEpochs{ii}.thetanonrun.thetapower_nonTheta_nonRipples_zmod = nonrun_thetapowerindb_zmod;
    thetaEpochs{ii}.thetanonrun.lgammapower_nonTheta_nonRipples_zmod = nonrun_lgammapowerindb_zmod;
    thetaEpochs{ii}.thetanonrun.hgammapower_nonTheta_nonRipples_zmod = nonrun_hgammapowerindb_zmod;


    % 6. Join all data:

    % allrec(ii, :) = [mean(thetapower) mean(lgamma_wave) mean(hgamma_wave)];
    % run(ii, :) = [runevents totalruntime meanruntime mean(run_thetapower) mean(run_lgammapower) mean(run_hgammapower)];
    % nonrun(ii, :) = [nonrunevents totalnonruntime meannonruntime mean(nonrun_thetapower) mean(nonrun_lgammapower) mean(nonrun_hgammapower)];
    % Median
    % medianallrec(ii, :) = [median(thetapower) median(lgamma_wave) median(hgamma_wave)];
    % medianrun(ii, :) = [runevents totalruntime meanruntime median(run_thetapowermedian) median(run_lgammapowermedian) median(run_hgammapowermedian)];
    % mediannonrun(ii, :) = [nonrunevents totalnonruntime meannonruntime mean(nonrun_thetapowermedian) mean(nonrun_lgammapowermedian) mean(nonrun_hgammapowermedian)];

    zallrec(ii, :) = [median(median(thetapowerindb,1),2) median(median(lgammapowerindb,1),2) median(median(hgammapowerindb,1),2)];
    zrun(ii, :) = [runevents totalruntime meanruntime median(run_rawthetapowermedian)...
        median(run_rawlgammapowermedian) median(run_rawhgammapowermedian)...
        median(run_thetapowerindb_median) median(run_lgammapowerindb_median) median(run_hgammapowerindb_median)...
        mean(run_thetapowerindb_zscore) mean(run_lgammapowerindb_zscore)  mean(run_hgammapowerindb_zscore)...
        median(run_thetapowerindb_zmod) median(run_lgammapowerindb_zmod)  median(run_hgammapowerindb_zmod)];
    znonrun(ii, :) = [nonrunevents totalnonruntime meannonruntime median(nonrun_rawthetapowermedian)...
        median(nonrun_rawlgammapowermedian) median(nonrun_rawhgammapowermedian)...
        median(nonrun_thetapowerindb_median) median(nonrun_lgammapowerindb_median) median(nonrun_hgammapowerindb_median)...
        mean(nonrun_thetapowerindb_zscore) mean(nonrun_lgammapowerindb_zscore)  mean(nonrun_hgammapowerindb_zscore)...
        median(nonrun_thetapowerindb_zmod) median(nonrun_lgammapowerindb_zmod)  median(nonrun_hgammapowerindb_zmod)];
    
    % 7. Save thisthetaEpochs data in a new file:

    thisthetaEpochs = thetaEpochs{ii};
    newfilename = strcat(thisfoldername, "\", lfppaths(ii), ".thetaEpochs6to12.states.mat");
    save(newfilename, 'thisthetaEpochs')
    clearvars -except lfppaths thisthetaEpochs zallrec zrun znonrun mainFolder SStable_allrec SStable_run SStable_nonrun OStable_allrec OStable_run OStable_nonrun ODtable_allrec ODtable_run ODtable_nonrun
    close all
end

%% Classify and join data for each group:

% varnames = ["events (n)", "totaldur (s)", "meandur (s)", "thetapower (db)", "lgammapower (db)", "hgammapower (db)"];
varnames = ["events (n)", "totaldur (s)", "meandur (s)", ...
    "medianTheta (uv)", "medianLgamma (uv)", "medianHgamma (uv)",...
    "medianTheta (db)", "medianLgamma (db)", "medianHgamma (db)",...
    "zTheta (au)", "zLgamma (au)", "zHgamma (au)",...
    "zmodTheta (au)", "zmodLgamma (au)", "zmodHgamma (au)"];
%% SS DATA:

ssid = lfppaths.';
ssexpcond = repmat("SS", length(lfppaths), 1);

% Median
SSdata_allrec = array2table(zallrec, 'VariableNames', ["thetapower (db)", "lgammapower (db)", "hgammapower (db)"]);
SSdata_run = array2table(zrun, 'VariableNames', varnames);
SSdata_nonrun = array2table(znonrun, 'VariableNames', varnames);
SSidtable = table(ssid, ssexpcond, 'VariableNames',["ID", "ExpCond"]);
SStable_allrec = [SSidtable SSdata_allrec];
SStable_run  = [SSidtable SSdata_run];
SStable_nonrun = [SSidtable SSdata_nonrun];

mainFolder = uigetdir();
ssallfile = strcat(mainFolder, "\theta_sstableallrec_4to12median.xlsx");
ssrunfile = strcat(mainFolder, "\theta_sstablerun_4to12median.xlsx");
ssnonrunfile = strcat(mainFolder, "\theta_sstablenonrun_4to12median.xlsx");
writetable(SStable_allrec, ssallfile)
writetable(SStable_run, ssrunfile)
writetable(SStable_nonrun, ssnonrunfile)

% SS NEW TABLE:
% 
% SSdata_run = array2table(zrun, 'VariableNames', varnames);
% SSidtable = table(ssid, ssexpcond, 'VariableNames',["ID", "ExpCond"]);
% SStable_run  = [SSidtable SSdata_run];

% mainFolder = uigetdir();
% ssrunfile = strcat(mainFolder, "\theta_NEWtablerun_6to9median.xlsx");
% writetable(SStable_run, ssrunfile)

%% OS DATA:
osid = lfppaths.';
osexpcond = repmat("OS", length(lfppaths), 1);
OSdata_allrec = array2table(zallrec, 'VariableNames', ["thetapower (db)", "lgammapower (db)", "hgammapower (db)"]);
OSdata_run = array2table(zrun, 'VariableNames', varnames);
OSdata_nonrun = array2table(znonrun, 'VariableNames', varnames);
OSidtable = table(osid, osexpcond, 'VariableNames',["ID", "ExpCond"]);
OStable_allrec = [OSidtable OSdata_allrec];
OStable_run  = [OSidtable OSdata_run];
OStable_nonrun = [OSidtable OSdata_nonrun];

mainFolder = uigetdir();
osallfile = strcat(mainFolder, "\newparamstheta_ostableallrec_4TO12.xlsx");
osrunfile = strcat(mainFolder, "\newparamstheta_ostablerun_4TO12.xlsx");
osnonrunfile = strcat(mainFolder, "\newparamstheta_ostablenonrun_4TO12.xlsx");
writetable(OStable_allrec, osallfile)
writetable(OStable_run, osrunfile)
writetable(OStable_nonrun, osnonrunfile)


% OS TABLE:

% OSdata_run = array2table(zrun, 'VariableNames', varnames);
% OSidtable = table(osid, osexpcond, 'VariableNames',["ID", "ExpCond"]);
% OStable_run  = [OSidtable OSdata_run];
% 
% mainFolder = uigetdir();
% osrunfile = strcat(mainFolder, "\theta_NEWtablerun_6to9median.xlsx");
% writetable(OStable_run, osrunfile)


%% OD DATA:
% lfppaths = lfppaths(1:end-1);
odid = lfppaths.';
odexpcond = repmat("OD", length(lfppaths), 1);
ODdata_allrec = array2table(zallrec, 'VariableNames', ["thetapower (db)", "lgammapower (db)", "hgammapower (db)"]);
ODdata_run = array2table(zrun, 'VariableNames', varnames);
ODdata_nonrun = array2table(znonrun, 'VariableNames', varnames);
ODidtable = table(odid, odexpcond, 'VariableNames',["ID", "ExpCond"]);
ODtable_allrec = [ODidtable ODdata_allrec];
ODtable_run  = [ODidtable ODdata_run];
ODtable_nonrun = [ODidtable ODdata_nonrun];

mainFolder = uigetdir();
odallfile = strcat(mainFolder, "\newparamstheta_odtableallrec_4to12.xlsx");
odrunfile = strcat(mainFolder, "\newparamstheta_odtablerun_4to12.xlsx");
odnonrunfile = strcat(mainFolder, "\newparamstheta_odtablenonrun_4to12.xlsx");
writetable(ODtable_allrec, odallfile)
writetable(ODtable_run, odrunfile)
writetable(ODtable_nonrun, odnonrunfile)

% OD TABLE:

% ODdata_run = array2table(zrun, 'VariableNames', varnames);
% ODidtable = table(odid, odexpcond, 'VariableNames',["ID", "ExpCond"]);
% ODtable_run  = [ODidtable ODdata_run];
% 
% mainFolder = uigetdir();
% odrunfile = strcat(mainFolder, "\theta_NEWtablerun_6to9median.xlsx");
% writetable(Odtable_run, odrunfile)

%% 
% Join all data:

allrec = [SStable_allrec; OStable_allrec; ODtable_allrec];
rundata = [SStable_run; OStable_run; ODtable_run];
nonrundata = [SStable_nonrun; OStable_nonrun; ODtable_nonrun];

datapath = uigetdir();
% allfile = strcat(datapath, "\newparamstheta_allrec.xlsx");
% runfile = strcat(datapath, "\newparamstheta_run.xlsx");
% nonrunfile = strcat(datapath, "\newparamstheta_nonrun.xlsx");
allfile = strcat(datapath, "\newparamstheta_allrec_6to12median.xlsx");
runfile = strcat(datapath, "\newparamstheta_run_6to12median.xlsx");
nonrunfile = strcat(datapath, "\newparamstheta_nonrun_6to12median.xlsx");

% 
% newrunfile = strcat(datapath, "\NEWzsc_run_6to9median.xlsx");
% writetable(rundata, newrunfile)
% 

writetable(allrec, allfile)
writetable(rundata, runfile)
writetable(nonrundata, nonrunfile)
% % 
% SStable_allrec.Properties.VariableNames = ODtable_allrec.Properties.VariableNames;
% SStable_run.Properties.VariableNames = ODtable_run.Properties.VariableNames;
% SStable_nonrun.Properties.VariableNames = ODtable_nonrun.Properties.VariableNames;
% OStable_allrec.Properties.VariableNames = ODtable_allrec.Properties.VariableNames;
% OStable_run.Properties.VariableNames = ODtable_run.Properties.VariableNames;
% OStable_nonrun.Properties.VariableNames = ODtable_nonrun.Properties.VariableNames;
% % 
% % 
% SStable_allrec.ExpCond = string(SStable_allrec.ExpCond);
% SStable_run.ExpCond = string(SStable_run.ExpCond);
% SStable_nonrun.ExpCond = string(SStable_nonrun.ExpCond);
% OStable_allrec.ExpCond = string(OStable_allrec.ExpCond);
% OStable_run.ExpCond = string(OStable_run.ExpCond);
% OStable_nonrun.ExpCond = string(Stable_nonrun.ExpCond);


%% 2-WAY ANOVA FOR PEAK AND POWER TO PEAK DATA:

% Load your data table (with columns: experimental_condition, or, pyr, rad, slm)
savepath = uigetdir();
mylabels = ["SS", "OS", "OD"];
layers = ["or", "pyr", "rad", "slm"];
[stats_results, fig_handle] = analyze_repeated_measures_anova(allrec_peakfq_Theta, mylabels, layers,'allrec_peakfq_Theta', savepath);

%% Locomotion events: compute total duration, mean duration, % of total time spent in locomotion

varnames = ["ID" "expcond" "events(n)" "total duration (s)" "mean duration (s)" "reclength (s)" "timespent in loc (%)"];
% 1. Get paths SS:

mainFolder = uigetdir(); % Change this to your main folder path
folderInfo = dir(mainFolder);

% Filter only directories and exclude '.' and '..'
subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

% Main folders:
SSlfppaths = {subfolders.name};
% lfppaths = lfppaths(1:end-1);

SSevents = zeros(length(SSlfppaths), 1);
SStotaldur = SSevents; SSmeandur = SSevents;
SSreclength = SSevents; SSlocperc = SSevents;

for ii = 1:length(SSlfppaths)
    thisfoldername = strcat(mainFolder, "\", SSlfppaths(ii));
    % Go to folder:
    cd(thisfoldername)
    filename = strcat(thisfoldername, "\", SSlfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(filename)
    temprunint = thisthetaEpochs.thetaRun.newintervals;

    SSevents(ii, 1) = size(temprunint, 1);
    SStotaldur(ii, 1) = sum(temprunint(:, 2)-temprunint(:, 1)); % in sec
    SSmeandur(ii, 1) = mean(temprunint(:, 2)-temprunint(:, 1)); % in sec
    SSreclength(ii, 1) = thisthetaEpochs.thetaRun.timestamps(end); % in sec
    SSlocperc(ii, 1) = SStotaldur(ii, 1)./SSreclength(ii, 1);

end
ssid = SSlfppaths.';
ssexpcond = repmat("SS", length(SSlfppaths), 1);

SStable_locevent = table(ssid, ssexpcond, SSevents, SStotaldur, SSmeandur, SSreclength, SSlocperc, 'VariableNames', varnames); 

% 2. Get paths OS:

varnames = ["ID" "expcond" "events(n)" "total duration (s)" "mean duration (s)" "reclength (s)" "timespent in loc (%)"];

mainFolder = uigetdir(); % Change this to your main folder path
folderInfo = dir(mainFolder);

% Filter only directories and exclude '.' and '..'
subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

% Main folders:
OSlfppaths = {subfolders.name};
% lfppaths = lfppaths(1:end-1);

OSevents = zeros(length(OSlfppaths), 1);
OStotaldur = OSevents; OSmeandur = OSevents;
OSreclength = OSevents; OSlocperc = OSevents;

for ii = 1:length(OSlfppaths)
    thisfoldername = strcat(mainFolder, "\", OSlfppaths(ii));
    % Go to folder:
    cd(thisfoldername)
    filename = strcat(thisfoldername, "\", OSlfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(filename)
    temprunint = thisthetaEpochs.thetaRun.newintervals;

    OSevents(ii, 1) = size(temprunint, 1);
    OStotaldur(ii, 1) = sum(temprunint(:, 2)-temprunint(:, 1)); % in sec
    OSmeandur(ii, 1) = mean(temprunint(:, 2)-temprunint(:, 1)); % in sec
    OSreclength(ii, 1) = thisthetaEpochs.thetaRun.timestamps(end); % in sec
    OSlocperc(ii, 1) = OStotaldur(ii, 1)./OSreclength(ii, 1);

end

osid = OSlfppaths.';
osexpcond = repmat("OS", length(OSlfppaths), 1);

OStable_locevent = table(osid, osexpcond, OSevents, OStotaldur, OSmeandur, OSreclength, OSlocperc, 'VariableNames', varnames); 

% 3. Get path OD: 

varnames = ["ID" "expcond" "events(n)" "total duration (s)" "mean duration (s)" "reclength (s)" "timespent in loc (%)"];

mainFolder = uigetdir(); % Change this to your main folder path
folderInfo = dir(mainFolder);

% Filter only directories and exclude '.' and '..'
subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

% Main folders:
ODlfppaths = {subfolders.name};
% lfppaths = lfppaths(1:end-1);

ODevents = zeros(length(ODlfppaths), 1);
ODtotaldur = ODevents; ODmeandur = ODevents;
ODreclength = ODevents; ODlocperc = ODevents;

for ii = 1:length(ODlfppaths)
    thisfoldername = strcat(mainFolder, "\", ODlfppaths(ii));
    % Go to folder:
    cd(thisfoldername)
    filename = strcat(thisfoldername, "\", ODlfppaths(ii), ".thetaEpochs6to12.states.mat");
    load(filename)
    temprunint = thisthetaEpochs.thetaRun.newintervals;

    ODevents(ii, 1) = size(temprunint, 1);
    ODtotaldur(ii, 1) = sum(temprunint(:, 2)-temprunint(:, 1)); % in sec
    ODmeandur(ii, 1) = mean(temprunint(:, 2)-temprunint(:, 1)); % in sec
    ODreclength(ii, 1) = thisthetaEpochs.thetaRun.timestamps(end); % in sec
    ODlocperc(ii, 1) = ODtotaldur(ii, 1)./ODreclength(ii, 1);

end

odid = ODlfppaths.';
odexpcond = repmat("OD", length(ODlfppaths), 1);

ODtable_locevent = table(odid, odexpcond, ODevents, ODtotaldur, ODmeandur, ODreclength, ODlocperc, 'VariableNames', varnames); 

% Join all tables:

locevent_table = [SStable_locevent; OStable_locevent; ODtable_locevent];

% Plot all data:

plot_cell_measuresII(locevent_table, 'LocomotionData')

% Save corresponding data table:

savepath = uigetdir();
savename = strcat(savepath, "\locevent_table.xlsx");
writetable(locevent_table, savename)

