
% CFC = CrossFrequencyCoupling_adapted('plt', true);

%% >>>>>>>>>>>>>>> CROSS FREQUENCY COUPLING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%%  1. LFP paths

mainFolder = uigetdir(); % Change this to your main folder path
folderInfo = dir(mainFolder);

% Filter only directories and exclude '.' and '..'
subfolders = folderInfo([folderInfo.isdir]); % Keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

% Main folders:
lfppaths = {subfolders.name};
% lfppaths = lfppaths(1:end-1);

%% 2. Extract data:

 for ii = 1:length(lfppaths)
        thisfoldername = strcat(mainFolder, "\", lfppaths(ii));
        % Go to folder and load data:
        cd(thisfoldername)
        layersname = strcat(thisfoldername, "\", lfppaths(ii), ".hippocampalLayers.channelinfo.mat");
        load(layersname)
        orch = hippocampalLayers.layers{1, 1}.oriens;
        pyrch = hippocampalLayers.layers{1, 1}.pyramidal;
        radch = hippocampalLayers.layers{1, 1}.radiatum;
        slmch = hippocampalLayers.layers{1, 1}.slm;

        thisCFC_or = CrossFrequencyCoupling_adapted('channel', orch,'saveMat', false,'plt', false);
        thisCFC_pyr = CrossFrequencyCoupling_adapted('channel', pyrch,'saveMat', false,'plt', false);
        thisCFC_rad = CrossFrequencyCoupling_adapted('channel', radch,'saveMat', false,'plt', false);
        thisCFC_slm = CrossFrequencyCoupling_adapted('channel', slmch,'saveMat', false,'plt', false);

        thisID_CFC.or = thisCFC_or; thisID_CFC.pyr = thisCFC_pyr;
        thisID_CFC.rad = thisCFC_rad; thisID_CFC.slm = thisCFC_slm;

        savename = strcat(thisfoldername ,"\", lfppaths(ii), ".thisIDCFC.mat");
        save(savename, "thisID_CFC")
        group_CFC{ii} = thisID_CFC;

 end

SS_CFC = group_CFC;
OS_CFC = group_CFC;
OD_CFC = group_CFC;

%% 

SSidtable = SSdata.idtable; OSidtable = OSdata.idtable; ODidtable = ODdata.idtable;

% 1. Extract modulation index and arrange for stats:

% SS data:

SSmiceid = SSidtable.ID; SSmice_n = size(SSmiceid, 1);
layersname = ["or", "pyr", "rad", "slm"];

for ii = 1:SSmice_n
    tempCFC = SS_CFC{1, ii};
    for jj = 1:length(layersname)
        SS_lowgamma(ii,jj) = tempCFC.(layersname(jj)).lGamma_MI;
        SS_midgamma(ii,jj) = tempCFC.(layersname(jj)).mGamma_MI;
        SS_highgamma(ii,jj) = tempCFC.(layersname(jj)).hGamma_MI;

        SS_lowgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.lGamma_MI;
        SS_midgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.mGamma_MI;
        SS_highgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.hGamma_MI;

        SS_lowgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.lGamma_MI;
        SS_midgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.mGamma_MI;
        SS_highgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.hGamma_MI;

        SS_lowgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.lGamma_MI;
        SS_midgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.mGamma_MI;
        SS_highgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.hGamma_MI;
    end
end

SS_CFCtable_lowg = array2table(SS_lowgamma, 'VariableNames',layersname);
SS_CFCtable_midg = array2table(SS_midgamma, 'VariableNames',layersname);
SS_CFCtable_highg = array2table(SS_highgamma, 'VariableNames',layersname);

SS_CFCtable_lowg_event = array2table(SS_lowgamma_event, 'VariableNames',layersname);
SS_CFCtable_midg_event = array2table(SS_midgamma_event, 'VariableNames',layersname);
SS_CFCtable_highg_event = array2table(SS_highgamma_event, 'VariableNames',layersname);

SS_CFCtable_lowg_nonloc = array2table(SS_lowgamma_nonloc, 'VariableNames',layersname);
SS_CFCtable_midg_nonloc = array2table(SS_midgamma_nonloc, 'VariableNames',layersname);
SS_CFCtable_highg_nonloc = array2table(SS_highgamma_nonloc, 'VariableNames',layersname);

SS_CFCtable_lowg_nonlocnonripple = array2table(SS_lowgamma_nonlocnonripple, 'VariableNames',layersname);
SS_CFCtable_midg_nonlocnonripple  = array2table(SS_midgamma_nonlocnonripple, 'VariableNames',layersname);
SS_CFCtable_highg_nonlocnonripple  = array2table(SS_highgamma_nonlocnonripple, 'VariableNames',layersname);


SS_CFCtable_lowg = [SSidtable SS_CFCtable_lowg];
SS_CFCtable_midg = [SSidtable SS_CFCtable_midg];
SS_CFCtable_highg = [SSidtable SS_CFCtable_highg];

SS_CFCtable_lowg_event = [SSidtable SS_CFCtable_lowg_event];
SS_CFCtable_midg_event = [SSidtable SS_CFCtable_midg_event];
SS_CFCtable_highg_event = [SSidtable SS_CFCtable_highg_event];

SS_CFCtable_lowg_nonloc = [SSidtable SS_CFCtable_lowg_nonloc];
SS_CFCtable_midg_nonloc = [SSidtable SS_CFCtable_midg_nonloc];
SS_CFCtable_highg_nonloc = [SSidtable SS_CFCtable_highg_nonloc];

SS_CFCtable_lowg_nonlocnonripple  = [SSidtable SS_CFCtable_lowg_nonlocnonripple];
SS_CFCtable_midg_nonlocnonripple  = [SSidtable SS_CFCtable_midg_nonlocnonripple];
SS_CFCtable_highg_nonlocnonripple  = [SSidtable SS_CFCtable_highg_nonlocnonripple];


% OS data:.


OSmiceid = OSidtable.ID; OSmice_n = size(OSmiceid, 1);

for ii = 1:OSmice_n
    tempCFC = OS_CFC{1, ii};
    for jj = 1:length(layersname)
        OS_lowgamma(ii,jj) = tempCFC.(layersname(jj)).lGamma_MI;
        OS_midgamma(ii,jj) = tempCFC.(layersname(jj)).mGamma_MI;
        OS_highgamma(ii,jj) = tempCFC.(layersname(jj)).hGamma_MI;

        OS_lowgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.lGamma_MI;
        OS_midgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.mGamma_MI;
        OS_highgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.hGamma_MI;

        OS_lowgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.lGamma_MI;
        OS_midgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.mGamma_MI;
        OS_highgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.hGamma_MI;

        OS_lowgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.lGamma_MI;
        OS_midgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.mGamma_MI;
        OS_highgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.hGamma_MI;
    end
end

OS_CFCtable_lowg = array2table(OS_lowgamma, 'VariableNames',layersname);
OS_CFCtable_midg = array2table(OS_midgamma, 'VariableNames',layersname);
OS_CFCtable_highg = array2table(OS_highgamma, 'VariableNames',layersname);

OS_CFCtable_lowg_event = array2table(OS_lowgamma_event, 'VariableNames',layersname);
OS_CFCtable_midg_event = array2table(OS_midgamma_event, 'VariableNames',layersname);
OS_CFCtable_highg_event = array2table(OS_highgamma_event, 'VariableNames',layersname);

OS_CFCtable_lowg_nonloc = array2table(OS_lowgamma_nonloc, 'VariableNames',layersname);
OS_CFCtable_midg_nonloc = array2table(OS_midgamma_nonloc, 'VariableNames',layersname);
OS_CFCtable_highg_nonloc = array2table(OS_highgamma_nonloc, 'VariableNames',layersname);

OS_CFCtable_lowg_nonlocnonripple = array2table(OS_lowgamma_nonlocnonripple, 'VariableNames',layersname);
OS_CFCtable_midg_nonlocnonripple  = array2table(OS_midgamma_nonlocnonripple, 'VariableNames',layersname);
OS_CFCtable_highg_nonlocnonripple  = array2table(OS_highgamma_nonlocnonripple, 'VariableNames',layersname);


OS_CFCtable_lowg = [OSidtable OS_CFCtable_lowg];
OS_CFCtable_midg = [OSidtable OS_CFCtable_midg];
OS_CFCtable_highg = [OSidtable OS_CFCtable_highg];

OS_CFCtable_lowg_event = [OSidtable OS_CFCtable_lowg_event];
OS_CFCtable_midg_event = [OSidtable OS_CFCtable_midg_event];
OS_CFCtable_highg_event = [OSidtable OS_CFCtable_highg_event];

OS_CFCtable_lowg_nonloc = [OSidtable OS_CFCtable_lowg_nonloc];
OS_CFCtable_midg_nonloc = [OSidtable OS_CFCtable_midg_nonloc];
OS_CFCtable_highg_nonloc = [OSidtable OS_CFCtable_highg_nonloc];

OS_CFCtable_lowg_nonlocnonripple  = [OSidtable OS_CFCtable_lowg_nonlocnonripple];
OS_CFCtable_midg_nonlocnonripple  = [OSidtable OS_CFCtable_midg_nonlocnonripple];
OS_CFCtable_highg_nonlocnonripple  = [OSidtable OS_CFCtable_highg_nonlocnonripple];



% OD data:


ODmiceid = ODidtable.ID; ODmice_n = size(ODmiceid, 1);

for ii = 1:ODmice_n
    tempCFC = OD_CFC{1, ii};
    for jj = 1:length(layersname)
        OD_lowgamma(ii,jj) = tempCFC.(layersname(jj)).lGamma_MI;
        OD_midgamma(ii,jj) = tempCFC.(layersname(jj)).mGamma_MI;
        OD_highgamma(ii,jj) = tempCFC.(layersname(jj)).hGamma_MI;

        OD_lowgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.lGamma_MI;
        OD_midgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.mGamma_MI;
        OD_highgamma_event(ii,jj) = tempCFC.(layersname(jj)).theta.hGamma_MI;


        OD_lowgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.lGamma_MI;
        OD_midgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.mGamma_MI;
        OD_highgamma_nonloc(ii,jj) = tempCFC.(layersname(jj)).nonlocinuv.hGamma_MI;

        OD_lowgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.lGamma_MI;
        OD_midgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.mGamma_MI;
        OD_highgamma_nonlocnonripple(ii,jj) = tempCFC.(layersname(jj)).nonlocnonripple.hGamma_MI;
    end
end

OD_CFCtable_lowg = array2table(OD_lowgamma, 'VariableNames',layersname);
OD_CFCtable_midg = array2table(OD_midgamma, 'VariableNames',layersname);
OD_CFCtable_highg = array2table(OD_highgamma, 'VariableNames',layersname);

OD_CFCtable_lowg_event = array2table(OD_lowgamma_event, 'VariableNames',layersname);
OD_CFCtable_midg_event = array2table(OD_midgamma_event, 'VariableNames',layersname);
OD_CFCtable_highg_event = array2table(OD_highgamma_event, 'VariableNames',layersname);

OD_CFCtable_lowg_nonloc = array2table(OD_lowgamma_nonloc, 'VariableNames',layersname);
OD_CFCtable_midg_nonloc = array2table(OD_midgamma_nonloc, 'VariableNames',layersname);
OD_CFCtable_highg_nonloc = array2table(OD_highgamma_nonloc, 'VariableNames',layersname);

OD_CFCtable_lowg_nonlocnonripple = array2table(OD_lowgamma_nonlocnonripple, 'VariableNames',layersname);
OD_CFCtable_midg_nonlocnonripple  = array2table(OD_midgamma_nonlocnonripple, 'VariableNames',layersname);
OD_CFCtable_highg_nonlocnonripple  = array2table(OD_highgamma_nonlocnonripple, 'VariableNames',layersname);


OD_CFCtable_lowg = [ODidtable OD_CFCtable_lowg];
OD_CFCtable_midg = [ODidtable OD_CFCtable_midg];
OD_CFCtable_highg = [ODidtable OD_CFCtable_highg];

OD_CFCtable_lowg_event = [ODidtable OD_CFCtable_lowg_event];
OD_CFCtable_midg_event = [ODidtable OD_CFCtable_midg_event];
OD_CFCtable_highg_event = [ODidtable OD_CFCtable_highg_event];

OD_CFCtable_lowg_nonloc = [ODidtable OD_CFCtable_lowg_nonloc];
OD_CFCtable_midg_nonloc = [ODidtable OD_CFCtable_midg_nonloc];
OD_CFCtable_highg_nonloc = [ODidtable OD_CFCtable_highg_nonloc];

OD_CFCtable_lowg_nonlocnonripple  = [ODidtable OD_CFCtable_lowg_nonlocnonripple];
OD_CFCtable_midg_nonlocnonripple  = [ODidtable OD_CFCtable_midg_nonlocnonripple];
OD_CFCtable_highg_nonlocnonripple  = [ODidtable OD_CFCtable_highg_nonlocnonripple];

% Concatenate all

all_CFCtable_lowg = [SS_CFCtable_lowg; OS_CFCtable_lowg; OD_CFCtable_lowg];
all_CFCtable_midg = [SS_CFCtable_midg; OS_CFCtable_midg; OD_CFCtable_midg];
all_CFCtable_highg = [SS_CFCtable_highg; OS_CFCtable_highg; OD_CFCtable_highg];

all_CFCtable_lowg_event = [SS_CFCtable_lowg_event; OS_CFCtable_lowg_event; OD_CFCtable_lowg_event];
all_CFCtable_midg_event = [SS_CFCtable_midg_event; OS_CFCtable_midg_event; OD_CFCtable_midg_event];
all_CFCtable_highg_event = [SS_CFCtable_highg_event; OS_CFCtable_highg_event; OD_CFCtable_highg_event];

all_CFCtable_lowg_nonloc = [SS_CFCtable_lowg_nonloc; OS_CFCtable_lowg_nonloc; OD_CFCtable_lowg_nonloc];
all_CFCtable_midg_nonloc = [SS_CFCtable_midg_nonloc; OS_CFCtable_midg_nonloc; OD_CFCtable_midg_nonloc];
all_CFCtable_highg_nonloc = [SS_CFCtable_highg_nonloc; OS_CFCtable_highg_nonloc; OD_CFCtable_highg_nonloc];

all_CFCtable_lowg_nonlocnonripple = [SS_CFCtable_lowg_nonlocnonripple; OS_CFCtable_lowg_nonlocnonripple; OD_CFCtable_lowg_nonlocnonripple];
all_CFCtable_midg_nonlocnonripple = [SS_CFCtable_midg_nonlocnonripple; OS_CFCtable_midg_nonlocnonripple; OD_CFCtable_midg_nonlocnonripple];
all_CFCtable_highg_nonlocnonripple = [SS_CFCtable_highg_nonlocnonripple; OS_CFCtable_highg_nonlocnonripple; OD_CFCtable_highg_nonlocnonripple];


savepath = uigetdir();
lowgname = strcat(savepath, "\", "CFC_lowgamma.xlsx");
midgname = strcat(savepath, "\", "CFC_midgamma.xlsx");
highgname = strcat(savepath, "\", "CFC_highgamma.xlsx");
lowgname_event = strcat(savepath, "\", "CFC_lowgamma_event.xlsx");
midgname_event = strcat(savepath, "\", "CFC_midgamma_event.xlsx");
highgname_event = strcat(savepath, "\", "CFC_highgamma_event.xlsx");
lowgname_nonloc = strcat(savepath, "\", "CFC_lowgamma_nonloc.xlsx");
midgname_nonloc = strcat(savepath, "\", "CFC_midgamma_nonloc.xlsx");
highgname_nonloc = strcat(savepath, "\", "CFC_highgamma_nonloc.xlsx");
lowgname_nonlocnonripple = strcat(savepath, "\", "CFC_lowgamma_nonlocnonripple.xlsx");
midgname_nonlocnonripple = strcat(savepath, "\", "CFC_midgamma_nonlocnonripple.xlsx");
highgname_nonlocnonripple = strcat(savepath, "\", "CFC_highgamma_nonlocnonripple.xlsx");

writetable(all_CFCtable_lowg, lowgname)
writetable(all_CFCtable_midg, midgname)
writetable(all_CFCtable_highg, highgname)

writetable(all_CFCtable_lowg_event, lowgname_event)
writetable(all_CFCtable_midg_event, midgname_event)
writetable(all_CFCtable_highg_event, highgname_event)

writetable(all_CFCtable_lowg_nonloc, lowgname_nonloc)
writetable(all_CFCtable_midg_nonloc, midgname_nonloc)
writetable(all_CFCtable_highg_nonloc, highgname_nonloc)

writetable(all_CFCtable_lowg_nonlocnonripple, lowgname_nonlocnonripple)
writetable(all_CFCtable_midg_nonlocnonripple, midgname_nonlocnonripple)
writetable(all_CFCtable_highg_nonlocnonripple, highgname_nonlocnonripple)

