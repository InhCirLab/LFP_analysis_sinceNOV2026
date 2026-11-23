function T = runThetaStage(groupFolders, groupLabels, P, figFolder, runOnly)
%% ========================================================================
%  RUNTHETASTAGE — theta stage (called twice: pass 1 and refresh)
%  runOnly = true  -> run only (generates newintervals), no nonrun/Excel/figs
%  runOnly = false -> full run + nonrun, returns tables and saves figures
%  Recordings are processed serially (one at a time).
%
%  Standalone function (extracted from LFPpipeline.m so sections can be run
%  one at a time). Depends on: processRecording, saveFigSession.
%% ========================================================================

nGroups = numel(groupFolders);
T.allrec = table();  T.run = table();  T.nonrun = table();  T.loc = table();

for g = 1:nGroups

    mainFolder = groupFolders{g};
    expCond    = groupLabels{g};

    entries  = dir(mainFolder);
    entries  = entries([entries.isdir]);
    entries  = entries(~ismember({entries.name}, {'.', '..'}));
    recNames = {entries.name};
    nRec     = numel(recNames);
    assert(nRec > 0, 'No subfolders found in: %s', mainFolder);
    fprintf('\n==== Group "%s": %d recording(s) (serial) ====\n', expCond, nRec);

    outs = cell(nRec, 1);
    for ii = 1:nRec   % SERIAL (was parfor)
        outs{ii} = processRecording(fullfile(mainFolder, recNames{ii}), recNames{ii}, P, runOnly);
    end

    zAllRec = nan(nRec, 3);   zRun = nan(nRec, 15);   zNonRun = nan(nRec, 15);
    locEvents = nan(nRec,1);  locTotalDur = nan(nRec,1);  locMeanDur = nan(nRec,1);
    locRecLen = nan(nRec,1);  locPercent = nan(nRec,1);

    for ii = 1:nRec
        o = outs{ii};
        fprintf('  [%d/%d] %s\n', ii, nRec, recNames{ii});
        for m = 1:numel(o.msg), fprintf('     %s\n', o.msg{m}); end

        if ~runOnly
            zAllRec(ii,:) = o.zAll;      zRun(ii,:) = o.zRun;      zNonRun(ii,:) = o.zNonRun;
            locEvents(ii) = o.nR;        locTotalDur(ii) = o.totDur;  locMeanDur(ii) = o.meanDur;
            locRecLen(ii) = o.recLen;    locPercent(ii)  = o.pct;
            try
                saveFigSession(recNames{ii}, expCond, o.runPow, o.nonrunPow, figFolder);
            catch ME
                warning('saveFigSession failed in %s: %s', recNames{ii}, ME.message);
            end
            close all;
        end
    end

    if runOnly, continue; end

    idCol   = recNames(:);
    condCol = repmat({expCond}, nRec, 1);
    idTab   = table(idCol, condCol, 'VariableNames', ["ID","ExpCond"]);

    T.allrec  = [T.allrec;  idTab  array2table(zAllRec, 'VariableNames', ["thetapower (db)","lgammapower (db)","hgammapower (db)"])];
    T.run     = [T.run;     idTab  array2table(zRun,    'VariableNames', P.powerVarNames)];
    T.nonrun  = [T.nonrun;  idTab  array2table(zNonRun, 'VariableNames', P.powerVarNames)];

    locVars = ["ID","ExpCond","events(n)","totaldur (s)","meandur (s)","reclength (s)","timespent (%)"];
    T.loc = [T.loc; table(idCol, condCol, locEvents, locTotalDur, locMeanDur, locRecLen, locPercent, 'VariableNames', locVars)];

    fprintf('  Group "%s" done.\n', expCond);
end
end
