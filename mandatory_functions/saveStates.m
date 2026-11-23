function saveStates(statesFile, thisthetaEpochs) %#ok<INUSD>
% SAVESTATES  Save the theta epochs struct to its .states.mat file.
%
% Standalone function (extracted from LFPpipeline.m). The saved variable is
% named 'thisthetaEpochs', which is what Part III and ps1to120_Paria expect.
save(statesFile, 'thisthetaEpochs');
end
