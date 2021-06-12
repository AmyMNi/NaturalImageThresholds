% t_readUseCalibrationFile

%% Clear
clear; close all;

%% Get preferences (set by project local hook file)
theProject = 'NaturalImageThresholds';
calDir = getpref(theProject,'CalDataFolder');
calFile = getpref(theProject,'CalDataFile');
fprintf('Calibration directory: %s, calibration file: %s\n',calDir,calFile);

%% Load the LMS struct
% LMSstructName = trialStruct.LMSstructName;
% pathToLMSStruct = fullfile(getpref(projectName,'stimulusInputBaseDir'),...
%     directoryName,[LMSstructName '.mat']);
% temp = load(pathToLMSStruct); LMSStruct = temp.LMSStruct; clear temp;

%% Load calibration file
cal = LoadCalFile(calFile,[],calDir);
if (isempty(cal))
    error('Could not find specified calibration file');
end
fprintf('Using calibration done by %s on %s\n', ...
    cal.describe.who,cal.describe.date);

%% Initialize calibration structure for the cones
cal = SetSensorColorSpace(cal, LMSStruct.T_cones, LMSStruct.S); % Fix the last option

%% Find the scale factor
if (scaleFactor == 0)
    scaleFactor = findScaleFactor(cal, LMSStruct);
end

%% Set Gamma Method
cal = SetGammaMethod(cal,0);