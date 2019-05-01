% this script can be used to find the scale factor for a given LMSStruct

nameOfCalibrationFile = 'VirtualWorldCalibration';
whichCalibration = Inf;
dir = fullfile(getpref('VirtualWorldPsychophysics','calibrationDir'));


%% Load the LMS struct
pathToLMSStruct = '/Users/colorlab/Dropbox (Aguirre-Brainard Lab)/CNST_materials/VirtualWorldPsychophysics/Experiment3/StimuliCondition2_covScaleFactor_0_1/LMSStruct.mat';
temp = load(pathToLMSStruct); 
LMSStruct = temp.LMSStruct; 
clear temp;

%% Load the cal file
cal = LoadCalFile(nameOfCalibrationFile,whichCalibration, dir);

%% Initialize calibration structure for the cones
cal = SetSensorColorSpace(cal, LMSStruct.T_cones, LMSStruct.S); % Fix the last option

%% Find Scale factor

scaleFactor = findScaleFactor(cal, LMSStruct);

%% Experiment 3

% CovScaleFactor = 0_0 ; monitor scalefactor = 6.3094 % This is the same as Condition 1 of Experiment 2
% CovScaleFactor = 0_01 ; monitor scalefactor = 11.7440
% CovScaleFactor = 0_05 ; monitor scalefactor = 9.4623
% CovScaleFactor = 0_1 ; monitor scalefactor = 9.6767
% CovScaleFactor = 0_5 ; monitor scalefactor = 5.4411
% CovScaleFactor = 1 ; monitor scalefactor = 4.9302

% CovScaleFactor = 5 ; monitor scalefactor = 5.0180
% CovScaleFactor = 10 ; monitor scalefactor = 4.9173


%% Experiment 2

% Condition = 1 ; monitor scalefactor = 6.3094
% Condition = 2 ; monitor scalefactor = 4.9302
% Condition = 2a ; monitor scalefactor = 4.6957
% Condition = 3 ; monitor scalefactor = 4.9302
% Condition = 3a ; monitor scalefactor = 4.6957