function runElectrophysRGB(varargin)
%runElectrophysRGB
%
% Usage:
%   runElectrophysRGB()
%
% Description:
%   For each ISET3d scene in the specified folder, convert the ISET3d 
%   image to a calibrated RGB image (not gamma corrected) for the Natural
%   Image Thresholds electrophysiology experiment. Save each RGB image in 
%   the specified output folder.
%
% Optional parameter/value:
%   'experimentName' : (string) Name of experiment folder (default: 'Experiment100')
% 
% History:
%   10/13/21  amn  Wrote it.

%% Parse the input
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set path to input folder
%
% Specify project name.
projectName = 'NaturalImageThresholds';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set to desktop.
pathToFolder = fullfile('/Users','amy','Desktop','scenes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% Set path to input folder.
pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageScenesElectrophys');
%}

%% Set path to electrophysiology setup calibration file, and load file
%
% PittVPixx.mat is a calibration file based on the electrophysiological experiment setup.
% It was created using CalFileXYZNoGamma (in the Calibration folder of BrainardLabToolbox).

% Set path to calibration folder.
calDir = getpref(projectName,'CalDataFolder');

% Set calibration file name.
calFile = 'PittVPixx';

% Load calibration file.
cal = LoadCalFile(calFile,[],calDir);
if (isempty(cal))
    error('Could not find specified calibration file');
end
fprintf('Using calibration done by %s on %s\n', ...
    cal.describe.who,cal.describe.date);

%% Set path to output folder
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set to desktop.
pathToOutput = fullfile('/Users','amy','Desktop','scenesSRGB');

%{
pathToOutput = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageSRGBsElectrophys');

% Create output folder if it doesn't exist.
if ~exist(pathToOutput, 'dir')
    mkdir(pathToOutput);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get names of all scene files in the input folder
%
% List .mat files in the folder.
fileInfo = dir([pathToFolder '/*.mat']);

%% Convert each ISET3d scene to a metameric RGB image and save the RGB image
%
% Convert each .mat file in the input folder.
for ii = 1:length(fileInfo)
    
    % Specify the .mat file.
    fileToLoad = fullfile(pathToFolder,fileInfo(ii).name);
    fileToSave = fullfile(pathToOutput,fileInfo(ii).name);
    
    % If the file hasn't already been converted to an RGB image, convert and save it.
    if ~isfile(fileToSave)
        
        % Load the 'theImage' variable contained in this .mat file.
        temp = load(fileToLoad,'theImage'); theImage = temp.theImage; clear temp;
        
        % Convert this ISET3d scene to a metameric RGB image.
        RGBImage = electrophysRGB(theImage,cal,'showRGB',false);
        
        % Save the RGB image in the output folder.
        theImage = RGBImage;
        save(fileToSave,'theImage');
        fprintf('RGB image saved as %s\n', fileInfo(ii).name);
    end
end

%% End