function runRenderISET3dHyperspectral(varargin)
%runRenderISET3dHyperspectral
%
% Usage:
%   runRenderISET3dHyperspectral()
%
% Description:
%   For each ISET3d scene in the specified folder, convert the ISET3d
%   scene to a metameric RGB image. Save each RGB image in the specified
%   output folder.
%
% Optional parameter/value:
%   'experimentName' : (string) Name of experiment folder (default: 'Experiment100')
% 
% History:
%   06/07/21  amn  Wrote it.
%   06/14/21  amn  Updated to use calibration file for this experimental
%                  machine, based on t_readUseCalFile.m (written by dhb).

%% Parse the input
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set paths to folders and calibration file
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to input folder.
pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageScenes');

% Set path to calibration folder.
calDir = getpref(projectName,'CalDataFolder');

% Set path to calibration file
% (set for the local experiment machine by the project local hook file).
calFile = getpref(projectName,'CalDataFile');

%% Load calibration file
cal = LoadCalFile(calFile,[],calDir);
if (isempty(cal))
    error('Could not find specified calibration file');
end
fprintf('Using calibration done by %s on %s\n', ...
    cal.describe.who,cal.describe.date);

%% Set path to output folder
if strcmp(calFile,'NaturalImageThresholdsCal_Amy')
    pathToOutput = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBsAmy');
else
    pathToOutput = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBs');
end

% Create output folder if it doesn't exist.
if ~exist(pathToOutput, 'dir')
    mkdir(pathToOutput);
end

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
        
        % Load the 'scene' variable contained in this .mat file.
        temp = load(fileToLoad,'scene'); scene = temp.scene; clear temp;
        
        % Convert this ISET3d scene to a metameric RGB image.
        % Input whether to display the RGBImage and/or the sRGBImage.
        RGBImage = renderISET3dHyperspectral(scene,cal,'showRGB',false,'showSRGB',false);
        
        % Save the RGB image in the output folder.
        save(fileToSave,'RGBImage');
        fprintf('RGB image saved as %s\n', fileInfo(ii).name);
    end
end

%% End