function maxIntensityAllScenes = runSetScale(varargin)
%runSetScale
%
% Usage:
%   maxIntensityAllScenes = runSetScale()
%
% Description:
%   For each ISET3d scene in the specified folder, get the maximum
%   intensity within the scene. Output the maximum intensity across all of
%   the scenes, to use for scaling all of the images within the gamut of 
%   the monitor.
%
% Optional parameter/value:
%   'experimentName' : (char) Name of experiment folder (default: 'Experiment100')
% 
% History:
%   10/27/21  amn  Wrote it.

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

%% Get names of all scene files in the input folder
%
% List .mat files in the folder.
fileInfo = dir([pathToFolder '/*.mat']);

%% Calculate the maximum intensity within each scene
%
% Calculate maximum intensity for each .mat file in the input folder.
maxIntensity = nan(length(fileInfo),1);
for ii = 1:length(fileInfo)
    
    % Specify the .mat file.
    fileToLoad = fullfile(pathToFolder,fileInfo(ii).name);
    
    % Load the 'theImage' variable contained in this .mat file.
    temp = load(fileToLoad,'theImage'); theImage = temp.theImage; clear temp;
    
    % Calculate the maximum intensity for this image.
    maxIntensity(ii,1) = calcMaxIntensity(theImage,cal);
end

%% Calculate the maximum intensity across all of the scenes
maxIntensityAllScenes = max(maxIntensity);

%% End