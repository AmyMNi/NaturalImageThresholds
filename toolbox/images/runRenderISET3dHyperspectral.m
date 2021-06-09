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
%   'experimentName' : (string) Name of experiment folder (default: 'Experiment000')
% 
% History:
%   06/07/21  amn  Wrote it.

%% Parse the input
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set paths to folders
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to input folder.
pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageScenes');

% Set path to output folder.
pathToOutput = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBs');
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
    
    % Load the 'scene' variable contained in this .mat file.
    temp = load(fileToLoad,'scene'); scene = temp.scene; clear temp;
    
    % Convert this ISET3d scene to a metameric RGB image.
    % Input whether to display the RGBImage and/or the sRGBImage.
    RGBImage = renderISET3dHyperspectral(scene,'showRGB',true,'showSRGB',false);

    % Save the RGB image in the output folder.
    savedFile = fullfile(pathToOutput,fileInfo(ii).name);
    save(savedFile,'RGBImage');
    fprintf('RGB image saved as %s\n', fileInfo(ii).name);
end

%% End