function maxIntensityAllScenes = runSetScale(pathToFolder,cal)
%runSetScale
%
% Usage:
%   maxIntensityAllScenes = runSetScale(pathToFolder,cal)
%
% Description:
%   For each ISET3d scene in the specified folder, get the maximum
%   intensity within the scene. Output the maximum intensity across all of
%   the scenes, to use for scaling all of the images within the gamut of 
%   the monitor.
%
% Inputs:
%   pathToFolder : (char) Path to the input folder with all of the images
%   cal          : (struct) Calibration file for the electrophysiology experimental machine
%
% History:
%   10/27/21  amn  Wrote it.

%% Parse the input
parser = inputParser();
parser.addRequired('pathToFolder',@(x)(ischar(x)));
parser.addRequired('cal',@(x)(isstruct(x)));
parser.parse(pathToFolder,cal);

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