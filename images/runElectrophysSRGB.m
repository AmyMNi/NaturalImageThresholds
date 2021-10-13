function runElectrophysSRGB(varargin)
%runElectrophysSRGB
%
% Usage:
%   runElectrophysSRGB()
%
% Description:
%   For each ISET3d scene in the specified folder, convert the ISET3d 
%   image to an sRGB image (not gamma corrected) for the Natural Image 
%   Thresholds electrophysiology experiment. Save each sRGB image in the 
%   specified output folder.
%
% Optional parameter/value:
%   'experimentName' : (string) Name of experiment folder (default: 'Experiment100')
% 
% History:
%   10/12/21  amn  Wrote it.

%% Parse the input
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set path to input folder
%
% Set to desktop.
pathToFolder = fullfile('/Users','amy','Desktop','scenes');

%% Set path to output folder
%
% Set to desktop.
pathToOutput = fullfile('/Users','amy','Desktop','scenesSRGB');

%% Get names of all scene files in the input folder
%
% List .mat files in the folder.
fileInfo = dir([pathToFolder '/*.mat']);

%% Convert each ISET3d scene to an sRGB image and save the sRGB image
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

        % Convert this ISET3d scene to an sRGB image.
        sRGBImage = electrophysSRGB(theImage,'showSRGB',false);
        
        % Save the sRGB image in the output folder.
        theImage = sRGBImage;
        save(fileToSave,'theImage');
        fprintf('sRGB image saved as %s\n', fileInfo(ii).name);
    end
end

%% End