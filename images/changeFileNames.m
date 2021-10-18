function changeFileNames(varargin)
%changeFileNames
%
% Usage:
%   changeFileNames()
%
% Description:
%   Change the name of each file in the specified folder, and save the
%   file with the revised name in a new folder.
%
% Optional parameter/value:
%   'experimentName' : (char) Name of experiment folder (default: 'Experiment000')
% 
% History: 
%   06/26/21  amn  Wrote it.

%% Parse the input
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set paths to input and output folders
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to input folder.
pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageScenes');

% Set path to output folder.
pathToOutput = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageScenesRENAMED');
if ~exist(pathToOutput, 'dir')
    mkdir(pathToOutput);
end

%% Get names of all files in the input folder
%
% List .mat files in the folder.
fileInfo = dir([pathToFolder '/*.mat']);

%% Change names of specified files, then save in output folder
%
% Check each .mat file in the input folder.
for ii = 1:length(fileInfo)
    
    % Get the file name.
    fileName = fileInfo(ii).name;
    % Specify the file path.
    fileToLoad = fullfile(pathToFolder,fileName);
    
    % Break down file name.
    p  = strfind(fileName,'_');
    s1 = fileName(1:p(1)-1);
    s2 = fileName(p(1):p(2)-1);
    s4 = fileName(p(3):end);
    
    % Check whether to change file name.
    imageNoiseAmount = str2double(regexp(s4,'[+-]?\d*','Match'));
    if ismember(imageNoiseAmount,[-30 -27 -24 -21 -18 18 21 24 27 30])
        
        % Specify new file name.
        fileName = sprintf('%s%s%s%s',s1,s2,'_noise2',s4);
    end
 
    % Load the 'scene' variable contained in this .mat file.
    temp = load(fileToLoad,'scene'); scene = temp.scene; clear temp;
    
    % Specify new file path.
    fileToSave = fullfile(pathToOutput,fileName);
    
    % Save renamed file in the new output folder.
    save(fileToSave,'scene');
    fprintf('Renamed file saved as %s\n', fileToSave);
end

%% End