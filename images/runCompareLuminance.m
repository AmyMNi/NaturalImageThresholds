function runCompareLuminance(varargin)
%runCompareLuminance
%
% Usage:
%   runCompareLuminance()
%
% Description:
%   For each RGB image in the specified folder, calculate the mean
%   luminance of the image, and compare between the images in the different
%   noise levels.
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

%% Set paths to folders and calibration file
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to input folder.
pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBs');

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

%% Get names of all image files in the input folder
%
% List .mat files in the folder.
fileInfo = dir([pathToFolder '/*.mat']);

%% Get the noise level per image
%
% Get the number of images in the folder.
nImages = numel(fileInfo);

% Get the image file names. Images will be called by their index here.
imageNames = {fileInfo(:).name}';

% Create vectors with image info.
imageCondition  = nan(nImages,1);
imageComparison = nan(nImages,1);
imageNoiseLevel = nan(nImages,1);

% Break down image file names to get image info.
for ii = 1:nImages
    name = imageNames{ii};
    p    = strfind(name,'_');
    s1 = name(1:p(1));
    imageCondition(ii)  = str2double(regexp(s1,'[+-]?\d*','Match'));
    s2 = name(p(1):p(2));
    imageComparison(ii) = str2double(regexp(s2,'[+-]?\d*','Match'));
    s3 = name(p(2):p(3));
    imageNoiseLevel(ii) = str2double(regexp(s3,'[+-]?\d*','Match'));
end

% Get the unique image values.
conditions  = unique(imageCondition);
comparisons = unique(imageComparison);
noiseLevels = unique(imageNoiseLevel);

%% Calculate the mean luminance per image
meanLuminance = nan(nImages,1);
for ii = 1:length(fileInfo)
    
    % Specify the .mat file.
    fileToLoad = fullfile(pathToFolder,fileInfo(ii).name);
    
    % Load the 'RGBImage' variable contained in this .mat file.
    temp = load(fileToLoad,'RGBImage'); RGBImage = temp.RGBImage; clear temp;
    
    % Calculate the mean luminance of this image.
    meanLuminance(ii) = mean2(rgb2gray(RGBImage));
end

%% Calculate average luminance at different banana positions for each noise level
%
% Set up vectors of n banana positions for each noise level.
luminanceNoiseLevel0 = nan(numel(conditions)*numel(comparisons),1);
luminanceNoiseLevel1 = nan(numel(conditions)*numel(comparisons),1);
luminanceNoiseLevel2 = nan(numel(conditions)*numel(comparisons),1);

row = 1;
% Calculate average luminance at each banana position.
for ii = 1:numel(conditions)
    conditionThis = conditions(ii);
    for jj = 1:numel(comparisons)
        comparisonThis = comparisons(jj);
        
        % Get indices of images at this banana position for each noise level.
        noiseLevel0This = imageNoiseLevel==0 & imageCondition==conditionThis & imageComparison==comparisonThis;
        noiseLevel1This = imageNoiseLevel==1 & imageCondition==conditionThis & imageComparison==comparisonThis;
        noiseLevel2This = imageNoiseLevel==2 & imageCondition==conditionThis & imageComparison==comparisonThis;
        
        % Calculate average luminance at this banana position for each noise level.
        luminanceNoiseLevel0(row) = mean(meanLuminance(noiseLevel0This));
        luminanceNoiseLevel1(row) = mean(meanLuminance(noiseLevel1This));
        luminanceNoiseLevel2(row) = mean(meanLuminance(noiseLevel2This));
        row = row+1;
    end
end

%% End