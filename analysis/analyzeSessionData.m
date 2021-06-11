function sessionAnalysis = analyzeSessionData(varargin)
%analyzeSessionData
%
% Usage:
%   analyzeSessionData();
%
% Description:
%   Analyze psychophysical data from a single data collection session. Save
%   the results (struct 'sessionAnalysis') in the specified output folder.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'subjectName'    : (string)  Name of subject (default: 'AN000')
%   'sessionNumber'  : (scalar)  Number of session (default: 1)
%   'plotFigures'    : (logical) Plot figures if option is on (default: true)
%
% History:
%   06/10/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.addParameter('subjectName', 'AN000', @ischar);
parser.addParameter('sessionNumber', 1, @isscalar);
parser.addParameter('plotFigures', true, @islogical);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
subjectName    = parser.Results.subjectName;
sessionNumber  = parser.Results.sessionNumber;
plotFigures    = parser.Results.plotFigures;

%% Set paths to input and output files
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to data file.
subjectFolder = sprintf('%s%s','subject',subjectName);
dataFile      = sprintf('%s%s_%d.mat','data',subjectName,sessionNumber);
pathToFile    = fullfile(getpref(projectName,'BaseDir'),experimentName, ...
                        'PsychophysicalData',subjectFolder,dataFile);

% Set path to output folder.
pathToOutput = fullfile(getpref(projectName,'BaseDir'), ...
                        experimentName,'PsychophysicalDataAnalysis');
if ~exist(pathToOutput, 'dir')
    mkdir(pathToOutput);
end

% Set path to specific output folder where data will be saved.
pathToOutputFolder = fullfile(pathToOutput,sprintf('%s%s','subject',subjectName));
if ~exist(pathToOutputFolder, 'dir')
    mkdir(pathToOutputFolder);
end
    
% Set path to the file to save.
fileName = sprintf('%s%s_%d.mat','sessionAnalysis',subjectName,sessionNumber);
pathToOutputFile = fullfile(pathToOutputFolder,fileName);
    
%% Load data
%
% Load specified data file.
temp = load(pathToFile,'data'); data = temp.data; clear temp;
    
%% Create vectors with image info
%
% Preallocate vector of target offset amount per image.
targetChange = nan(length(data.imageNames),1);

% Preallocate cell array of target offset units per image.
targetChangeUnits = cell(length(data.imageNames),1);

% Break down image file names to get image info.
for ii = 1:length(data.imageNames)
    name = data.imageNames{ii};
    p1   = strfind(name,'banana');
    p2   = strfind(name,'_');
    p2b  = p2(2);
    
    % Store banana offset amount.
    targetChange(ii) = sscanf(name(p1+6:p2b-1),'%f');
    
    % Store banana offset amount units.
    targetChangeUnits{ii} = sscanf(name(p2b+1:end-4),'%s');
end

%% Create vector of target difference between the two images per trial
%
% Get image indices used per trial (col 1: 1st interval; col 2: 2nd interval).
imageIndicesPerTrial = data.trialOrder;

% Convert each image index to the image's target offset amount.
targetOffsetsPerTrial = targetChange(imageIndicesPerTrial);

% Calculate difference in target offset amount between the two images per trial.
targetDiffPerTrial = diff(targetOffsetsPerTrial,1,2);

%% Calculate performance of observer per target difference
%
% Get observer response per trial.
% Key:
%   1: 2nd target was to the left
%   2: 2nd target was to the right
responsePerTrial = data.selectedResponse;

% Get unique target difference amounts.
targetDiffs = unique(targetDiffPerTrial);

% Calculate performance per target difference.
performancePerTargetDiff = nan(numel(targetDiffs),1);
for ii = 1:numel(targetDiffs)
    
    % Target difference to analyze.
    thisD = targetDiffs(ii);
    
    % Get observer responses for this target difference.
    thisO = responsePerTrial(targetDiffPerTrial==thisD);
    
    % Calculate performance for this target difference.
    % For Key above: calculate proportion trials observer reported that 
    %                the 2nd target was to the right.
    performancePerTargetDiff(ii,1) = sum(thisO==2)/numel(thisO);
end

%% Plot performance of observer per target difference
if plotFigures
    figure; box off;
    plot(targetDiffs,performancePerTargetDiff,'ok');
    hold on;
    title(sprintf('%s:%s%s%d',experimentName,subjectName,'\_',sessionNumber));
    xlabel(sprintf('Comparison offset rightward'));
    ylabel('Proportion selected rightward');
    axis([-Inf Inf 0 1]);
    hold off;
end

%% Save data analysis results
%
% Save data in 'sessionAnalysis' struct.
sessionAnalysis = struct;
sessionAnalysis.experimentName           = experimentName;
sessionAnalysis.subjectName              = subjectName;
sessionAnalysis.sessionNumber            = sessionNumber;
sessionAnalysis.fileNamePerImage         = data.imageNames;
sessionAnalysis.targetOffsetPerImage     = targetChange;
sessionAnalysis.targetUnitsPerImage      = targetChangeUnits;
sessionAnalysis.imageIndicesPerTrial     = imageIndicesPerTrial;
sessionAnalysis.targetOffsetsPerTrial    = targetOffsetsPerTrial;
sessionAnalysis.targetDiffPerTrial       = targetDiffPerTrial;
sessionAnalysis.responsePerTrial         = responsePerTrial;
sessionAnalysis.uniqueTargetDiffs        = targetDiffs;
sessionAnalysis.performancePerTargetDiff = performancePerTargetDiff;

% Save data analysis results file.
save(pathToOutputFile,'sessionAnalysis');
fprintf('\nData was saved in:\n%s\n', pathToOutputFile);

%% End