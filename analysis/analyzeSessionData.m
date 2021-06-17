function data = analyzeSessionData(varargin)
%analyzeSessionData
%
% Usage:
%   data = analyzeSessionData('experimentName', 'Experiment000', ...
%                             'subjectName', 'AN', ...
%                             'sessionNumber', 1);
%
% Description:
%   Analyze psychophysical data from a single data collection session. Save
%   the results (struct 'data') in the specified output folder.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'subjectName'    : (string)  Name of subject (default: 'AN')
%   'sessionNumber'  : (scalar)  Number of session (default: 1)
%   'plotFigures'    : (logical) Plot figures if option is on (default: true)
%   'saveData'       : (logical) Save data if option is on (default: true)
%
% History:
%   06/10/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.addParameter('subjectName', 'AN', @ischar);
parser.addParameter('sessionNumber', 1, @isscalar);
parser.addParameter('plotFigures', true, @islogical);
parser.addParameter('saveData', true, @islogical);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
subjectName    = parser.Results.subjectName;
sessionNumber  = parser.Results.sessionNumber;
plotFigures    = parser.Results.plotFigures;
saveData       = parser.Results.saveData;

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

% Save session number in 'data' struct.
data.sessionNumber = sessionNumber;
    
%% Analyze performance on each condition separately
%
% Get the identifies and number of conditions, and save.
conditions      = unique(data.imageCondition);
nConditions     = numel(conditions);
data.conditions = conditions;

% Get the identities and number of comparisons (same per condition).
comparisons      = unique(data.imageComparison);
nComparisons     = numel(comparisons);
data.comparisons = comparisons;

for ii = 1:nConditions
    % Find image index for the image with the center position for this condition.
    centerpos = conditions(ii);
    centerIdx = intersect(find(data.imageCondition==centerpos),find(data.imageComparison==0));
    
    % Get trials in this condition.
    trials = any(data.trialOrder==centerIdx,2);
    
    % Get the comparison amount and observer response per trial.
    offsetsC = data.trialOrderComparison(trials,:);
    comparisonsC = sum(offsetsC,2);
    responsesC = data.selectedResponse(trials);
    
    % Calculate performance per comparison amount.
    performanceC = nan(nComparisons,1);
    for jj = 1:nComparisons
        comparisonThis = comparisons(jj);
        offsetsThis    = offsetsC(comparisonsC==comparisonThis,:);
        responsesThis  = responsesC(comparisonsC==comparisonThis);
        
        % Calculate proportion observer chose the comparison as rightward.
        if comparisonThis==0
            choseRight = find(responsesThis==2);
        else
            choseRight = [intersect(find(offsetsThis(:,1)==0), find(responsesThis==2)); ...
                          intersect(find(offsetsThis(:,2)==0), find(responsesThis==1))];
        end
        performanceC(jj) = numel(choseRight)/numel(responsesThis);
    end
    
    % Save performance for this condition.
    data.performancePerCondition{ii,1} = performanceC;
end

%% Plot performance on each condition separately
if plotFigures
    for ii = 1:nConditions
        figure;
        plot(comparisons,data.performancePerCondition{ii},'ok','MarkerFace','k');
        hold on;
        title({sprintf('%s%s%s%d: %s %d',experimentName,subjectName,'\_',sessionNumber,'condition',ii),''});
        xlabel(sprintf('Comparison offset rightward'));
        ylabel('Proportion chose comparison as rightward');
        axis([-Inf Inf 0 1]);
        set(gca,'tickdir','out');
        box off; hold off;
    end
end

%% Plot performance for all conditions combined
if plotFigures
    % Average performance across all conditions.
    performanceAll = nan(nComparisons,nConditions);
    for ii = 1:nConditions
        performanceAll(:,ii) = data.performancePerCondition{ii};
    end
    performanceAll = mean(performanceAll,2);
    
    % Plot average across conditions.
    figure;
    plot(comparisons,performanceAll,'ok','MarkerFace','k');
    hold on;
    title({sprintf('%s%s%s%d: %s',experimentName,subjectName,'\_',sessionNumber,'all conditions'),''});
    xlabel(sprintf('Comparison offset rightward'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    box off; hold off;
end

%% Save data analysis results

if saveData 
    % Save data struct (with analysis result additions).
    save(pathToOutputFile,'data');
    fprintf('\nData was saved in:\n%s\n', pathToOutputFile);
end

%% End