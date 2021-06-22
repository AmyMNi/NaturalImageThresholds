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

%% Get experiment info
%
% Get the identifies and number of noise levels, and save.
noiseLevels  = unique(data.imageNoiseLevel);
nNoiseLevels = numel(noiseLevels);
data.noiseLevels = noiseLevels;

% Get the identifies and number of conditions, and save (same per noise level).
conditions      = unique(data.imageCondition);
nConditions     = numel(conditions);
data.conditions = conditions;

% Get the identities and number of comparisons (same per condition).
comparisons      = unique(data.imageComparison);
nComparisons     = numel(comparisons);
data.comparisons = comparisons;

%% Analyze performance on each noise level and condition separately
%
% Analyze each noise level separately.
for nn = 1:nNoiseLevels
    
    % Get trial indices for this noise level.
    noiseLevelthis = noiseLevels(nn);
    trialsNoise    = data.trialNoiseLevel==noiseLevelthis;

    % Get the image indices, target offset amount, noise amount, and 
    % observer response per trial.
    imagesN       = data.trialOrder(trialsNoise,:);
    offsetsN      = data.trialOrderComparison(trialsNoise,:);
    noiseAmountsN = data.trialNoiseAmount(trialsNoise,1);
    responsesN    = data.selectedResponse(trialsNoise);
    
    % Analyze performance on each condition separately.
    for ii = 1:nConditions
        
        % Get the center position for this condition.
        centerpos = conditions(ii);
            
        % For the center position of this condition, get the pool of
        % images (of the various noise amounts for this noise level,
        % including the noise amount of 0 from Noise Level 0).
        centerPool = find(data.imageComparison==0 & ...
            data.imageCondition==centerpos & ...
            (data.imageNoiseLevel==noiseLevelthis | data.imageNoiseLevel==0));
        
        % Get the target offset amounts, comparison amount, noise amount,
        % and observer response per trial.
        trialsC       = any(ismember(imagesN,centerPool),2);
        offsetsC      = offsetsN(trialsC,:);
        comparisonsC  = sum(offsetsC,2);
        noiseAmountsC = noiseAmountsN(trialsC); 
        responsesC    = responsesN(trialsC);
        
        % Calculate performance per comparison amount.
        performanceC = nan(nComparisons,1);
        for jj = 1:nComparisons
            comparisonthis  = comparisons(jj);
            offsetsthis     = offsetsC      (comparisonsC==comparisonthis,:);
            responsesthis    = responsesC   (comparisonsC==comparisonthis);
            
            % Calculate proportion observer chose the comparison as rightward.
            if comparisonthis==0
                choseRight = find(responsesthis==2);
            else
                choseRight = [find(offsetsthis(:,1)==0 & responsesthis==2); ...
                              find(offsetsthis(:,2)==0 & responsesthis==1)];
            end
            performanceC(jj) = numel(choseRight)/numel(responsesthis);
        end
        
        % Save noise amounts and performance for this condition (in this noise level).
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevelthis);
        conditionName  = sprintf('%s%d','condition',ii);
        data.noiseAmounts.(noiseLevelName).(conditionName) = noiseAmountsC;
        data.performance.(noiseLevelName).(conditionName) = performanceC;
    end
end

%% Plot performance on each condition separately, for each noise level
if plotFigures
    for nn = 1:nNoiseLevels
        for ii = 1:nConditions
            figure; hold on;
            noiseLevelName  = sprintf('%s%d','noiseLevel',noiseLevels(nn));
            conditionName   = sprintf('%s%d','condition',ii);
            performancethis = data.performance.(noiseLevelName).(conditionName);
            
            % Plot data and psychometric function fit.
            [xOffset,FittedCurve,threshold] = plotPsychometric(noiseLevelName,comparisons,performancethis);
            plot(xOffset,performancethis,'ok','MarkerFace','k');
            plot(xOffset,FittedCurve,'-k','LineWidth',1);
        
            % Plot parameters.
            title({sprintf('%s%s%s%d %s %s: %s %0.1f',experimentName,subjectName,'\_', sessionNumber, ...
                           noiseLevelName,conditionName,'threshold =',threshold),''});
            xlabel(sprintf('Comparison offset rightward (mm)'));
            ylabel('Proportion chose comparison as rightward');
            axis([-Inf Inf 0 1]);
            set(gca,'tickdir','out');
            set(gca,'XTick',xOffset);
            set(gca,'XTickLabel',comparisons);
            box off; hold off;
        end
    end
end

%% Plot performance for all conditions combined, for each noise level
if plotFigures
    for nn = 1:nNoiseLevels
        figure; hold on;
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        conditionName  = sprintf('%s%d','condition',ii);
        
        % Average performance across all conditions.
        performanceAll = nan(nComparisons,nConditions);
        for ii = 1:nConditions
            performanceAll(:,ii) = data.performance.(noiseLevelName).(conditionName);
        end
        performanceAll = mean(performanceAll,2);
        
        % Plot data and psychometric function fit.
        [xOffset,FittedCurve,threshold] = plotPsychometric(noiseLevelName,comparisons,performanceAll);
        plot(xOffset,performanceAll,'ok','MarkerFace','k');
        plot(xOffset,FittedCurve,'-k','LineWidth',1);
            
        % Plot parameters.
        title({sprintf('%s%s%s%d %s: %s %0.1f',experimentName,subjectName,'\_', sessionNumber, ...
                           noiseLevelName,'threshold =',threshold),''});
        xlabel(sprintf('Comparison offset rightward (mm)'));
        ylabel('Proportion chose comparison as rightward');
        axis([-Inf Inf 0 1]);
        set(gca,'tickdir','out');
        set(gca,'XTick',xOffset);
        set(gca,'XTickLabel',comparisons);
        box off; hold off;
    end
end

%% Save data analysis results

if saveData 
    % Save data struct (with analysis result additions).
    save(pathToOutputFile,'data');
    fprintf('\nData was saved in:\n%s\n', pathToOutputFile);
end
end
%% Helper functions

%% Calculate data to plot psychometric function fit

function [xOffset,FittedCurve,threshold] = plotPsychometric(noiseLevelName,comparisons,performance)

% Offset x-axis values so that there are no negative values.
xOffset = comparisons + abs(min(comparisons));

% Calculate Weibull fit.
[estimates,model] = fitWeibull(xOffset,performance);
[~,FittedCurve,p] = model(estimates);
if isnan(p)
	fprintf(2,'Weibull function fit fail for %s\n',noiseLevelName);
end

% Offset threshold according to offset used for x-axis above.
threshold = p(1) - abs(min(comparisons));
end

%% End