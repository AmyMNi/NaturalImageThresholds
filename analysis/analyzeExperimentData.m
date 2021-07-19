function data = analyzeExperimentData(varargin)
%analyzeExperimentData
%
% Usage:
%   data = analyzeExperimentData('experimentName', 'Experiment100', ...
%                                'subjectName', 'test');
%
% Description:
%   Analyze psychophysical data for a single experiment, across all data
%   collection sessions for that experiment. Save the results (struct
%   'dataExperiment') in the specified output folder.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment100')
%   'subjectName'    : (string)  Name of subject (default: 'test')
%   'plotFigures'    : (logical) Plot figures if option is on (default: true)
%   'saveData'       : (logical) Save data if option is on (default: true)
%
% History:
%   07/02/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.addParameter('subjectName', 'test', @ischar);
parser.addParameter('plotFigures', true, @islogical);
parser.addParameter('saveData', true, @islogical);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
subjectName    = parser.Results.subjectName;
plotFigures    = parser.Results.plotFigures;
saveData       = parser.Results.saveData;

%% Set paths to data folder and output file
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to data folder.
subjectFolder = sprintf('%s%s','subject',subjectName);
pathToFolder  = fullfile(getpref(projectName,'BaseDir'),experimentName, ...
                        'PsychophysicalDataAnalysis',subjectFolder);

% Set path to the file to save.
fileName = sprintf('%s%s.mat','experimentAnalysis',subjectName);
pathToOutputFile = fullfile(pathToFolder,fileName);
    
%% Get names of all data files in the data folder
%
% List .mat files in the folder.
fileInfo = dir([pathToFolder '/sessionAnalysis*.mat']);

%% Set up 'dataExperiment' struct
%
% Assumes experiment parameters are the same for all sessions of the experiment.
dataExperiment = struct;

% Load the first .mat file and get all experiment parameters from this first file.
fileToLoad = fullfile(pathToFolder,fileInfo(1).name);
temp = load(fileToLoad,'data'); data = temp.data; clear temp;
fn = fieldnames(data);
fn = fn([1:20 37:41]);
for ii=1:numel(fn)
    dataExperiment.(fn{ii}) = data.(fn{ii});
end

% Get experiment parameters to plot.
noiseLevels = dataExperiment.noiseLevels;
conditions  = dataExperiment.conditions;
comparisons = dataExperiment.comparisons;
nNoiseLevels = numel(noiseLevels);
nConditions  = numel(conditions);
nComparisons = numel(comparisons);

%% Save performance data from each session
for ii = 1:length(fileInfo)
    
    % Specify the .mat file for a session.
    fileToLoad = fullfile(pathToFolder,fileInfo(ii).name);
    
    % Load the data variable contained in this .mat file.
    temp = load(fileToLoad,'data'); data = temp.data; clear temp;
    
    % Save the performance data from this session.
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        NumPos       = nan(nComparisons,nConditions);
        OutOfNum     = nan(nComparisons,nConditions);
        reactionTime = nan(nComparisons,nConditions);
        NumPosSMALL       = nan(nComparisons,nConditions);
        OutOfNumSMALL     = nan(nComparisons,nConditions);
        reactionTimeSMALL = nan(nComparisons,nConditions);
        NumPosLARGE       = nan(nComparisons,nConditions);
        OutOfNumLARGE     = nan(nComparisons,nConditions);
        reactionTimeLARGE = nan(nComparisons,nConditions);
        for jj = 1:nConditions
            conditionName  = sprintf('%s%d','condition',jj);
            NumPos      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPos;
            OutOfNum    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNum;
            reactionTime(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTime;
            NumPosSMALL      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosSMALL;
            OutOfNumSMALL    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL;
            reactionTimeSMALL(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeSMALL;
            NumPosLARGE      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosLARGE;
            OutOfNumLARGE    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE;
            reactionTimeLARGE(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeLARGE;
        end
        
        % For this session, calculate across the conditions of this noise level.
        NumPosAll       = sum(NumPos,  2);
        OutOfNumAll     = sum(OutOfNum,2);
        reactionTimeAll = nanmean(reactionTime,2);
        NumPosAllSMALL       = sum(NumPosSMALL,  2);
        OutOfNumAllSMALL     = sum(OutOfNumSMALL,2);
        reactionTimeAllSMALL = nanmean(reactionTimeSMALL,2);
        NumPosAllLARGE       = sum(NumPosLARGE,  2);
        OutOfNumAllLARGE     = sum(OutOfNumLARGE,2);
        reactionTimeAllLARGE = nanmean(reactionTimeLARGE,2);
        
        % Store data.
        dataExperiment.performance.(noiseLevelName).NumPos      (:,ii) = NumPosAll;
        dataExperiment.performance.(noiseLevelName).OutOfNum    (:,ii) = OutOfNumAll;
        dataExperiment.performance.(noiseLevelName).reactionTime(:,ii) = reactionTimeAll;
        dataExperiment.performance.(noiseLevelName).NumPosSMALL      (:,ii) = NumPosAllSMALL;
        dataExperiment.performance.(noiseLevelName).OutOfNumSMALL    (:,ii) = OutOfNumAllSMALL;
        dataExperiment.performance.(noiseLevelName).reactionTimeSMALL(:,ii) = reactionTimeAllSMALL;
        dataExperiment.performance.(noiseLevelName).NumPosLARGE      (:,ii) = NumPosAllLARGE;
        dataExperiment.performance.(noiseLevelName).OutOfNumLARGE    (:,ii) = OutOfNumAllLARGE;
        dataExperiment.performance.(noiseLevelName).reactionTimeLARGE(:,ii) = reactionTimeAllLARGE;
    end
end

%% Convert comparison amounts from mm to degrees of visual angle for plotting
monitorDistance   = 75; %cm
monitorDistancemm = monitorDistance*10;
comparisonsDeg    = atand(comparisons/monitorDistancemm);

%% Plot performance for all sessions combined, for each noise level
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
threshold = nan(nNoiseLevels,1);
if plotFigures
    figure; hold on;
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Average performance across all sessions.
        performanceAllNumPos   = dataExperiment.performance.(noiseLevelName).NumPos;
        performanceAllOutOfNum = dataExperiment.performance.(noiseLevelName).OutOfNum;
        NumPos   = sum(performanceAllNumPos,  2);
        OutOfNum = sum(performanceAllOutOfNum,2); 
        performanceAll = NumPos./OutOfNum;
        
        % Plot data.
        plot(comparisonsDeg,performanceAll,'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        
        % Plot psychometric function fit.
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPos,OutOfNum);
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(nn) = thresholdthis;
    end
    % Plot parameters.
    if nNoiseLevels==2
        title({sprintf('%s%s%s%0.2f%s%0.2f',experimentName,subjectName, ...
            ': threshold0=',threshold(1),' threshold1=',threshold(2)),''});
        legend('Noise0 data','Noise0 fit','Noise1 data','Noise1 fit','Location','northwest')
    elseif nNoiseLevels==3
        title({sprintf('%s%s%s%0.2f%s%0.2f%s%0.2f',experimentName,subjectName, ...
            ': threshold0=',threshold(1),' threshold1=',threshold(2),' threshold2=',threshold(3)),''});
        legend('Noise0 data','Noise0 fit','Noise1 data','Noise1 fit','Noise2 data','Noise2 fit','Location','northwest')
    end
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Plot performance for each noise level as above, but separated by noise amount combination level (smaller, larger)
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
threshold = nan(nNoiseLevels+nNoiseLevels-1,1);
row = 1;
if plotFigures
    figure; hold on;
    
    % Average performance across all sessions, for the first noise level (NoiseLevel0).
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(1));
    performanceAllNumPos   = dataExperiment.performance.(noiseLevelName).NumPos;
    performanceAllOutOfNum = dataExperiment.performance.(noiseLevelName).OutOfNum;
    NumPos   = sum(performanceAllNumPos,  2);
    OutOfNum = sum(performanceAllOutOfNum,2);
    performanceAll = NumPos./OutOfNum;
    
    % Plot data.
    plot(comparisonsDeg,performanceAll,'o','MarkerFace',colors{1},'MarkerEdge',colors{1});
    
    % Plot psychometric function fit.
    [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPos,OutOfNum);
    plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{1});
    threshold(row) = thresholdthis;
    row = row+1;
    
    % Plot average performance across all sessions for the other noise levels, per noise amount combination level.   
    for nn = 2:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Average performance across all sessions.
        performanceAllNumPosSMALL   = dataExperiment.performance.(noiseLevelName).NumPosSMALL;
        performanceAllOutOfNumSMALL = dataExperiment.performance.(noiseLevelName).OutOfNumSMALL;
        performanceAllNumPosLARGE   = dataExperiment.performance.(noiseLevelName).NumPosLARGE;
        performanceAllOutOfNumLARGE = dataExperiment.performance.(noiseLevelName).OutOfNumLARGE;
        NumPosSMALL   = sum(performanceAllNumPosSMALL,  2);
        OutOfNumSMALL = sum(performanceAllOutOfNumSMALL,2); 
        performanceAllSMALL = NumPosSMALL./OutOfNumSMALL;
        NumPosLARGE   = sum(performanceAllNumPosLARGE,  2);
        OutOfNumLARGE = sum(performanceAllOutOfNumLARGE,2); 
        performanceAllLARGE = NumPosLARGE./OutOfNumLARGE;
        
        % Plot data: SMALL noise amount combination level.
        plot(comparisonsDeg,performanceAllSMALL,'o','MarkerFace','w','MarkerEdge',colors{nn});
        % Plot psychometric function fit: SMALL noise amount combination level.
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPosSMALL,OutOfNumSMALL);
        plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdthis;
        row = row+1;
        
        % Plot data: LARGE noise amount combination level.
        plot(comparisonsDeg,performanceAllLARGE,'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        % Plot psychometric function fit: LARGE noise amount combination level.
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPosLARGE,OutOfNumLARGE);
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdthis;
        row = row+1;  
    end
    % Plot parameters.
    if nNoiseLevels==2
        title({sprintf('%s%s%s%0.2f%s%0.2f%s%0.2f%s',experimentName,subjectName, ...
            ': threshold0=',threshold(1),' threshold1=(',threshold(2),',',threshold(3),')'),''});
        legend('Noise0 data','Noise0 fit','Noise1 SMALLER data','Noise1 SMALLER fit','Noise1 LARGER data','Noise1 LARGER fit','Location','northwest')
    elseif nNoiseLevels==3
        title({sprintf('%s%s%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName,subjectName, ...
            ': threshold0=',threshold(1),' threshold1=(',threshold(2),',',threshold(3),')', ...
            ' threshold2=(',threshold(4),',',threshold(5),')'),''});
        legend('Noise0 data','Noise0 fit','Noise1 SMALLER data','Noise1 SMALLER fit','Noise1 LARGER data','Noise1 LARGER fit', ...
            'Noise2 SMALLER data','Noise2 SMALLER fit','Noise2 LARGER data','Noise2 LARGER fit','Location','northwest')
    end
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Plot reaction times for all sessions combined, for each noise level
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
meanRT = nan(nNoiseLevels,1);
if plotFigures
    figure; hold on;
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Average reaction time across all sessions.
        reactionTime    = dataExperiment.performance.(noiseLevelName).reactionTime;
        reactionTimeAll = nanmean(reactionTime,2);
        
        % Plot data.
        plot(comparisonsDeg,reactionTimeAll,'Color',colors{nn});
        meanRT(nn) = round(nanmean(reactionTimeAll));
    end
    % Plot parameters.
    if nNoiseLevels==2
        title({sprintf('%s%s%s%d%s%d',experimentName,subjectName, ...
            ': mean0=',meanRT(1),' mean1=',meanRT(2)),''});
        legend('Noise0 data','Noise1 data','Location','northwest')
    elseif nNoiseLevels==3
        title({sprintf('%s%s%s%d%s%d%s%d',experimentName,subjectName, ...
            ': mean0=',meanRT(1),' mean1=',meanRT(2),' mean2=',meanRT(3)),''});
        legend('Noise0 data','Noise1 data','Noise2 data','Location','northwest')
    end
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Reaction time (ms)');
    axis([-Inf Inf -Inf Inf]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Save data analysis results

if saveData 
    % Save data struct.
    save(pathToOutputFile,'dataExperiment');
    fprintf('\nData was saved in:\n%s\n', pathToOutputFile);
end
end
%% End