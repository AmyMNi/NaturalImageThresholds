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
%   'experimentName' : (char)    Name of experiment folder (default: 'Experiment100')
%   'subjectName'    : (char)    Name of subject (default: 'test')
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
pathToFolder  = fullfile(getpref(projectName,'BaseDirAnalysis'),experimentName, ...
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
fn = fn([1:17 36:41]);
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
row = 1;
for ii = 1:length(fileInfo)
    
    % Exclude session for specified participant (see CM session notes for
    % details).
    if strcmp(subjectName,'CNSN_0009')
        if ii==3 || ii==5
            continue
        end
    end
    
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
        NumPosZERO1        = nan(nComparisons,nConditions);
        OutOfNumZERO1      = nan(nComparisons,nConditions);
        reactionTimeZERO1  = nan(nComparisons,nConditions);
        NumPosSMALL1       = nan(nComparisons,nConditions);
        OutOfNumSMALL1     = nan(nComparisons,nConditions);
        reactionTimeSMALL1 = nan(nComparisons,nConditions);
        NumPosLARGE1       = nan(nComparisons,nConditions);
        OutOfNumLARGE1     = nan(nComparisons,nConditions);
        reactionTimeLARGE1 = nan(nComparisons,nConditions);
        NumPosZERO2        = nan(nComparisons,nConditions);
        OutOfNumZERO2      = nan(nComparisons,nConditions);
        reactionTimeZERO2  = nan(nComparisons,nConditions);
        NumPosSMALL2       = nan(nComparisons,nConditions);
        OutOfNumSMALL2     = nan(nComparisons,nConditions);
        reactionTimeSMALL2 = nan(nComparisons,nConditions);
        NumPosLARGE2       = nan(nComparisons,nConditions);
        OutOfNumLARGE2     = nan(nComparisons,nConditions);
        reactionTimeLARGE2 = nan(nComparisons,nConditions);
        for jj = 1:nConditions
            conditionName  = sprintf('%s%d','condition',jj);
            NumPos      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPos;
            OutOfNum    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNum;
            reactionTime(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTime;
            NumPosZERO1       (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosZERO1;
            OutOfNumZERO1     (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumZERO1;
            reactionTimeZERO1 (:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeZERO1;
            NumPosSMALL1      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosSMALL1;
            OutOfNumSMALL1    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL1;
            reactionTimeSMALL1(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeSMALL1;
            NumPosLARGE1      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosLARGE1;
            OutOfNumLARGE1    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE1;
            reactionTimeLARGE1(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeLARGE1;
            NumPosZERO2       (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosZERO2;
            OutOfNumZERO2     (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumZERO2;
            reactionTimeZERO2 (:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeZERO2;
            NumPosSMALL2      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosSMALL2;
            OutOfNumSMALL2    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL2;
            reactionTimeSMALL2(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeSMALL2;
            NumPosLARGE2      (:,jj) = data.performance.(noiseLevelName).(conditionName).NumPosLARGE2;
            OutOfNumLARGE2    (:,jj) = data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE2;
            reactionTimeLARGE2(:,jj) = data.performance.(noiseLevelName).(conditionName).reactionTimeLARGE2;
        end
        
        % For this session, calculate across the conditions of this noise level.
        NumPosAll       = sum(NumPos,  2);
        OutOfNumAll     = sum(OutOfNum,2);
        reactionTimeAll = nanmean(reactionTime,2);
        NumPosAllZERO1        = sum(NumPosZERO1,  2);
        OutOfNumAllZERO1      = sum(OutOfNumZERO1,2);
        reactionTimeAllZERO1  = nanmean(reactionTimeZERO1,2);
        NumPosAllSMALL1       = sum(NumPosSMALL1,  2);
        OutOfNumAllSMALL1     = sum(OutOfNumSMALL1,2);
        reactionTimeAllSMALL1 = nanmean(reactionTimeSMALL1,2);
        NumPosAllLARGE1       = sum(NumPosLARGE1,  2);
        OutOfNumAllLARGE1     = sum(OutOfNumLARGE1,2);
        reactionTimeAllLARGE1 = nanmean(reactionTimeLARGE1,2);
        NumPosAllZERO2        = sum(NumPosZERO2,  2);
        OutOfNumAllZERO2      = sum(OutOfNumZERO2,2);
        reactionTimeAllZERO2  = nanmean(reactionTimeZERO2,2);
        NumPosAllSMALL2       = sum(NumPosSMALL2,  2);
        OutOfNumAllSMALL2     = sum(OutOfNumSMALL2,2);
        reactionTimeAllSMALL2 = nanmean(reactionTimeSMALL2,2);
        NumPosAllLARGE2       = sum(NumPosLARGE2,  2);
        OutOfNumAllLARGE2     = sum(OutOfNumLARGE2,2);
        reactionTimeAllLARGE2 = nanmean(reactionTimeLARGE2,2);
        
        % Store data.
        dataExperiment.performance.(noiseLevelName).NumPos            (:,row) = NumPosAll;
        dataExperiment.performance.(noiseLevelName).OutOfNum          (:,row) = OutOfNumAll;
        dataExperiment.performance.(noiseLevelName).reactionTime      (:,row) = reactionTimeAll;
        dataExperiment.performance.(noiseLevelName).NumPosZERO1       (:,row) = NumPosAllZERO1;
        dataExperiment.performance.(noiseLevelName).OutOfNumZERO1     (:,row) = OutOfNumAllZERO1;
        dataExperiment.performance.(noiseLevelName).reactionTimeZERO1 (:,row) = reactionTimeAllZERO1;
        dataExperiment.performance.(noiseLevelName).NumPosSMALL1      (:,row) = NumPosAllSMALL1;
        dataExperiment.performance.(noiseLevelName).OutOfNumSMALL1    (:,row) = OutOfNumAllSMALL1;
        dataExperiment.performance.(noiseLevelName).reactionTimeSMALL1(:,row) = reactionTimeAllSMALL1;
        dataExperiment.performance.(noiseLevelName).NumPosLARGE1      (:,row) = NumPosAllLARGE1;
        dataExperiment.performance.(noiseLevelName).OutOfNumLARGE1    (:,row) = OutOfNumAllLARGE1;
        dataExperiment.performance.(noiseLevelName).reactionTimeLARGE1(:,row) = reactionTimeAllLARGE1;
        dataExperiment.performance.(noiseLevelName).NumPosZERO2       (:,row) = NumPosAllZERO2;
        dataExperiment.performance.(noiseLevelName).OutOfNumZERO2     (:,row) = OutOfNumAllZERO2;
        dataExperiment.performance.(noiseLevelName).reactionTimeZERO2 (:,row) = reactionTimeAllZERO2;
        dataExperiment.performance.(noiseLevelName).NumPosSMALL2      (:,row) = NumPosAllSMALL2;
        dataExperiment.performance.(noiseLevelName).OutOfNumSMALL2    (:,row) = OutOfNumAllSMALL2;
        dataExperiment.performance.(noiseLevelName).reactionTimeSMALL2(:,row) = reactionTimeAllSMALL2;
        dataExperiment.performance.(noiseLevelName).NumPosLARGE2      (:,row) = NumPosAllLARGE2;
        dataExperiment.performance.(noiseLevelName).OutOfNumLARGE2    (:,row) = OutOfNumAllLARGE2;
        dataExperiment.performance.(noiseLevelName).reactionTimeLARGE2(:,row) = reactionTimeAllLARGE2;
    end
    row = row+1;
end

%% Convert comparison amounts from mm to degrees of visual angle for plotting
monitorDistance   = 1.2; % meters (distance of the scene in iset3d)
monitorDistancemm = monitorDistance*1000; % convert to mm
comparisonsDeg    = atand(comparisons/monitorDistancemm);

%% Plot performance for all sessions combined, for each noise level
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
threshold = nan(nNoiseLevels,1);
pse       = nan(nNoiseLevels,1);
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
        [xx,FittedCurve,thresholdthis,psethis] = fitPsychometric(comparisonsDeg,NumPos,OutOfNum);
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(nn) = thresholdthis;
        pse(nn) = psethis;
    end
    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f',experimentName,subjectName, ...
        ': threshold0=',threshold(1),' threshold1=',threshold(2),' threshold2=',threshold(3)), ...
        sprintf('%s%0.2f%s%0.2f%s%0.2f','pse0=',pse(1),' pse1=',pse(2),' pse2=',pse(3))},'Interpreter','none');
    legend('Noise0 data','Noise0 fit','Noise1 data','Noise1 fit','Noise2 data','Noise2 fit','Location','northwest')
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Plot performance for each noise level as above, but separated by NoiseAmounts1 combination level (smaller, larger)
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
        performanceAllNumPosSMALL1   = dataExperiment.performance.(noiseLevelName).NumPosSMALL1;
        performanceAllOutOfNumSMALL1 = dataExperiment.performance.(noiseLevelName).OutOfNumSMALL1;
        performanceAllNumPosLARGE1   = dataExperiment.performance.(noiseLevelName).NumPosLARGE1;
        performanceAllOutOfNumLARGE1 = dataExperiment.performance.(noiseLevelName).OutOfNumLARGE1;
        NumPosSMALL   = sum(performanceAllNumPosSMALL1,  2);
        OutOfNumSMALL = sum(performanceAllOutOfNumSMALL1,2); 
        performanceAllSMALL = NumPosSMALL./OutOfNumSMALL;
        NumPosLARGE   = sum(performanceAllNumPosLARGE1,  2);
        OutOfNumLARGE = sum(performanceAllOutOfNumLARGE1,2); 
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
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName,subjectName, ...
        ': threshold0=',threshold(1),' threshold1=(',threshold(2),',',threshold(3),')', ...
        ' threshold2=(',threshold(4),',',threshold(5),')'),''},'Interpreter','none');
    legend('Noise0 data','Noise0 fit','Noise1 SMALLER data','Noise1 SMALLER fit','Noise1 LARGER data','Noise1 LARGER fit', ...
        'Noise2 SMALLER data','Noise2 SMALLER fit','Noise2 LARGER data','Noise2 LARGER fit','Location','northwest')
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Plot performance for each noise level as above, but separated by NoiseAmounts2 combination level (smaller, larger)
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
threshold = nan(nNoiseLevels+nNoiseLevels-2,1);
row = 1;
if plotFigures
    figure; hold on;
    
    % Average performance across all sessions, for NoiseLevel0 and NoiseLevel1.
    for nn = 1:2
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
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
        threshold(row) = thresholdthis;
        row = row+1;
    end
    
    % Plot average performance across all sessions for NoiseLevel2, per noise amount combination level.   
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(3));
    
    % Average performance across all sessions.
    performanceAllNumPosSMALL2   = dataExperiment.performance.(noiseLevelName).NumPosSMALL2;
    performanceAllOutOfNumSMALL2 = dataExperiment.performance.(noiseLevelName).OutOfNumSMALL2;
    performanceAllNumPosLARGE2   = dataExperiment.performance.(noiseLevelName).NumPosLARGE2;
    performanceAllOutOfNumLARGE2 = dataExperiment.performance.(noiseLevelName).OutOfNumLARGE2;
    NumPosSMALL   = sum(performanceAllNumPosSMALL2,  2);
    OutOfNumSMALL = sum(performanceAllOutOfNumSMALL2,2);
    performanceAllSMALL = NumPosSMALL./OutOfNumSMALL;
    NumPosLARGE   = sum(performanceAllNumPosLARGE2,  2);
    OutOfNumLARGE = sum(performanceAllOutOfNumLARGE2,2);
    performanceAllLARGE = NumPosLARGE./OutOfNumLARGE;
    
    % Plot data: SMALL noise amount combination level.
    plot(comparisonsDeg,performanceAllSMALL,'o','MarkerFace','w','MarkerEdge',colors{3});
    % Plot psychometric function fit: SMALL noise amount combination level.
    [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPosSMALL,OutOfNumSMALL);
    plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{3});
    threshold(row) = thresholdthis;
    row = row+1;
    
    % Plot data: LARGE noise amount combination level.
    plot(comparisonsDeg,performanceAllLARGE,'o','MarkerFace',colors{3},'MarkerEdge',colors{3});
    % Plot psychometric function fit: LARGE noise amount combination level.
    [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPosLARGE,OutOfNumLARGE);
    plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{3});
    threshold(row) = thresholdthis;

    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s',experimentName,subjectName, ...
        ': threshold0=',threshold(1),' threshold1=',threshold(2), ...
        ' threshold2=(',threshold(3),',',threshold(4),')'),''},'Interpreter','none');
    legend('Noise0 data','Noise0 fit','Noise1 data','Noise1 fit', ...
        'Noise2 SMALLER data','Noise2 SMALLER fit','Noise2 LARGER data','Noise2 LARGER fit','Location','northwest')
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Plot performance for each noise level as above, but separated by NoiseAmounts1 combination level (0, all)
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
        performanceAllNumPosZERO1   = dataExperiment.performance.(noiseLevelName).NumPosZERO1;
        performanceAllOutOfNumZERO1 = dataExperiment.performance.(noiseLevelName).OutOfNumZERO1;
        performanceAllNumPos        = dataExperiment.performance.(noiseLevelName).NumPos;
        performanceAllOutOfNum      = dataExperiment.performance.(noiseLevelName).OutOfNum;
        NumPosZERO   = sum(performanceAllNumPosZERO1,  2);
        OutOfNumZERO = sum(performanceAllOutOfNumZERO1,2); 
        performanceAllZERO = NumPosZERO./OutOfNumZERO;
        NumPos   = sum(performanceAllNumPos,  2);
        OutOfNum = sum(performanceAllOutOfNum,2); 
        performanceAll = NumPos./OutOfNum;
        
        % Plot data: ZERO noise amount combination level.
        plot(comparisonsDeg,performanceAllZERO,'o','MarkerFace','w','MarkerEdge',colors{nn});
        % Plot psychometric function fit: ZERO noise amount combination level.
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPosZERO,OutOfNumZERO);
        plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdthis;
        row = row+1;
        
        % Plot data: ALL noise amount combination level.
        plot(comparisonsDeg,performanceAll,'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        % Plot psychometric function fit: ALL noise amount combination level.
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPos,OutOfNum);
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdthis;
        row = row+1;  
    end
    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName,subjectName, ...
        ': threshold0=',threshold(1),' threshold1=(',threshold(2),',',threshold(3),')', ...
        ' threshold2=(',threshold(4),',',threshold(5),')'),''},'Interpreter','none');
    legend('Noise0 data','Noise0 fit','Noise1 ZERO NOISE data','Noise1 ZERO NOISE fit','Noise1 ALL data','Noise1 ALL fit', ...
        'Noise2 ZERO NOISE data','Noise2 ZERO NOISE fit','Noise2 ALL data','Noise2 ALL fit','Location','northwest')
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Proportion chose comparison as rightward');
    axis([-Inf Inf 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',comparisonsDeg);
    set(gca,'XTickLabel',comparisonsDeg);
    xtickformat('%.1f');
    box off; hold off;
end

%% Plot performance for each noise level as above, but separated by NoiseAmounts2 combination level (0, all)
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
threshold = nan(nNoiseLevels+nNoiseLevels-2,1);
row = 1;
if plotFigures
    figure; hold on;
    
    % Average performance across all sessions, for NoiseLevel0 and NoiseLevel1.
    for nn = 1:2
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
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
        threshold(row) = thresholdthis;
        row = row+1;
    end
    
    % Plot average performance across all sessions for NoiseLevel2, per noise amount combination level.   
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(3));
    
    % Average performance across all sessions.
    performanceAllNumPosZERO2   = dataExperiment.performance.(noiseLevelName).NumPosZERO2;
    performanceAllOutOfNumZERO2 = dataExperiment.performance.(noiseLevelName).OutOfNumZERO2;
    performanceAllNumPos        = dataExperiment.performance.(noiseLevelName).NumPos;
    performanceAllOutOfNum      = dataExperiment.performance.(noiseLevelName).OutOfNum;
    NumPosZERO   = sum(performanceAllNumPosZERO2,  2);
    OutOfNumZERO = sum(performanceAllOutOfNumZERO2,2);
    performanceAllZERO = NumPosZERO./OutOfNumZERO;
    NumPos   = sum(performanceAllNumPos,  2);
    OutOfNum = sum(performanceAllOutOfNum,2);
    performanceAll = NumPos./OutOfNum;
    
    % Plot data: ZERO noise amount combination level.
    plot(comparisonsDeg,performanceAllZERO,'o','MarkerFace','w','MarkerEdge',colors{3});
    % Plot psychometric function fit: ZERO noise amount combination level.
    [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPosZERO,OutOfNumZERO);
    plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{3});
    threshold(row) = thresholdthis;
    row = row+1;
    
    % Plot data: ALL noise amount combination level.
    plot(comparisonsDeg,performanceAll,'o','MarkerFace',colors{3},'MarkerEdge',colors{3});
    % Plot psychometric function fit: ALL noise amount combination level.
    [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,NumPos,OutOfNum);
    plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{3});
    threshold(row) = thresholdthis;

    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s',experimentName,subjectName, ...
        ': threshold0=',threshold(1),' threshold1=',threshold(2), ...
        ' threshold2=(',threshold(3),',',threshold(4),')'),''},'Interpreter','none');
    legend('Noise0 data','Noise0 fit','Noise1 data','Noise1 fit', ...
        'Noise2 ZERO NOISE data','Noise2 ZERO NOISE fit','Noise2 ALL data','Noise2 ALL fit','Location','northwest')
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
    title({sprintf('%s_%s%s%d%s%d%s%d',experimentName,subjectName, ...
        ': mean0=',meanRT(1),' mean1=',meanRT(2),' mean2=',meanRT(3)),''},'Interpreter','none');
    legend('Noise0 data','Noise1 data','Noise2 data','Location','southwest')
    xlabel(sprintf('Comparison offset rightward (deg)'));
    ylabel('Reaction time (ms)');
    axis([-Inf Inf 0 Inf]);
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