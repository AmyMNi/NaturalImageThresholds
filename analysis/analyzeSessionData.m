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
% Get the identifies of noise levels.
noiseLevels      = unique(data.imageNoiseLevel);
nNoiseLevels     = numel(noiseLevels);
data.noiseLevels = noiseLevels;

% Per noise level, get the identities of noise amounts and median noise amount combination.
noiseAmountsDiffMed = zeros(nNoiseLevels,1);
noiseAmounts = data.imageNoiseAmount;
% Skip first NoiseLevel (NoiseLevel0) because the median noise amount combination is zero.
for nn = 2:nNoiseLevels
    % Get image indices for this noise level.
    noiseLevelthis = noiseLevels(nn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: NoiseLevel1 and NoiseLevel2 are separate (though both include noise 0)
    
    % Always include noise level 0 in each noise level.
    imagesNoise = data.imageNoiseLevel==noiseLevelthis | data.imageNoiseLevel==0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the noise amounts for this noise level.
    noiseAmountsthis = unique(noiseAmounts(imagesNoise));
    data.noiseAmounts{nn,1} = noiseAmountsthis;
    noiseAmountsCombo = nchoosek(noiseAmountsthis,2);
    noiseAmountsAll = [noiseAmountsCombo; [noiseAmountsthis noiseAmountsthis]];
    noiseAmountsDiff  = abs(noiseAmountsAll(:,1)-noiseAmountsAll(:,2));
    noiseAmountsDiffMed(nn) = median(noiseAmountsDiff);
end

% Get the identifies of conditions (same per noise level).
conditions      = unique(data.imageCondition);
nConditions     = numel(conditions);
data.conditions = conditions;

% Get the identities of comparisons (same per condition).
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

    % Get the image indices, target offset amount, noise amount, observer
    % response, and reaction time per trial.
    imagesN       = data.trialOrder(trialsNoise,:);
    offsetsN      = data.trialOrderComparison(trialsNoise,:);
    noiseAmountsN = data.trialNoiseAmount(trialsNoise,:);
    responsesN    = data.selectedResponse(trialsNoise);
    rtbeg = datetime(data.reactionTimeStart(trialsNoise),'InputFormat','MM/dd/yyyy HH:mm:ss.SSS','Format','HH:mm:ss.SSS');
    rtend = datetime(data.reactionTimeEnd(trialsNoise),  'InputFormat','MM/dd/yyyy HH:mm:ss.SSS','Format','HH:mm:ss.SSS');
    reactionTimeN = milliseconds(diff([rtbeg rtend],1,2));
    
    % Analyze each condition separately.
    for ii = 1:nConditions
        
        % Get the center position for this condition.
        centerpos = conditions(ii);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: NoiseLevel1 and NoiseLevel2 are separate (though both include noise 0)
        
        % For the center position of this condition, get the pool of
        % images (of the various noise amounts for this noise level,
        % including the noise amount of 0 from Noise Level 0).
        centerPool = find(data.imageComparison==0 & ...
            data.imageCondition==centerpos & ...
            (data.imageNoiseLevel==noiseLevelthis | data.imageNoiseLevel==0));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get parameters per trial.
        trialsC       = any(ismember(imagesN,centerPool),2);
        offsetsC      = offsetsN(trialsC,:);
        comparisonsC  = sum(offsetsC,2);
        noiseAmountsC = noiseAmountsN(trialsC,:); 
        responsesC    = responsesN(trialsC);
        reactionTimeC = reactionTimeN(trialsC);
        
        % Calculate per comparison amount: 
        %   The number of trials a positive response was given
        % 	The total number of trials
        %   The mean reaction time (outliers excluded with Tukey method)
        performanceNumPos   = nan(nComparisons,1);
        performanceOutOfNum = nan(nComparisons,1);
        reactionTime        = nan(nComparisons,1);
        % Calculate same as above, but for each noise amount combination level (smaller, larger).
        performanceNumPosSMALL   = nan(nComparisons,1);
        performanceNumPosLARGE   = nan(nComparisons,1);
        performanceOutOfNumSMALL = nan(nComparisons,1);
        performanceOutOfNumLARGE = nan(nComparisons,1);
        reactionTimeSMALL        = nan(nComparisons,1);
        reactionTimeLARGE        = nan(nComparisons,1);
        for jj = 1:nComparisons
            comparisonthis   = comparisons(jj);
            offsetsthis      = offsetsC     (comparisonsC==comparisonthis,:);
            noiseAmountsthis = noiseAmountsC(comparisonsC==comparisonthis,:);
            responsesthis    = responsesC   (comparisonsC==comparisonthis);
            reactionTimethis = reactionTimeC(comparisonsC==comparisonthis);
            
            % Calculate proportion observer chose the comparison as rightward.
            if comparisonthis==0
                choseRight = find(responsesthis==2);
            else
                choseRight = [find(offsetsthis(:,1)==0 & responsesthis==2); ...
                              find(offsetsthis(:,2)==0 & responsesthis==1)];
            end
            performanceNumPos(jj)   = numel(choseRight);
            performanceOutOfNum(jj) = numel(responsesthis);
            reactionTimeTukey       = calcTukey(reactionTimethis);
            reactionTime(jj)        = nanmean(reactionTimeTukey);   
            
            % Calculate for each noise amount combination level (smaller, larger).
            if nn>1
                noiseAmountsDiffthis = abs(noiseAmountsthis(:,1)-noiseAmountsthis(:,2));
                idxSMALL = find(noiseAmountsDiffthis<=noiseAmountsDiffMed(nn));
                idxLARGE = find(noiseAmountsDiffthis> noiseAmountsDiffMed(nn));
                performanceNumPosSMALL(jj)   = numel(intersect(choseRight,idxSMALL));
                performanceNumPosLARGE(jj)   = numel(intersect(choseRight,idxLARGE));
                performanceOutOfNumSMALL(jj) = numel(idxSMALL);
                performanceOutOfNumLARGE(jj) = numel(idxLARGE);
                reactionTimeSMALL(jj)        = nanmean(reactionTimeTukey(idxSMALL));
                reactionTimeLARGE(jj)        = nanmean(reactionTimeTukey(idxLARGE));
            end
        end
        
        % Save noise amounts, performance and reaction times for this condition (in this noise level).
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevelthis);
        conditionName  = sprintf('%s%d','condition',ii);
        data.performance.(noiseLevelName).(conditionName).NumPos       = performanceNumPos;
        data.performance.(noiseLevelName).(conditionName).OutOfNum     = performanceOutOfNum;
        data.performance.(noiseLevelName).(conditionName).reactionTime = reactionTime;
        data.performance.(noiseLevelName).(conditionName).NumPosSMALL       = performanceNumPosSMALL;
        data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL     = performanceOutOfNumSMALL;
        data.performance.(noiseLevelName).(conditionName).reactionTimeSMALL = reactionTimeSMALL;
        data.performance.(noiseLevelName).(conditionName).NumPosLARGE       = performanceNumPosLARGE;
        data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE     = performanceOutOfNumLARGE;
        data.performance.(noiseLevelName).(conditionName).reactionTimeLARGE = reactionTimeLARGE;
    end
end

%% Convert comparison amounts from mm to degrees of visual angle for plotting

monitorDistance   = 75; %cm
monitorDistancemm = monitorDistance*10;
comparisonsDeg    = atand(comparisons/monitorDistancemm);

%% Plot performance for all conditions combined, for each noise level
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
threshold = nan(nNoiseLevels,1);
if plotFigures
    figure; hold on;
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Calculate performance across all conditions.
        performanceAllNumPos   = nan(nComparisons,nConditions);
        performanceAllOutOfNum = nan(nComparisons,nConditions);
        for ii = 1:nConditions
            conditionName  = sprintf('%s%d','condition',ii);
            performanceAllNumPos  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPos;
            performanceAllOutOfNum(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNum;
        end
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
        title({sprintf('%s%s%s%d%s%0.2f%s%0.2f',experimentName,subjectName,'\_', sessionNumber, ...
            ': threshold0=',threshold(1),' threshold1=',threshold(2)),''});
        legend('Noise0 data','Noise0 fit','Noise1 data','Noise1 fit','Location','northwest')
    elseif nNoiseLevels==3
        title({sprintf('%s%s%s%d%s%0.2f%s%0.2f%s%0.2f',experimentName,subjectName,'\_', sessionNumber, ...
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
    
    % Plot performance across all conditions as above, for the first noise level (NoiseLevel0).
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(1));
    performanceAllNumPos   = nan(nComparisons,nConditions);
    performanceAllOutOfNum = nan(nComparisons,nConditions);
    for ii = 1:nConditions
        conditionName  = sprintf('%s%d','condition',ii);
        performanceAllNumPos  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPos;
        performanceAllOutOfNum(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNum;
    end
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
      
    % Plot performance for the other noise levels, per noise amount combination level.        
    for nn = 2:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));

        % Calculate performance across all conditions, per noise amount combination level.
        performanceAllNumPosSMALL   = nan(nComparisons,nConditions);
        performanceAllOutOfNumSMALL = nan(nComparisons,nConditions);
        performanceAllNumPosLARGE   = nan(nComparisons,nConditions);
        performanceAllOutOfNumLARGE = nan(nComparisons,nConditions);
        for ii = 1:nConditions
            conditionName  = sprintf('%s%d','condition',ii);
            performanceAllNumPosSMALL  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPosSMALL;
            performanceAllOutOfNumSMALL(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL;
            performanceAllNumPosLARGE  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPosLARGE;
            performanceAllOutOfNumLARGE(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE;
        end
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
        title({sprintf('%s%s%s%d%s%0.2f%s%0.2f%s%0.2f%s',experimentName,subjectName,'\_', sessionNumber, ...
            ': threshold0=',threshold(1),' threshold1=(',threshold(2),',',threshold(3),')'),''});
        legend('Noise0 data','Noise0 fit','Noise1 SMALLER data','Noise1 SMALLER fit','Noise1 LARGER data','Noise1 LARGER fit','Location','northwest')
    elseif nNoiseLevels==3
        title({sprintf('%s%s%s%d%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName,subjectName,'\_', sessionNumber, ...
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

%% Plot reaction times for all condition combined, for each noise level
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% Plot all noise levels.
meanRT = nan(nNoiseLevels,1);
if plotFigures
    figure; hold on;
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Calculate mean reaction time across all conditions.
        reactionTimeAll = nan(nComparisons,nConditions);
        for ii = 1:nConditions
            conditionName = sprintf('%s%d','condition',ii);
            reactionTimeAll(:,ii) = data.performance.(noiseLevelName).(conditionName).reactionTime;
        end
        reactionTime = nanmean(reactionTimeAll,2);
        
        % Plot data.
        plot(comparisonsDeg,reactionTime,'Color',colors{nn});
        meanRT(nn) = round(nanmean(reactionTime));
    end
    % Plot parameters.
    if nNoiseLevels==2
        title({sprintf('%s%s%s%d%s%d%s%d',experimentName,subjectName,'\_', sessionNumber, ...
            ': mean0=',meanRT(1),' mean1=',meanRT(2)),''});
        legend('Noise0 data','Noise1 data','Location','northwest')
    elseif nNoiseLevels==3
        title({sprintf('%s%s%s%d%s%d%s%d%s%d',experimentName,subjectName,'\_', sessionNumber, ...
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
    % Save data struct (with analysis result additions).
    save(pathToOutputFile,'data');
    fprintf('\nData was saved in:\n%s\n', pathToOutputFile);
end
end
%% End