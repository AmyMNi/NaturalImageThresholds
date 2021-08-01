function data = analyzeSessionData(varargin)
%analyzeSessionData
%
% Usage:
%   data = analyzeSessionData('experimentName', 'Experiment100', ...
%                             'subjectName', 'test', ...
%                             'sessionNumber', 1);
%
% Description:
%   Analyze psychophysical data from a single data collection session. Save
%   the results (struct 'data') in the specified output folder.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment100')
%   'subjectName'    : (string)  Name of subject (default: 'test')
%   'sessionNumber'  : (scalar)  Number of session (default: 1)
%   'plotFigures'    : (logical) Plot figures if option is on (default: true)
%   'saveData'       : (logical) Save data if option is on (default: true)
%
% History:
%   06/10/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.addParameter('subjectName', 'test', @ischar);
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
noiseLevels      = unique(data.trialNoiseLevel);
nNoiseLevels     = numel(noiseLevels);
data.noiseLevels = noiseLevels;

% Per noise level, get the identities of noise amounts and median noise amount combination
% based on the noise amounts available in all of the images (this should be
% fine - the noise amounts tested are not user-controlled; also, this
% ensures that all possible noise amounts are included).
noiseAmounts1DiffMed = zeros(nNoiseLevels,1);
noiseAmounts2DiffMed = zeros(nNoiseLevels,1);
noiseAmounts1 = data.imageNoiseAmount1;
noiseAmounts2 = data.imageNoiseAmount2;
% Skip first NoiseLevel (NoiseLevel0) because the median noise amount combination is zero.
for nn = 2:nNoiseLevels
    % Get image indices for this noise level.
    noiseLevelthis = noiseLevels(nn);

    % Always include noise level 0 in each noise level.
    imagesNoise = data.imageNoiseLevel==noiseLevelthis | data.imageNoiseLevel==0;

    % Get noiseAmounts1 for this noise level.
    noiseAmounts1this = unique(noiseAmounts1(imagesNoise));
    data.noiseAmounts1{nn,1} = noiseAmounts1this;
    noiseAmounts1Combo = nchoosek(noiseAmounts1this,2);
    noiseAmounts1All = [noiseAmounts1Combo; [noiseAmounts1this noiseAmounts1this]];
    noiseAmounts1Diff = abs(noiseAmounts1All(:,1)-noiseAmounts1All(:,2));
    noiseAmounts1DiffMed(nn) = median(noiseAmounts1Diff);
    
    if noiseLevelthis == 2
        % Get noiseAmounts2 for this noise level.
        noiseAmounts2this = unique(noiseAmounts2(imagesNoise));
        data.noiseAmounts2{nn,1} = noiseAmounts2this;
        noiseAmounts2Combo = nchoosek(noiseAmounts2this,2);
        noiseAmounts2All = [noiseAmounts2Combo; [noiseAmounts2this noiseAmounts2this]];
        noiseAmounts2Diff = abs(noiseAmounts2All(:,1)-noiseAmounts2All(:,2));
        noiseAmounts2DiffMed(nn) = median(noiseAmounts2Diff);
    end
end

% Get the identifies of conditions (same per noise level) based on the
% conditions available in all of the images (this should be fine - the
% conditions tested are not user-controlled).
conditions      = unique(data.imageCondition);
nConditions     = numel(conditions);
data.conditions = conditions;

% Get the identities of comparisons (same per condition).
comparisons      = unique(data.trialOrderComparison);
nComparisons     = numel(comparisons);
data.comparisons = comparisons;

%% Analyze performance on each noise level and condition separately
%
% Analyze each noise level separately.
for nn = 1:nNoiseLevels
    
    % Get trial indices for this noise level.
    noiseLevelthis = noiseLevels(nn);
    trialsNoise    = data.trialNoiseLevel==noiseLevelthis;

    % Get the image indices, target offset amount, noise amounts, observer
    % response, and reaction time per trial.
    imagesN        = data.trialOrder(trialsNoise,:);
    offsetsN       = data.trialOrderComparison(trialsNoise,:);
    noiseAmounts1N = data.trialNoiseAmount1(trialsNoise,:);
    noiseAmounts2N = data.trialNoiseAmount2(trialsNoise,:);
    responsesN     = data.selectedResponse(trialsNoise);
    rtbeg = datetime(data.reactionTimeStart(trialsNoise),'InputFormat','MM/dd/yyyy HH:mm:ss.SSS','Format','HH:mm:ss.SSS');
    rtend = datetime(data.reactionTimeEnd(trialsNoise),  'InputFormat','MM/dd/yyyy HH:mm:ss.SSS','Format','HH:mm:ss.SSS');
    reactionTimeN = milliseconds(diff([rtbeg rtend],1,2));
    
    % Analyze each condition separately.
    for ii = 1:nConditions
        
        % Get the center position for this condition.
        centerpos = conditions(ii);

        % For the center position of this condition, get the pool of
        % images (of the various noise amounts for this noise level,
        % including the noise amount of 0 from Noise Level 0).
        centerPool = find(data.imageComparison==0 & ...
            data.imageCondition==centerpos & ...
            (data.imageNoiseLevel==noiseLevelthis | data.imageNoiseLevel==0));

        % Get parameters per trial.
        trialsC        = any(ismember(imagesN,centerPool),2);
        offsetsC       = offsetsN(trialsC,:);
        comparisonsC   = sum(offsetsC,2);
        noiseAmounts1C = noiseAmounts1N(trialsC,:); 
        noiseAmounts2C = noiseAmounts2N(trialsC,:);
        responsesC     = responsesN(trialsC);
        reactionTimeC  = reactionTimeN(trialsC);
        
        % Calculate per comparison amount: 
        %   The number of trials a positive response was given
        % 	The total number of trials
        %   The mean reaction time (outliers excluded with Tukey method)
        performanceNumPos   = nan(nComparisons,1);
        performanceOutOfNum = nan(nComparisons,1);
        reactionTime        = nan(nComparisons,1);
        % Calculate same as above, but for each noise amount combination level (smaller, larger) for NoiseAmounts1.
        performanceNumPosSMALL1   = nan(nComparisons,1);
        performanceNumPosLARGE1   = nan(nComparisons,1);
        performanceOutOfNumSMALL1 = nan(nComparisons,1);
        performanceOutOfNumLARGE1 = nan(nComparisons,1);
        reactionTimeSMALL1        = nan(nComparisons,1);
        reactionTimeLARGE1        = nan(nComparisons,1);
        % Calculate same as above, but for each noise amount combination level (smaller, larger) for NoiseAmounts2.
        performanceNumPosSMALL2   = nan(nComparisons,1);
        performanceNumPosLARGE2   = nan(nComparisons,1);
        performanceOutOfNumSMALL2 = nan(nComparisons,1);
        performanceOutOfNumLARGE2 = nan(nComparisons,1);
        reactionTimeSMALL2        = nan(nComparisons,1);
        reactionTimeLARGE2        = nan(nComparisons,1);
        for jj = 1:nComparisons
            comparisonthis    = comparisons(jj);
            offsetsthis       = offsetsC      (comparisonsC==comparisonthis,:);
            noiseAmounts1this = noiseAmounts1C(comparisonsC==comparisonthis,:);
            noiseAmounts2this = noiseAmounts2C(comparisonsC==comparisonthis,:);
            responsesthis     = responsesC    (comparisonsC==comparisonthis);
            reactionTimethis  = reactionTimeC (comparisonsC==comparisonthis);
            
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
                noiseAmounts1Diffthis = abs(noiseAmounts1this(:,1)-noiseAmounts1this(:,2));
                idxSMALL = find(noiseAmounts1Diffthis<=noiseAmounts1DiffMed(nn));
                idxLARGE = find(noiseAmounts1Diffthis> noiseAmounts1DiffMed(nn));
                performanceNumPosSMALL1(jj)   = numel(intersect(choseRight,idxSMALL));
                performanceNumPosLARGE1(jj)   = numel(intersect(choseRight,idxLARGE));
                performanceOutOfNumSMALL1(jj) = numel(idxSMALL);
                performanceOutOfNumLARGE1(jj) = numel(idxLARGE);
                reactionTimeSMALL1(jj)        = nanmean(reactionTimeTukey(idxSMALL));
                reactionTimeLARGE1(jj)        = nanmean(reactionTimeTukey(idxLARGE));
            end
            if nn==3
                noiseAmounts2Diffthis = abs(noiseAmounts2this(:,1)-noiseAmounts2this(:,2));
                idxSMALL = find(noiseAmounts2Diffthis<=noiseAmounts2DiffMed(nn));
                idxLARGE = find(noiseAmounts2Diffthis> noiseAmounts2DiffMed(nn));
                performanceNumPosSMALL2(jj)   = numel(intersect(choseRight,idxSMALL));
                performanceNumPosLARGE2(jj)   = numel(intersect(choseRight,idxLARGE));
                performanceOutOfNumSMALL2(jj) = numel(idxSMALL);
                performanceOutOfNumLARGE2(jj) = numel(idxLARGE);
                reactionTimeSMALL2(jj)        = nanmean(reactionTimeTukey(idxSMALL));
                reactionTimeLARGE2(jj)        = nanmean(reactionTimeTukey(idxLARGE));     
            end
        end
        
        % Save noise amounts, performance and reaction times for this condition (in this noise level).
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevelthis);
        conditionName  = sprintf('%s%d','condition',ii);
        data.performance.(noiseLevelName).(conditionName).NumPos       = performanceNumPos;
        data.performance.(noiseLevelName).(conditionName).OutOfNum     = performanceOutOfNum;
        data.performance.(noiseLevelName).(conditionName).reactionTime = reactionTime;
        data.performance.(noiseLevelName).(conditionName).NumPosSMALL1       = performanceNumPosSMALL1;
        data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL1     = performanceOutOfNumSMALL1;
        data.performance.(noiseLevelName).(conditionName).reactionTimeSMALL1 = reactionTimeSMALL1;
        data.performance.(noiseLevelName).(conditionName).NumPosLARGE1       = performanceNumPosLARGE1;
        data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE1     = performanceOutOfNumLARGE1;
        data.performance.(noiseLevelName).(conditionName).reactionTimeLARGE1 = reactionTimeLARGE1;
        
        data.performance.(noiseLevelName).(conditionName).NumPosSMALL2       = performanceNumPosSMALL2;
        data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL2     = performanceOutOfNumSMALL2;
        data.performance.(noiseLevelName).(conditionName).reactionTimeSMALL2 = reactionTimeSMALL2;
        data.performance.(noiseLevelName).(conditionName).NumPosLARGE2       = performanceNumPosLARGE2;
        data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE2     = performanceOutOfNumLARGE2;
        data.performance.(noiseLevelName).(conditionName).reactionTimeLARGE2 = reactionTimeLARGE2;
    end
end

%% Convert comparison amounts from mm to degrees of visual angle for plotting

monitorDistance   = 1.2; % meters (distance of the scene in iset3d)
monitorDistancemm = monitorDistance*1000; % convert to mm
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
    title({sprintf('%s%s%s%d%s%0.2f%s%0.2f%s%0.2f',experimentName,subjectName,'\_', sessionNumber, ...
        ': threshold0=',threshold(1),' threshold1=',threshold(2),' threshold2=',threshold(3)),''});
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
      
    % Plot performance across all conditions as above, for the other noise levels.        
    for nn = 2:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));

        % Calculate performance across all conditions, per noise amount combination level.
        performanceAllNumPosSMALL1   = nan(nComparisons,nConditions);
        performanceAllOutOfNumSMALL1 = nan(nComparisons,nConditions);
        performanceAllNumPosLARGE1   = nan(nComparisons,nConditions);
        performanceAllOutOfNumLARGE1 = nan(nComparisons,nConditions);
        for ii = 1:nConditions
            conditionName  = sprintf('%s%d','condition',ii);
            performanceAllNumPosSMALL1  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPosSMALL1;
            performanceAllOutOfNumSMALL1(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL1;
            performanceAllNumPosLARGE1  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPosLARGE1;
            performanceAllOutOfNumLARGE1(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE1;
        end
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
    title({sprintf('%s%s%s%d%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName,subjectName,'\_', sessionNumber, ...
        ': threshold0=',threshold(1),' threshold1=(',threshold(2),',',threshold(3),')', ...
        ' threshold2=(',threshold(4),',',threshold(5),')'),''});
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
    
    % Plot performance across all conditions as above, for NoiseLevel0 and NoiseLevel1.
    for nn = 1:2
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
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
        threshold(row) = thresholdthis;
        row = row+1;
    end
      
    % Plot performance across all conditions as above, for NoiseLevel2.        
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(3));
    
    % Calculate performance across all conditions, per noise amount combination level.
    performanceAllNumPosSMALL2   = nan(nComparisons,nConditions);
    performanceAllOutOfNumSMALL2 = nan(nComparisons,nConditions);
    performanceAllNumPosLARGE2   = nan(nComparisons,nConditions);
    performanceAllOutOfNumLARGE2 = nan(nComparisons,nConditions);
    for ii = 1:nConditions
        conditionName  = sprintf('%s%d','condition',ii);
        performanceAllNumPosSMALL2  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPosSMALL2;
        performanceAllOutOfNumSMALL2(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNumSMALL2;
        performanceAllNumPosLARGE2  (:,ii) = data.performance.(noiseLevelName).(conditionName).NumPosLARGE2;
        performanceAllOutOfNumLARGE2(:,ii) = data.performance.(noiseLevelName).(conditionName).OutOfNumLARGE2;
    end
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
    row = row+1;

    % Plot parameters.
    title({sprintf('%s%s%s%d%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s',experimentName,subjectName,'\_', sessionNumber, ...
        ': threshold0=',threshold(1),' threshold1=',threshold(2), ...
        ' threshold2=(',threshold(3),',',threshold(4),')'),''});
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
    title({sprintf('%s%s%s%d%s%d%s%d%s%d',experimentName,subjectName,'\_', sessionNumber, ...
        ': mean0=',meanRT(1),' mean1=',meanRT(2),' mean2=',meanRT(3)),''});
    legend('Noise0 data','Noise1 data','Noise2 data','Location','northwest')
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