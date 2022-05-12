function data = analyzeSummaryData(varargin)
%analyzeSummaryData
%
% Usage:
%   data = analyzeSummaryData('experimentName', 'Experiment100');
%
% Description:
%   Analyze summary psychophysical data for all participants for the 
%   specified experiment.
%
% Optional parameters/values:
%   'experimentName' : (char) Name of experiment folder (default: 'Experiment100')
%
% History:
%   05/11/22  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set path to data folder
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to data folder.
pathToFolder = fullfile(getpref(projectName,'BaseDirAnalysis'), ...
               experimentName,'PsychophysicalDataAnalysis');
    
%% Load participant list
%
% Set path to participant list document.
participantList = fullfile(pathToFolder,sprintf('%s_ParticipantList.rtf',experimentName));

% Open, read, close document.
fileID = fopen(participantList,'r');
tmp = textscan(fileID,'%s','Delimiter','\n');
list = tmp{1};
fclose(fileID);

% Get participant list.
idx = contains(list,'subject');
list = list(idx);
participants = cellfun(@(x){x(1:16)}, list);

%% Get experiment analysis results for each participant above
%
% Save each participant's results in the struct below.
data = struct;

% Get each participant's results.
for ii = 1:numel(participants)
    participant = participants{ii};
    resultsFile = fullfile(pathToFolder,participant,sprintf('experimentAnalysis%s.mat',participant(8:end)));
    temp = load(resultsFile,'dataExperiment'); dataExperiment = temp.dataExperiment; clear temp;
    fn = fieldnames(dataExperiment);
    for jj=1:numel(fn)
        data(ii).(fn{jj}) = dataExperiment.(fn{jj});
    end
end

%% Plot performance for each participant and save analysis results
%
% Plot colors for each noise level.
colors{1}='k'; colors{2}=[255 165 0]/255; colors{3}='r';

% For each participant, plot performance across all sessions, per noise level.
for ii = 1:numel(data)
    
    % Get experiment parameters to plot.
    noiseLevels = data(ii).noiseLevels;
    comparisons = data(ii).comparisons;
    nNoiseLevels = numel(noiseLevels);
    
    % Convert comparison amounts from mm to degrees of visual angle.
    monitorDistance   = 1.2; % meters (distance of the scene in iset3d)
    monitorDistancemm = monitorDistance*1000; % convert to mm
    comparisonsDeg    = atand(comparisons/monitorDistancemm);
    
    % Save analysis results for this participant.
    data(ii).results.comparisonsDeg = comparisonsDeg;
        
    % Plot each noise level.
    threshold = nan(nNoiseLevels,1);
    pse       = nan(nNoiseLevels,1);
    
    % Plot individual participant figure.
    figure; hold on;
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % For this individual participant, average performance across all sessions.
        performanceAllNumPos   = data(ii).performance.(noiseLevelName).NumPos;
        performanceAllOutOfNum = data(ii).performance.(noiseLevelName).OutOfNum;
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
        
        % Save performance analysis results for this participant.
        data(ii).results.(noiseLevelName).NumPos      = NumPos;
        data(ii).results.(noiseLevelName).OutOfNum    = OutOfNum;
        data(ii).results.(noiseLevelName).performance = performanceAll;
        data(ii).results.(noiseLevelName).threshold   = thresholdthis;
        data(ii).results.(noiseLevelName).pse         = psethis; 
    end
    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f',experimentName,data(ii).subjectName, ...
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

%% Plot summary plot of above for all participants: psychometric function
%
% NOTE: for Experiment100, the same noise levels and comparison offsets
% were used for all participants.
noiseLevels    = data(1).noiseLevels;
nNoiseLevels   = numel(noiseLevels);
comparisonsDeg = data(1).results.comparisonsDeg;

% Calculate average performance across all participants.
performanceAll = nan(numel(data),numel(comparisonsDeg),nNoiseLevels);
numPosAll      = nan(numel(data),numel(comparisonsDeg),nNoiseLevels);
outOfNumAll    = nan(numel(data),numel(comparisonsDeg),nNoiseLevels);
for ii = 1:numel(data)
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        performanceAll(ii,:,nn) = data(ii).results.(noiseLevelName).performance;
        numPosAll     (ii,:,nn) = data(ii).results.(noiseLevelName).NumPos;
        outOfNumAll   (ii,:,nn) = data(ii).results.(noiseLevelName).OutOfNum;
    end  
end
performanceAvg = nanmean(performanceAll,1);
numPosSum      = nansum(numPosAll,1);
outOfNumSum    = nansum(outOfNumAll,1);

% Plot average performance across all participants.
threshold = nan(nNoiseLevels,1);
pse       = nan(nNoiseLevels,1);
figure; hold on;
for nn = 1:nNoiseLevels
    
    % Plot data.
    plot(comparisonsDeg,performanceAvg(1,:,nn),'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});

    % Plot psychometric function fit.
    [xx,FittedCurve,thresholdthis,psethis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,nn)',outOfNumSum(1,:,nn)');
    plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
    threshold(nn) = thresholdthis;
    pse(nn) = psethis;
end

% Plot parameters.
title({sprintf('%s%s%0.2f%s%0.2f%s%0.2f',experimentName, ...
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

%% Plot summary plots of above for all participants: x-y scatterplots of thresholds
%
% Get all participants' thresholds, per noise level.
thresholdAll = nan(numel(data),nNoiseLevels);
for ii = 1:numel(data)
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        thresholdAll(ii,nn) = data(ii).results.(noiseLevelName).threshold;
    end
end

% Run paired t-tests.
pvalues = nan(3,1);
[~,pvalues(1)] = ttest(thresholdAll(:,1),thresholdAll(:,2));
[~,pvalues(2)] = ttest(thresholdAll(:,1),thresholdAll(:,3));
[~,pvalues(3)] = ttest(thresholdAll(:,2),thresholdAll(:,3));

% Plot summary x-y scatterplot: NoiseLevel0 vs. NoiseLevel1
figure; hold on; axis square;
plot(thresholdAll(:,1),thresholdAll(:,2),'ok');
title(sprintf('%s: NoiseLevel0 vs. NoiseLevel1, p=%0.3f',experimentName,pvalues(1)));
xlabel('Noise 0 threshold');
ylabel('Noise 1 threshold');
axis([0 1 0 1]);
plot([0 1], [0 1], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 .2 .4 .6 .8 1]);
set(gca,'YTick',[0 .2 .4 .6 .8 1]);
xtickformat('%.1f');
ytickformat('%.1f');
box off; hold off;

% Plot summary x-y scatterplot: NoiseLevel0 vs. NoiseLevel2
figure; hold on; axis square;
plot(thresholdAll(:,1),thresholdAll(:,3),'ok');
title(sprintf('%s: NoiseLevel0 vs. NoiseLevel2, p=%0.3f',experimentName,pvalues(2)));
xlabel('Noise 0 threshold');
ylabel('Noise 2 threshold');
axis([0 1 0 1]);
plot([0 1], [0 1], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 .2 .4 .6 .8 1]);
set(gca,'YTick',[0 .2 .4 .6 .8 1]);
xtickformat('%.1f');
ytickformat('%.1f');
box off; hold off;

% Plot summary x-y scatterplot: NoiseLevel1 vs. NoiseLevel2
figure; hold on; axis square;
plot(thresholdAll(:,2),thresholdAll(:,3),'ok');
title(sprintf('%s: NoiseLevel1 vs. NoiseLevel2, p=%0.3f',experimentName,pvalues(3)));
xlabel('Noise 1 threshold');
ylabel('Noise 2 threshold');
axis([0 1 0 1]);
plot([0 1], [0 1], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 .2 .4 .6 .8 1]);
set(gca,'YTick',[0 .2 .4 .6 .8 1]);
xtickformat('%.1f');
ytickformat('%.1f');
box off; hold off;

%% Plot and save performance for each participant, separated by NoiseAmounts1 combination level (smaller, larger)
for ii = 1:numel(data)
    
    % Get experiment parameters to plot.
    noiseLevels = data(ii).noiseLevels;
    nNoiseLevels = numel(noiseLevels);
    comparisonsDeg = data(ii).results.comparisonsDeg;
    
    % Plot all noise levels.
    threshold = nan(nNoiseLevels+nNoiseLevels-1,1);
    row = 1;
    figure; hold on;
    
    % Average performance across all sessions, for the first noise level (NoiseLevel0).
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(1));
    performanceAllNumPos   = data(ii).performance.(noiseLevelName).NumPos;
    performanceAllOutOfNum = data(ii).performance.(noiseLevelName).OutOfNum;
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
        performanceAllNumPosSMALL1   = data(ii).performance.(noiseLevelName).NumPosSMALL1;
        performanceAllOutOfNumSMALL1 = data(ii).performance.(noiseLevelName).OutOfNumSMALL1;
        performanceAllNumPosLARGE1   = data(ii).performance.(noiseLevelName).NumPosLARGE1;
        performanceAllOutOfNumLARGE1 = data(ii).performance.(noiseLevelName).OutOfNumLARGE1;
        NumPosSMALL   = sum(performanceAllNumPosSMALL1,  2);
        OutOfNumSMALL = sum(performanceAllOutOfNumSMALL1,2);
        performanceAllSMALL = NumPosSMALL./OutOfNumSMALL;
        NumPosLARGE   = sum(performanceAllNumPosLARGE1,  2);
        OutOfNumLARGE = sum(performanceAllOutOfNumLARGE1,2);
        performanceAllLARGE = NumPosLARGE./OutOfNumLARGE;
        
        % Plot data: SMALL noise amount combination level.
        plot(comparisonsDeg,performanceAllSMALL,'o','MarkerFace','w','MarkerEdge',colors{nn});
        % Plot psychometric function fit: SMALL noise amount combination level.
        [xx,FittedCurve,thresholdSMALL] = fitPsychometric(comparisonsDeg,NumPosSMALL,OutOfNumSMALL);
        plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdSMALL;
        row = row+1;
        
        % Plot data: LARGE noise amount combination level.
        plot(comparisonsDeg,performanceAllLARGE,'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        % Plot psychometric function fit: LARGE noise amount combination level.
        [xx,FittedCurve,thresholdLARGE] = fitPsychometric(comparisonsDeg,NumPosLARGE,OutOfNumLARGE);
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdLARGE;
        row = row+1;
        
        % Save performance analysis results for this participant.
        data(ii).results.(noiseLevelName).Noise1combosSMALLnumPos      = NumPosSMALL;
        data(ii).results.(noiseLevelName).Noise1combosSMALLoutOfNum    = OutOfNumSMALL;
        data(ii).results.(noiseLevelName).Noise1combosSMALLperformance = performanceAllSMALL;
        data(ii).results.(noiseLevelName).Noise1combosSMALLthreshold   = thresholdSMALL;
        data(ii).results.(noiseLevelName).Noise1combosLARGEnumPos      = NumPosLARGE;
        data(ii).results.(noiseLevelName).Noise1combosLARGEoutOfNum    = OutOfNumLARGE;
        data(ii).results.(noiseLevelName).Noise1combosLARGEperformance = performanceAllLARGE;
        data(ii).results.(noiseLevelName).Noise1combosLARGEthreshold   = thresholdLARGE;
    end
    
    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName,data(ii).subjectName, ...
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

%% Plot summary plot of above for all participants: psychometric function
%
% NOTE: for Experiment100, the same noise levels and comparison offsets
% were used for all participants.
noiseLevels    = data(1).noiseLevels;
nNoiseLevels   = numel(noiseLevels);
comparisonsDeg = data(1).results.comparisonsDeg;

% Calculate average performance across all participants.
performanceAll = nan(numel(data),numel(comparisonsDeg),nNoiseLevels+nNoiseLevels-1);
numPosAll      = nan(numel(data),numel(comparisonsDeg),nNoiseLevels+nNoiseLevels-1);
outOfNumAll    = nan(numel(data),numel(comparisonsDeg),nNoiseLevels+nNoiseLevels-1);
row = 2;

% Calculate average performance across all participants.
for nn = 1:nNoiseLevels
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
    for ii = 1:numel(data)
        if nn==1
            performanceAll(ii,:,1) = data(ii).results.(noiseLevelName).performance;
            numPosAll     (ii,:,1) = data(ii).results.(noiseLevelName).NumPos;
            outOfNumAll   (ii,:,1) = data(ii).results.(noiseLevelName).OutOfNum;
        else
            % Smaller NoiseAmounts1 combination levels.
            performanceAll(ii,:,row) = data(ii).results.(noiseLevelName).Noise1combosSMALLperformance;
            numPosAll     (ii,:,row) = data(ii).results.(noiseLevelName).Noise1combosSMALLnumPos;
            outOfNumAll   (ii,:,row) = data(ii).results.(noiseLevelName).Noise1combosSMALLoutOfNum;
            % Larger NoiseAmounts1 combination levels.
            performanceAll(ii,:,row+1) = data(ii).results.(noiseLevelName).Noise1combosLARGEperformance;
            numPosAll     (ii,:,row+1) = data(ii).results.(noiseLevelName).Noise1combosLARGEnumPos;
            outOfNumAll   (ii,:,row+1) = data(ii).results.(noiseLevelName).Noise1combosLARGEoutOfNum; 
        end
    end
    if nn==2
        row = row+2;
    end
end
performanceAvg = nanmean(performanceAll,1);
numPosSum      = nansum(numPosAll,1);
outOfNumSum    = nansum(outOfNumAll,1);

% Plot average performance across all participants.
threshold = nan(nNoiseLevels+nNoiseLevels-1,1);
figure; hold on;
row = 2;
for nn = 1:nNoiseLevels
    if nn==1
        
        % Plot data and psychometric function fit.
        plot(comparisonsDeg,performanceAvg(1,:,nn),'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,nn)',outOfNumSum(1,:,nn)');
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(1) = thresholdthis;
    else
        
        % Plot data and psychometric function fit: smaller NoiseAmounts1 combination levels.
        plot(comparisonsDeg,performanceAvg(1,:,row),'o','MarkerFace','w','MarkerEdge',colors{nn});
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,row)',outOfNumSum(1,:,row)');
        plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdthis;
        row = row+1;
        
        % Plot data and psychometric function fit: larger NoiseAmounts1 combination levels.
        plot(comparisonsDeg,performanceAvg(1,:,row),'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,row)',outOfNumSum(1,:,row)');
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(row) = thresholdthis;
        row = row+1;
    end
end

% Plot parameters.
title({sprintf('%s%s%0.2f%s%0.2f%s%0.2f%s%s%0.2f%s%0.2f%s',experimentName, ...
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

%% Plot summary plots of above for all participants: x-y scatterplots of thresholds
%
% Get all participants' thresholds, per noise level and noise amount combo level.
thresholdAll = nan(numel(data),nNoiseLevels+nNoiseLevels-1);
row = 2;
for nn = 1:nNoiseLevels
    for ii = 1:numel(data)
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Get thresholds for NoiseLevel0.
        if nn==1
            thresholdAll(ii,1) = data(ii).results.(noiseLevelName).threshold;
        else
            
            %Get thresholds for smaller NoiseAmounts1 combination levels.
            thresholdAll(ii,row) = data(ii).results.(noiseLevelName).Noise1combosSMALLthreshold;
   
            %Get thresholds for larger NoiseAmounts1 combination levels.
            thresholdAll(ii,row+1) = data(ii).results.(noiseLevelName).Noise1combosLARGEthreshold;
        end  
    end
    if nn==2
        row = row+2;
    end
end

% Run paired t-tests.
pvalues = nan(2,1);
[~,pvalues(1)] = ttest(thresholdAll(:,2),thresholdAll(:,3));
[~,pvalues(2)] = ttest(thresholdAll(:,4),thresholdAll(:,5));

% Plot summary x-y scatterplot: NoiseLevel1, smaller vs. larger NoiseAmounts1 combos.
figure; hold on; axis square;
plot(thresholdAll(:,2),thresholdAll(:,3),'ok');
title(sprintf('%s: NoiseLevel1, smaller vs. larger NoiseAmounts1 combos, p=%0.3f',experimentName,pvalues(1)));
xlabel('Smaller combos threshold');
ylabel('Larger combos threshold');
axis([0 1 0 1]);
plot([0 1], [0 1], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 .2 .4 .6 .8 1]);
set(gca,'YTick',[0 .2 .4 .6 .8 1]);
xtickformat('%.1f');
ytickformat('%.1f');
box off; hold off;

% Plot summary x-y scatterplot: NoiseLevel2, smaller vs. larger NoiseAmounts1 combos.
figure; hold on; axis square;
plot(thresholdAll(:,4),thresholdAll(:,5),'ok');
title(sprintf('%s: NoiseLevel2, smaller vs. larger NoiseAmounts1 combos, p=%0.3f',experimentName,pvalues(2)));
xlabel('Smaller combos threshold');
ylabel('Larger combos threshold');
axis([0 1 0 1]);
plot([0 1], [0 1], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 .2 .4 .6 .8 1]);
set(gca,'YTick',[0 .2 .4 .6 .8 1]);
xtickformat('%.1f');
ytickformat('%.1f');
box off; hold off;

%% Plot and save performance for each participant, separated by NoiseAmounts2 combination level (smaller, larger)
for ii = 1:numel(data)
    
    % Get experiment parameters to plot.
    noiseLevels = data(ii).noiseLevels;
    nNoiseLevels = numel(noiseLevels);
    comparisonsDeg = data(ii).results.comparisonsDeg;
    
    % Plot all noise levels.
    threshold = nan(nNoiseLevels+nNoiseLevels-2,1);
    row = 1;
    figure; hold on;
    
    % Average performance across all sessions, for NoiseLevel0 and NoiseLevel1.
    for nn = 1:2
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        performanceAllNumPos   = data(ii).performance.(noiseLevelName).NumPos;
        performanceAllOutOfNum = data(ii).performance.(noiseLevelName).OutOfNum;
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
    performanceAllNumPosSMALL2   = data(ii).performance.(noiseLevelName).NumPosSMALL2;
    performanceAllOutOfNumSMALL2 = data(ii).performance.(noiseLevelName).OutOfNumSMALL2;
    performanceAllNumPosLARGE2   = data(ii).performance.(noiseLevelName).NumPosLARGE2;
    performanceAllOutOfNumLARGE2 = data(ii).performance.(noiseLevelName).OutOfNumLARGE2;
    NumPosSMALL   = sum(performanceAllNumPosSMALL2,  2);
    OutOfNumSMALL = sum(performanceAllOutOfNumSMALL2,2);
    performanceAllSMALL = NumPosSMALL./OutOfNumSMALL;
    NumPosLARGE   = sum(performanceAllNumPosLARGE2,  2);
    OutOfNumLARGE = sum(performanceAllOutOfNumLARGE2,2);
    performanceAllLARGE = NumPosLARGE./OutOfNumLARGE;
    
    % Plot data: SMALL noise amount combination level.
    plot(comparisonsDeg,performanceAllSMALL,'o','MarkerFace','w','MarkerEdge',colors{3});
    % Plot psychometric function fit: SMALL noise amount combination level.
    [xx,FittedCurve,thresholdSMALL] = fitPsychometric(comparisonsDeg,NumPosSMALL,OutOfNumSMALL);
    plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{3});
    threshold(row) = thresholdSMALL;
    row = row+1;
    
    % Plot data: LARGE noise amount combination level.
    plot(comparisonsDeg,performanceAllLARGE,'o','MarkerFace',colors{3},'MarkerEdge',colors{3});
    % Plot psychometric function fit: LARGE noise amount combination level.
    [xx,FittedCurve,thresholdLARGE] = fitPsychometric(comparisonsDeg,NumPosLARGE,OutOfNumLARGE);
    plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{3});
    threshold(row) = thresholdLARGE;
    
    % Save performance analysis results for this participant.
    data(ii).results.(noiseLevelName).Noise2combosSMALLnumPos      = NumPosSMALL;
    data(ii).results.(noiseLevelName).Noise2combosSMALLoutOfNum    = OutOfNumSMALL;
    data(ii).results.(noiseLevelName).Noise2combosSMALLperformance = performanceAllSMALL;
    data(ii).results.(noiseLevelName).Noise2combosSMALLthreshold   = thresholdSMALL;
    data(ii).results.(noiseLevelName).Noise2combosLARGEnumPos      = NumPosLARGE;
    data(ii).results.(noiseLevelName).Noise2combosLARGEoutOfNum    = OutOfNumLARGE;
    data(ii).results.(noiseLevelName).Noise2combosLARGEperformance = performanceAllLARGE;
    data(ii).results.(noiseLevelName).Noise2combosLARGEthreshold   = thresholdLARGE;
    
    % Plot parameters.
    title({sprintf('%s_%s%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s',experimentName,data(ii).subjectName, ...
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

%% Plot summary plot of above for all participants: psychometric function
%
% NOTE: for Experiment100, the same noise levels and comparison offsets
% were used for all participants.
noiseLevels    = data(1).noiseLevels;
nNoiseLevels   = numel(noiseLevels);
comparisonsDeg = data(1).results.comparisonsDeg;

% Calculate average performance across all participants.
performanceAll = nan(numel(data),numel(comparisonsDeg),nNoiseLevels+nNoiseLevels-2);
numPosAll      = nan(numel(data),numel(comparisonsDeg),nNoiseLevels+nNoiseLevels-2);
outOfNumAll    = nan(numel(data),numel(comparisonsDeg),nNoiseLevels+nNoiseLevels-2);

% Calculate average performance across all participants.
for nn = 1:nNoiseLevels
    noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
    for ii = 1:numel(data)
        if nn==1 || nn==2
            performanceAll(ii,:,nn) = data(ii).results.(noiseLevelName).performance;
            numPosAll     (ii,:,nn) = data(ii).results.(noiseLevelName).NumPos;
            outOfNumAll   (ii,:,nn) = data(ii).results.(noiseLevelName).OutOfNum;
        else
            % Smaller NoiseAmounts2 combination levels.
            performanceAll(ii,:,nn) = data(ii).results.(noiseLevelName).Noise2combosSMALLperformance;
            numPosAll     (ii,:,nn) = data(ii).results.(noiseLevelName).Noise2combosSMALLnumPos;
            outOfNumAll   (ii,:,nn) = data(ii).results.(noiseLevelName).Noise2combosSMALLoutOfNum;
            % Larger NoiseAmounts2 combination levels.
            performanceAll(ii,:,nn+1) = data(ii).results.(noiseLevelName).Noise2combosLARGEperformance;
            numPosAll     (ii,:,nn+1) = data(ii).results.(noiseLevelName).Noise2combosLARGEnumPos;
            outOfNumAll   (ii,:,nn+1) = data(ii).results.(noiseLevelName).Noise2combosLARGEoutOfNum; 
        end
    end
end
performanceAvg = nanmean(performanceAll,1);
numPosSum      = nansum(numPosAll,1);
outOfNumSum    = nansum(outOfNumAll,1);

% Plot average performance across all participants.
threshold = nan(nNoiseLevels+nNoiseLevels-2,1);
figure; hold on;
for nn = 1:nNoiseLevels
    if nn==1 || nn==2
        
        % Plot data and psychometric function fit.
        plot(comparisonsDeg,performanceAvg(1,:,nn),'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,nn)',outOfNumSum(1,:,nn)');
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(nn) = thresholdthis;
    else
        
        % Plot data and psychometric function fit: smaller NoiseAmounts2 combination levels.
        plot(comparisonsDeg,performanceAvg(1,:,nn),'o','MarkerFace','w','MarkerEdge',colors{nn});
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,nn)',outOfNumSum(1,:,nn)');
        plot(xx,FittedCurve,'--','LineWidth',1,'Color',colors{nn});
        threshold(nn) = thresholdthis;
        
        % Plot data and psychometric function fit: larger NoiseAmounts2 combination levels.
        plot(comparisonsDeg,performanceAvg(1,:,nn+1),'o','MarkerFace',colors{nn},'MarkerEdge',colors{nn});
        [xx,FittedCurve,thresholdthis] = fitPsychometric(comparisonsDeg,numPosSum(1,:,nn+1)',outOfNumSum(1,:,nn+1)');
        plot(xx,FittedCurve,'-','LineWidth',1,'Color',colors{nn});
        threshold(nn+1) = thresholdthis;
    end
end

% Plot parameters.
title({sprintf('%s%s%0.2f%s%0.2f%s%0.2f%s%0.2f%s',experimentName, ...
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

%% Plot summary plots of above for all participants: x-y scatterplots of thresholds
%
% Get all participants' thresholds, per noise level and noise amount combo level.
thresholdAll = nan(numel(data),nNoiseLevels+nNoiseLevels-2);
for nn = 1:nNoiseLevels
    for ii = 1:numel(data)
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Get thresholds for NoiseLevel0 and NoiseLevel1.
        if nn==1 || nn==2
            thresholdAll(ii,nn) = data(ii).results.(noiseLevelName).threshold;
        else
            
            %Get thresholds for smaller NoiseAmounts2 combination levels.
            thresholdAll(ii,nn) = data(ii).results.(noiseLevelName).Noise2combosSMALLthreshold;
   
            %Get thresholds for larger NoiseAmounts2 combination levels.
            thresholdAll(ii,nn+1) = data(ii).results.(noiseLevelName).Noise2combosLARGEthreshold;
        end  
    end
end

% Run paired t-tests.
pvalues = nan;
[~,pvalues(1)] = ttest(thresholdAll(:,3),thresholdAll(:,4));

% Plot summary x-y scatterplot: NoiseLevel2, smaller vs. larger NoiseAmounts2 combos.
figure; hold on; axis square;
plot(thresholdAll(:,3),thresholdAll(:,4),'ok');
title(sprintf('%s: NoiseLevel2, smaller vs. larger NoiseAmounts2 combos, p=%0.3f',experimentName,pvalues(1)));
xlabel('Smaller combos threshold');
ylabel('Larger combos threshold');
axis([0 1 0 1]);
plot([0 1], [0 1], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 .2 .4 .6 .8 1]);
set(gca,'YTick',[0 .2 .4 .6 .8 1]);
xtickformat('%.1f');
ytickformat('%.1f');
box off; hold off;

%% Plot and save reaction times for each participant
for ii = 1:numel(data)
    
    % Get experiment parameters to plot.
    noiseLevels = data(ii).noiseLevels;
    comparisons = data(ii).comparisons;
    nNoiseLevels = numel(noiseLevels);
    
    % Convert comparison amounts from mm to degrees of visual angle.
    monitorDistance   = 1.2; % meters (distance of the scene in iset3d)
    monitorDistancemm = monitorDistance*1000; % convert to mm
    comparisonsDeg    = atand(comparisons/monitorDistancemm);
    
    % Plot each noise level.
    meanRT = nan(nNoiseLevels,1);
    
    % Plot individual participant figure.
    figure; hold on;
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        
        % Average reaction time across all sessions.
        reactionTime    = data(ii).performance.(noiseLevelName).reactionTime;
        reactionTimeAll = nanmean(reactionTime,2);
        
        % Plot data.
        plot(comparisonsDeg,reactionTimeAll,'Color',colors{nn});
        meanRT(nn) = round(nanmean(reactionTimeAll));
        
        % Save performance analysis results for this participant.
        data(ii).results.(noiseLevelName).reactionTime = reactionTimeAll;
    end
    
    % Plot parameters.
    title({sprintf('%s_%s%s%d%s%d%s%d',experimentName,data(ii).subjectName, ...
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

%% Plot summary plot of above for all participants: average across participants
%
% NOTE: for Experiment100, the same noise levels and comparison offsets
% were used for all participants.
noiseLevels    = data(1).noiseLevels;
nNoiseLevels   = numel(noiseLevels);
comparisonsDeg = data(1).results.comparisonsDeg;

% Calculate average reaction times across all participants.
reactionTimeAll = nan(numel(data),numel(comparisonsDeg),nNoiseLevels);
for ii = 1:numel(data)
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        reactionTimeAll(ii,:,nn) = data(ii).results.(noiseLevelName).reactionTime;
    end  
end
reactionTimeAvg = nanmean(reactionTimeAll,1);

% Plot average reaction times across all participants.
meanRT = nan(nNoiseLevels,1);
figure; hold on;
for nn = 1:nNoiseLevels
    plot(comparisonsDeg,reactionTimeAvg(1,:,nn),'Color',colors{nn});
    meanRT(nn) = round(nanmean(reactionTimeAvg(1,:,nn)));
end

% Plot parameters.
title({sprintf('%s%s%d%s%d%s%d',experimentName, ...
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

%% Plot summary plots of above for all participants: x-y scatterplots of mean reaction times
%
% Get all participants' mean reaction time, per noise level.
reactionTimeAll = nan(numel(data),nNoiseLevels);
for ii = 1:numel(data)
    for nn = 1:nNoiseLevels
        noiseLevelName = sprintf('%s%d','noiseLevel',noiseLevels(nn));
        reactionTimeAll(ii,nn) = nanmean(data(ii).results.(noiseLevelName).reactionTime);
    end
end

% Run paired t-tests.
pvalues = nan(3,1);
[~,pvalues(1)] = ttest(reactionTimeAll(:,1),reactionTimeAll(:,2));
[~,pvalues(2)] = ttest(reactionTimeAll(:,1),reactionTimeAll(:,3));
[~,pvalues(3)] = ttest(reactionTimeAll(:,2),reactionTimeAll(:,3));

% Plot summary x-y scatterplot: NoiseLevel0 vs. NoiseLevel1
figure; hold on; axis square;
plot(reactionTimeAll(:,1),reactionTimeAll(:,2),'ok');
title(sprintf('%s: NoiseLevel0 vs. NoiseLevel1, p=%0.3f',experimentName,pvalues(1)));
xlabel('Noise 0 mean reaction time (ms)');
ylabel('Noise 1 mean reaction time (ms)');
axis([0 700 0 700]);
plot([0 700], [0 700], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 200 400 600]);
set(gca,'YTick',[0 200 400 600]);
box off; hold off;

% Plot summary x-y scatterplot: NoiseLevel0 vs. NoiseLevel2
figure; hold on; axis square;
plot(reactionTimeAll(:,1),reactionTimeAll(:,3),'ok');
title(sprintf('%s: NoiseLevel0 vs. NoiseLevel2, p=%0.3f',experimentName,pvalues(2)));
xlabel('Noise 0 mean reaction time (ms)');
ylabel('Noise 2 mean reaction time (ms)');
axis([0 700 0 700]);
plot([0 700], [0 700], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 200 400 600]);
set(gca,'YTick',[0 200 400 600]);
box off; hold off;

% Plot summary x-y scatterplot: NoiseLevel1 vs. NoiseLevel2
figure; hold on; axis square;
plot(reactionTimeAll(:,2),reactionTimeAll(:,3),'ok');
title(sprintf('%s: NoiseLevel1 vs. NoiseLevel2, p=%0.3f',experimentName,pvalues(3)));
xlabel('Noise 1 mean reaction time (ms)');
ylabel('Noise 2 mean reaction time (ms)');
axis([0 700 0 700]);
plot([0 700], [0 700], '-k');
set(gca,'tickdir','out');
set(gca,'XTick',[0 200 400 600]);
set(gca,'YTick',[0 200 400 600]);
box off; hold off;

end
%% End