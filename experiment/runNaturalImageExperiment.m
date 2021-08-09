function acquisitionStatus = runNaturalImageExperiment(varargin)
%runNaturalImageExperiment
%
% Usage:
%   runNaturalImageExperiment('experimentName','Experiment100','subjectName','test');
%
% Description:
%   Run the natural image psychophysics experiment given the specified
%   experiment folder. Save the experimental parameters and psychophysical 
%   responses (struct 'data') in the specified output folder.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment100')
%   'subjectName'    : (string)  Name of subject (default: 'test')
%   'nIterations'    : (scalar)  Number of iterations per image comparison *must be an even number (default: 14)
%   'controlSignal'  : (string)  Input method for user response (options: 'gamePad', 'keyboard') (default: 'gamePad')
%   'option1Key'     : (string)  For gamePad either 'GP:UpperLeftTrigger'  or 'GP:X', for keyboard -> '1' (default: 'GP:UpperLeftTrigger')
%   'option2Key'     : (string)  For gamePad either 'GP:UpperRightTrigger' or 'GP:A', for keyboard -> '2' (default: 'GP:UpperRightTrigger')
%   'giveFeedback'   : (logical) Give feedback if option is on (default: true)
%   'isDemo'         : (logical) Data won't be saved if on (default: false)
%   'isFirstSession' : (logical) Warmup trials will be run if on (default: true)
%
% History:
%   06/07/21  amn  Adapted from BrainardLab/VirtualWorldPsychophysics
%   07/05/21  amn  Edits based on pilot tests
%   07/19/21  amn  Edits for main experiment 
%   08/03/21  amn  Final edits

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.addParameter('subjectName', 'test', @ischar);
parser.addParameter('nIterations', 14, @isscalar);
parser.addParameter('controlSignal', 'gamePad', @ischar);
parser.addParameter('option1Key', 'GP:UpperLeftTrigger', @ischar);
parser.addParameter('option2Key', 'GP:UpperRightTrigger', @ischar);
parser.addParameter('giveFeedback', true, @islogical);
parser.addParameter('isDemo', false, @islogical);
parser.addParameter('isFirstSession', true, @islogical);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
subjectName    = parser.Results.subjectName;
nIterations    = parser.Results.nIterations;
controlSignal  = parser.Results.controlSignal;
option1Key     = parser.Results.option1Key;
option2Key     = parser.Results.option2Key;
giveFeedback   = parser.Results.giveFeedback;
isDemo         = parser.Results.isDemo;
isFirstSession = parser.Results.isFirstSession;

% Give warning that experiment will not run if 'nIterations' is an odd number.
if rem(nIterations,2)==1
    error('nIterations must be an even number.');
end

%% Get user input to verify/change inputs
%
% Ask user if the value of 'experimentName' is correct.
fprintf(2,'The experiment name is: %s\n', experimentName);
str1 = input('Is this correct? Enter Y if Yes, N if No: ','s');
if strcmpi(str1,'N')
    fprintf(2,'Enter the correct experiment name below\n');
    experimentName = input('Enter here (without quotes): ','s');
    fprintf(2,'The experiment name has been updated to: %s\n', experimentName);
end

% Ask user if the value of 'subjectName' is correct.
fprintf(2,'The subject name is: %s\n', subjectName);
str1 = input('Is this correct? Enter Y if Yes, N if No: ','s');
if strcmpi(str1,'N')
    fprintf(2,'Enter the correct subject name below\n');
    subjectName = input('Enter here (without quotes): ','s');
    fprintf(2,'The subject name has been updated to: %s\n', subjectName);
end

% Ask user to input 'isFirstSession' value.
str1 = input('Is this the participant''s first session? Enter Y if Yes, N if No: ','s');
if strcmpi(str1,'Y')
    isFirstSession = true;
    fprintf(2,'Because this is the first session, warmup trials will be run.\n');
elseif strcmpi(str1,'N')
    isFirstSession = false;
    fprintf(2,'Warmup trials will not be run.\n');
end

%% Set paths to folders
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to image folder.
if strcmpi(subjectName,'ANauthor')
    pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBsAmy');
else
    pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBs');
end

% Set path to output folder (only experimental computers have Write access
% to this data folder).
pathToOutput = fullfile(getpref(projectName,'BaseDirData'),experimentName,'PsychophysicalData');

%% Set up experiment parameters
%
% Specify feedback sounds.
rightSound = sin(2*pi*(1:1000)/10)/10;
wrongSound = rand(1,1000).*ceil(sin(2*pi*(1:1000)/10))/10;

% Set acquisition status.
acquisitionStatus = 0;

% Set task parameters.
params.screenDimsCm = [59.67 33.57]; % cm
params.bgColor      = [128 128 128]/255; % to match electrophysiology task
params.textColor    = [0.6 0.2 0.2];
params.image1Loc  = [0 0];
params.image2Loc  = [0 0];
params.image1Size = [10.54 10.54]; % monitor distance=75cm: scene 8 deg vis angle (target 4 deg)
params.image2Size = [10.54 10.54];
params.ISI          = 0.40; % seconds
params.ITI          = 0.00; % seconds
params.stimDuration = 0.25; % seconds
params.option1Key = option1Key;
params.option2Key = option2Key;
params.nBlocks = 16; % number of blocks/image for mask

%% Get image info
%
% Get info about the images in the image folder.
fileInfo = dir([pathToFolder '/*.mat']);
    
% Get the number of images in the folder.
nImages = numel(fileInfo);

% Get the image file names. Images will be called by their index here.
imageNames = {fileInfo(:).name}';

%% EXPERIMENT ORGANIZATION
%
% CONDITION   : Per condition, 1 center position.
%               Per trial, a pseudorandom change position is compared to the center position.
%               5 (change left) + 5 (change right) + 1 (no change) = 11 comparisons per condition.
%               (range of comparisons defined by 'rangeMin' and 'rangeMax')
%
% COMPARISON  : Per comparison, center position presented randomly in either 1st/2nd interval.
%
% BLOCK       : Per block, trials will be interspersed from 2 conditions. 
%               11 comparisons per condition * 2 conditions = 22 total comparisons per block.
%
% ITERATION   : nIterations (default 14) of each block.
%               Run a complete block (22 total comparisons, pseudorandom order) 
%               before moving on to next block iteration.
%               22 total comparisons per block * 14 iterations of each block = 308 trials.
%
% NOISE LEVEL : 308 trials make up 1 noise level.
%               Noise is in a task-irrelevant feature/object.
%               Each noise level is made up of a pool of various task-irrelevant change amounts.
%               Per stimulus (center position & comparison) take a random draw from the pool.
%
% SESSION     : Per session, 3 noise levels (0, 1, and 2).
%               308 trials per noise level * 3 noise levels = 924 trials per session.
%
% RUN         : Each noise level will be divided into 2 runs.
%               Each run is made up of nIterations/2 blocks.
%               3 noise levels * 2 runs = 6 runs per session.
%               The order of runs will be randomly assigned per session.
%               1 minute break between runs.
%
% SUBJECT     : 6 sessions per subject.

%% Create vectors with image info
%
% Preallocate vectors of noise level, noise amount, condition, and 
% comparison amount per image.
imageNoiseLevel   = nan(nImages,1);
imageNoiseAmount1 = nan(nImages,1);
imageNoiseAmount2 = nan(nImages,1);
imageCondition    = nan(nImages,1);
imageComparison   = nan(nImages,1);

% Break down image file names to get image info.
for ii = 1:nImages
    name = imageNames{ii};
    p    = strfind(name,'_');

    % Label the 'imageNames' indices by their condition, comparison amount,
    % noise level, and noise amount.
    s1 = name(1:p(1));
    imageCondition(ii)  = str2double(regexp(s1,'[+-]?\d*','Match'));
    s2 = name(p(1):p(2));
    imageComparison(ii) = str2double(regexp(s2,'[+-]?\d*','Match'));
    s3 = name(p(2):p(3));
    imageNoiseLevel(ii) = str2double(regexp(s3,'[+-]?\d*','Match'));
    if numel(p)==3
        s4 = name(p(3):end);
        s5 = '0';
    elseif numel(p)==4
        s4 = name(p(3):p(4));   
        s5 = name(p(4):end);
    end
    imageNoiseAmount1(ii) = str2double(regexp(s4,'[+-]?\d*','Match'));
    imageNoiseAmount2(ii) = str2double(regexp(s5,'[+-]?\d*','Match'));
end

% Get the identities and number of noise levels, conditions, and comparison 
% amounts (same amounts for each condition).
noiseLevels   = unique(imageNoiseLevel);
nNoiseLevels  = numel(noiseLevels);
conditions    = unique(imageCondition);
nConditions   = numel(conditions);
comparisons   = unique(imageComparison);
nComparisons  = numel(comparisons); 

%% Create a trial order for this session
%
% Preallocate matrix of trial order:
%   col 1: image index for 1st interval
%   col 2: image index for 2nd interval
nTrialsBlock = nComparisons*nConditions;
nTrialsRun   = nTrialsBlock*nIterations/2;
nTrials      = nTrialsRun*2*nNoiseLevels;
trialOrder   = nan(nTrials,2);
srow         = 1;

% Preallocate vectors of noise level and noise amounts (within the level) per trial.
trialNoiseLevel   = nan(nTrials,1);
trialNoiseAmount1 = nan(nTrials,2);
trialNoiseAmount2 = nan(nTrials,2);

% Create random order of runs (each noise level divided into two runs).
allRunsOrdered = [noiseLevels; noiseLevels];
allRuns        = allRunsOrdered(randperm(numel(allRunsOrdered)));

% Create a trial order per run, then combine across runs.
for jj = 1:numel(allRuns)
    noiseLevelthis       = allRuns(jj);
    trialOrderRun        = nan(nTrialsRun,2);
    trialNoiseAmount1Run = nan(nTrialsRun,2);
    trialNoiseAmount2Run = nan(nTrialsRun,2);
    rrow = 1;
    
    % Include nIterations/2 blocks per run.
    % Complete a block before moving on to the next block.
    % Within each block, create a randomized order of trials.
    for ii = 1:nIterations/2
        trialsPerBlock        = nan(nTrialsBlock,2);
        noiseAmounts1PerBlock = nan(nTrialsBlock,2);
        noiseAmounts2PerBlock = nan(nTrialsBlock,2);
        brow = 1;
        
        % Create an ordered list of image index comparisons for trials of a single block.
        for jjj = 1:nConditions
            
            % Get the center position for this condition.
            centerpos = conditions(jjj);
            
            % For the center position of this condition, create a pool of
            % images (of the various noise amounts for this noise level,
            % including the noise amount of 0 from Noise Level 0).
            centerPool = find(imageComparison==0 & ...
                imageCondition==centerpos & ...
                (imageNoiseLevel==noiseLevelthis | imageNoiseLevel==0));

            % Create ordered list of image index comparisons for trials for this condition.
            trialsPerCondition       = nan(nComparisons,2);
            noiseAmount1PerCondition = nan(nComparisons,2);
            noiseAmount2PerCondition = nan(nComparisons,2);
            for iii = 1:nComparisons
                
                % Randomly draw from the pool of images of the center position.
                centerthis = centerPool(randi(numel(centerPool)));
                
                % For each comparison amount of this condition, create a pool
                % of images (of the various noise amounts for this noise level,
                % including the noise amount of 0 from Noise Level 0).
                comparisonPool = find(imageComparison==comparisons(iii) & ...
                    imageCondition==centerpos & ...
                    (imageNoiseLevel==noiseLevelthis | imageNoiseLevel==0));
                
                % Randomly draw the comparison image index from the above pool.
                comparisonthis = comparisonPool(randi(numel(comparisonPool)));
                
                % Randomize whether the center position is shown in the 1st/2nd interval.
                thisIndices       = [centerthis comparisonthis];
                thisNoiseAmounts1 = imageNoiseAmount1(thisIndices)';
                thisNoiseAmounts2 = imageNoiseAmount2(thisIndices)';
                orderrand         = randperm(2);
                trialsPerCondition(iii,:)       = thisIndices(orderrand);
                noiseAmount1PerCondition(iii,:) = thisNoiseAmounts1(orderrand);
                noiseAmount2PerCondition(iii,:) = thisNoiseAmounts2(orderrand);
            end
            
            % Combine ordered trials across all conditions for a single block.
            trialsPerBlock       (brow:brow+nComparisons-1,:) = trialsPerCondition;
            noiseAmounts1PerBlock(brow:brow+nComparisons-1,:) = noiseAmount1PerCondition;
            noiseAmounts2PerBlock(brow:brow+nComparisons-1,:) = noiseAmount2PerCondition;
            brow = brow+nComparisons;
        end
        
        % Randomize trial order within this single block.
        trialsrand = randperm(nTrialsBlock);
        trialsPerBlockRand        = trialsPerBlock       (trialsrand,:);
        noiseAmounts1PerBlockRand = noiseAmounts1PerBlock(trialsrand,:);
        noiseAmounts2PerBlockRand = noiseAmounts2PerBlock(trialsrand,:);
        
        % Combine iterations of blocks for this run.
        trialOrderRun       (rrow:rrow+nTrialsBlock-1,:) = trialsPerBlockRand;
        trialNoiseAmount1Run(rrow:rrow+nTrialsBlock-1,:) = noiseAmounts1PerBlockRand;
        trialNoiseAmount2Run(rrow:rrow+nTrialsBlock-1,:) = noiseAmounts2PerBlockRand;
        rrow = rrow+nTrialsBlock;
    end
    
    % Combine trials across runs.
    trialOrder       (srow:srow+nTrialsRun-1,:) = trialOrderRun;
    trialNoiseLevel  (srow:srow+nTrialsRun-1,1) = noiseLevelthis;
    trialNoiseAmount1(srow:srow+nTrialsRun-1,:) = trialNoiseAmount1Run;
    trialNoiseAmount2(srow:srow+nTrialsRun-1,:) = trialNoiseAmount2Run;
    srow = srow+nTrialsRun;
end

%% List the easy trials
%
% The practice trials and the first third of the warmup trials will be
% selected randomly from this list of easy trials.
% The easy trials will all be from NoiseLevel0.

% Preallocate matrix of trial order.
nEasyTrials = nConditions*2*2;
easyTrialOrder = nan(nEasyTrials,2);
trow = 1;

% Create a list of image index comparisons for easy comparisons.
for jjj = 1:nConditions
    centerpos = conditions(jjj);
    centerIdx = find(imageComparison==0 & imageCondition==centerpos & ...
                imageNoiseLevel==0);
    compIdx1  = find(imageComparison==comparisons(1) & imageCondition==centerpos & ...
                imageNoiseLevel==0);
    easyTrialOrder(trow,:) = [centerIdx compIdx1];
    trow = trow+1;
    easyTrialOrder(trow,:) = [compIdx1 centerIdx];
    trow = trow+1;
    compIdx2  = find(imageComparison==comparisons(end) & imageCondition==centerpos & ...
                imageNoiseLevel==0);
    easyTrialOrder(trow,:) = [centerIdx compIdx2];
    trow = trow+1;
    easyTrialOrder(trow,:) = [compIdx2 centerIdx];
    trow = trow+1;
end

%% Create a trial order for the warmup trials
%
% If this is the participant's first session, warmup trials will be run.
% The first third of the warmup trials will be selected randomly from the
% easy trials, the second third will be selected randomly from medium
% difficulty trials, and the final third will be selected randomly from the
% full list of trials.
% The warmup trials will all be from NoiseLevel0.

if isFirstSession
    % Preallocate matrix of trial order.
    nWarmupTrials = 30; %must be divisible by 3
    warmupTrialOrder = nan(nWarmupTrials,2);
    trow = 1;
    
    % Start warmup trials with randomly selected easy trials.
    for ii = 1:nWarmupTrials/3
        warmupTrialOrder(trow,:) = easyTrialOrder(randi(size(easyTrialOrder,1)),:);
        trow = trow+1;
    end
    
    % Randomly select medium difficulty trials.
    medComp = comparisons([2:3 end-2:end-1]);
    for ii = 1:nWarmupTrials/3
        centerpos = conditions(randi(nConditions));
        centerIdx = find(imageComparison==0 & imageCondition==centerpos & ...
            imageNoiseLevel==0);
        compIdx   = find(imageComparison==medComp(randi(numel(medComp))) & imageCondition==centerpos & ...
            imageNoiseLevel==0);
        thisIndices = [centerIdx compIdx];
        warmupTrialOrder(trow,:) = thisIndices(randperm(2));
        trow = trow+1;
    end
    
    % Randomly select trials from the full list of trials.
    for ii = 1:nWarmupTrials/3
        centerpos = conditions(randi(nConditions));
        centerIdx = find(imageComparison==0 & imageCondition==centerpos & ...
            imageNoiseLevel==0);
        compIdx   = find(imageComparison==comparisons(randi(nComparisons)) & imageCondition==centerpos & ...
            imageNoiseLevel==0);
        thisIndices = [centerIdx compIdx];
        warmupTrialOrder(trow,:) = thisIndices(randperm(2));
        trow = trow+1;
    end
end

%% Calculate the correct response for each trial
%
% Get comparison amount per image index.
trialOrderComparison = imageComparison(trialOrder);

% The correct response is 1 if the 2nd target is to the left of the 1st target.
% The correct response is 2 if the 2nd target is to the right of the 1st target.
trialDiff = diff(trialOrderComparison,1,2);
correctResponse = nan(nTrials,1);
correctResponse(trialDiff <0) = 1; % 2nd target is to the left
correctResponse(trialDiff >0) = 2; % 2nd target is to the right
correctResponse(trialDiff==0) = randi(2); % no difference: random assignment

%% Calculate the correct response for each easy trial
easyTrialOrderComparison = imageComparison(easyTrialOrder);
easyTrialDiff = diff(easyTrialOrderComparison,1,2);
easyCorrectResponse = nan(nEasyTrials,1);
easyCorrectResponse(easyTrialDiff <0) = 1;
easyCorrectResponse(easyTrialDiff >0) = 2;
easyCorrectResponse(easyTrialDiff==0) = randi(2);

%% Calculate the correct response for each warmup trial
if isFirstSession
    warmupTrialOrderComparison = imageComparison(warmupTrialOrder);
    warmupTrialDiff = diff(warmupTrialOrderComparison,1,2);
    warmupCorrectResponse = nan(nWarmupTrials,1);
    warmupCorrectResponse(warmupTrialDiff <0) = 1;
    warmupCorrectResponse(warmupTrialDiff >0) = 2;
    warmupCorrectResponse(warmupTrialDiff==0) = randi(2);
end

%% Set up vectors to keep track of trial info
%
% Set up vector for subject response per trial.
selectedResponse = nan(nTrials,1);
if isFirstSession
    warmupSelectedResponse = nan(nWarmupTrials,1);
end

% Set up cell array for second stimulus start time per trial
% (to calculate observer reaction time per trial).
reactionTimeStart = cell(nTrials,1);

% Set up cell array for observer response time per trial
% (to calculate observer reaction time per trial).
reactionTimeEnd = cell(nTrials,1);

%% Create mask pool from NoiseLevel0 images only (to be used for all noise levels)
%
% For each NoiseLevel0 image, take the average intensity (for each RGB 
% channel) per image block. Thus, the masks will have the same basic 
% luminance and color as the images.

%Get the number of blocks/image for the mask.
nBlocks = params.nBlocks;

% Load one image to get image size.
file1 = fullfile(pathToFolder,imageNames{1});
temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
imageSize = size(image1,1);

% The number of blocks must evenly divide the number of image pixels.
nPixels = imageSize;
if rem(nPixels,nBlocks)~=0
    error('params.nBlocks must evenly divide the number of image pixels.');
end

% Calculate the number of pixels per block.
blockPixels = nPixels/nBlocks;

% Set up matrices with mask pool info.
nMaskImages    = nConditions*nComparisons;
maskCondition  = nan(nMaskImages,1);
maskComparison = nan(nMaskImages,1);
maskPool       = nan(nBlocks,nBlocks,3,nMaskImages);

% Analyze each NoiseLevel0 image.
page = 1;
for ii = 1:nConditions
    for jj = 1:nComparisons
        imageThis = imageNoiseLevel==0 & ...
            imageCondition==conditions(ii) & ...
            imageComparison==comparisons(jj);
        
        % Load the image.
        file1 = fullfile(pathToFolder,imageNames{imageThis});
        temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
        
        % Flip image.
        image1 = image1(end:-1:1,:,:);
        
        % Take the average intensity (for each RGB channel) per block.
        for bii = 1:nBlocks
            for bjj = 1:nBlocks
                for kk = 1:3
                    theBlock = image1((bii-1)*blockPixels+1:bii*blockPixels,(bjj-1)*blockPixels+1:bjj*blockPixels,kk);
                    blockRGB = mean(theBlock(:));
                    maskPool(bii,bjj,kk,page) = blockRGB;
                end
            end
        end
        
        % Save other mask pool info.
        maskCondition(page)  = conditions(ii);
        maskComparison(page) = comparisons(jj);
        page = page+1;
    end
end

%% Begin task

% Note start time of experiment now.
startTime = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');

% Initialize task display.
[win, params] = initDisplay(params);

%% Start key capture and clear keyboard queue
ListenChar(2);
mglGetKeyEvent;

% Clear out any previous key presses.
FlushEvents;

%% Instantiate a gamePad object
if strcmp(controlSignal, 'gamePad')
    gamePad = GamePad();
else
    gamePad = [];
end

%% Enable start text
win.enableObject('instructions');
win.enableObject('keyOptions');
win.enableObject('startText');
win.draw;

%% Wait for key press ('2') to begin task (this is the 'GP:A' button on the gamepad)
key =[];
while isempty(key)
    % Get user response from keyboard.
    if strcmp(controlSignal, 'keyboard')
        key = mglGetKeyEvent;
        if ~isempty(key)
            switch key.charCode
                case {'2'}
                otherwise
                    key = [];
            end
        end
    % Get user response from gamePad.
    else
        key = gamePad.getKeyEvent();
        if ~isempty(key)
            switch key.charCode
                 % The instructions will request that the participant press
                 % the '2' button, because this button has the '2' symbol
                 % on the Philly gamepad.
                case {'GP:A'}
                otherwise
                    key = [];
            end
        end
    end
end

%% Turn off start text
win.disableObject('instructions');
win.disableObject('keyOptions');
win.disableObject('startText');
win.draw;

%% Run warmup trials: per trial, present images and wait for key press response
%
% If this is the first session, run warmup trials so participant can practice. The data will not be saved.
saveData    = 0;
keepLooping = 1;
iiTrial     = 0;
warmupquit  = false;

if isFirstSession
    
    % Display warmup text and wait for a key press to begin warmup trials.
    win.enableObject('warmupText');
    win.draw;
    FlushEvents;
    % Wait for key press.
    if strcmp(controlSignal, 'keyboard')
        key = [];
        while isempty(key)
            key = mglGetKeyEvent;
        end
    else
        key = [];
        while isempty(key)
            key = gamePad.getKeyEvent();
        end
    end
    win.disableObject('warmupText');
    win.draw;
    mglWaitSecs(params.ITI);
    % Reset the keyboard queue.
    mglGetKeyEvent;
    
    while keepLooping
        iiTrial = iiTrial + 1; %trial iteration
        
        % Get image index for 1st interval and for 2nd interval.
        idx1 = warmupTrialOrder(iiTrial,1);
        idx2 = warmupTrialOrder(iiTrial,2);
        
        % Get RGB image for 1st interval and for 2nd interval.
        file1 = fullfile(pathToFolder,imageNames{idx1});
        file2 = fullfile(pathToFolder,imageNames{idx2});
        temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
        temp = load(file2,'RGBImage'); image2 = temp.RGBImage; clear temp;
        
        % Flip images.
        image1 = image1(end:-1:1,:,:);
        image2 = image2(end:-1:1,:,:);
        
        % Create masks.
        image1condition  = imageCondition(idx1);
        image2condition  = imageCondition(idx2);
        image1comparison = imageComparison(idx1);
        image2comparison = imageComparison(idx2);
        maskIdx1 = find(maskCondition==image1condition & maskComparison==image1comparison);
        maskIdx2 = find(maskCondition==image2condition & maskComparison==image2comparison);
        mask1 = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels);
        mask2 = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels);
        
        % Write the images into the window and disable.
        win.addImage(params.image1Loc, params.image1Size, image1, 'Name', 'image1');
        win.addImage(params.image1Loc, params.image1Size, mask1,  'Name', 'mask1');
        win.addImage(params.image2Loc, params.image2Size, mask2,  'Name', 'mask2');
        win.addImage(params.image2Loc, params.image2Size, image2, 'Name', 'image2');
        win.disableObject('image1');
        win.disableObject('mask1');
        win.disableObject('mask2');
        win.disableObject('image2');
        
        % Enable 1st image and draw.
        win.enableObject('image1');
        win.draw;
        
        % Wait for stimulus duration.
        mglWaitSecs(params.stimDuration);
        win.disableObject('image1');
        
        % Enable 1st mask and draw.
        win.enableObject('mask1');
        win.draw;
        
        % Wait for ISI.
        mglWaitSecs(params.ISI);
        win.disableObject('mask1');
        
        % Enable 2nd mask and draw.
        win.enableObject('mask2');
        win.draw;
        
        % Wait for ISI.
        mglWaitSecs(params.ISI);
        win.disableObject('mask2');
        
        % Enable 2nd image and draw.
        win.enableObject('image2');
        win.draw;
        
        % Wait for stimulus duration.
        mglWaitSecs(params.stimDuration);
        win.disableObject('image2');
        win.draw;
        
        % Wait for key press response.
        FlushEvents;
        key =[];
        while isempty(key)
            % Get user response from keyboard.
            if strcmp(controlSignal, 'keyboard')
                key = mglGetKeyEvent;
                if ~isempty(key)
                    switch key.charCode
                        case {option1Key,option2Key}
                            warmupSelectedResponse(iiTrial) = getUserResponse(params,key);
                        case {'q'}
                            fprintf(2,'Do you want to quit? Type Y for Yes, otherwise give your response \n');
                            key2 = [];
                            while isempty(key2)
                                key2 = mglGetKeyEvent;
                                if ~isempty(key2)
                                    switch key2.charCode
                                        case {option1Key,option2Key}
                                            warmupSelectedResponse(iiTrial) = getUserResponse(params,key2);
                                        case {'y'}
                                            keepLooping = false;
                                        otherwise
                                            key2 = [];
                                    end
                                end
                            end
                        otherwise
                            key = [];
                    end
                end
                % Get user response from gamePad.
            else
                key = gamePad.getKeyEvent();
                if ~isempty(key)
                    switch key.charCode
                        case {option1Key,option2Key}
                            warmupSelectedResponse(iiTrial) = getUserResponse(params,key);
                        otherwise
                            key = [];
                    end
                end
                pressedKeyboard = mglGetKeyEvent;
                if ~isempty(pressedKeyboard)
                    switch pressedKeyboard.charCode
                        case {'q'}
                            fprintf(2,'Do you want to quit? Type Y for Yes, otherwise give your response using gamepad \n');
                            key2 = [];
                            keyG = [];
                            FlushEvents;
                            while isempty(key2) && isempty(keyG)
                                key2 = mglGetKeyEvent;
                                keyG = gamePad.getKeyEvent();
                                if ~isempty(key2)
                                    switch key2.charCode
                                        case {'y'}
                                            keepLooping = false;
                                            key = 0; % set key = 0 in order to exit this loop
                                        otherwise
                                            key2 = [];
                                    end
                                elseif ~isempty(keyG)
                                    switch keyG.charCode
                                        case {option1Key,option2Key}
                                            key = keyG;
                                            warmupSelectedResponse(iiTrial) = getUserResponse(params,key);
                                        otherwise
                                            keyG = [];
                                    end
                                end
                            end
                    end
                end
            end
        end
        
        % Check if the experiment continues, otherwise quit without saving data.
        if keepLooping
            fprintf('Selected interval in warmup trial: %d\n',warmupSelectedResponse(iiTrial));
            % Give feedback if option is on.
            if giveFeedback
                if warmupSelectedResponse(iiTrial) == warmupCorrectResponse(iiTrial)
                    sound(rightSound);
                else
                    sound(wrongSound);
                end
            end
            mglWaitSecs(params.ITI);
        else
            warmupquit = true;
            fprintf(2,'Quitting without saving any data.\n');
        end
        
        % Check if end of experiment is reached.
        if iiTrial == nWarmupTrials
            keepLooping = false;
        end
    end
    
    if ~warmupquit
        % Display warmup end text and wait for a key press to proceed.
        win.enableObject('warmupEndText');
        win.draw;
        FlushEvents;
        % Wait for key press.
        if strcmp(controlSignal, 'keyboard')
            key = [];
            while isempty(key)
                key = mglGetKeyEvent;
            end
        else
            key = [];
            while isempty(key)
                key = gamePad.getKeyEvent();
            end
        end
        win.disableObject('warmupEndText');
        win.draw;
        mglWaitSecs(params.ITI);
        % Reset the keyboard queue.
        mglGetKeyEvent;
    end
end

%% Run experiment: per trial, present images and wait for key press response
if ~warmupquit
    
    % Run practice trials to start (data not saved).
    nTrialsPractice     = 4; %number of practice trials to run
    keepLoopingPractice = 1;
    iiTrialPractice     = 0;
    while keepLoopingPractice
        iiTrialPractice = iiTrialPractice + 1; %practice trial iteration
        PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                      imageComparison,maskCondition,maskComparison,maskPool, ...
                      blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                      gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound);
        
        % Check if end of practice trials is reached.
        if iiTrialPractice == nTrialsPractice
            keepLoopingPractice = false;
        end
    end
    
    % Reset the keyboard queue.
    mglGetKeyEvent;
    
    saveData    = 1;
    keepLooping = 1;
    iiTrial     = 0;
    
    while keepLooping
        iiTrial = iiTrial + 1; %trial iteration
        
        % Get image index for 1st interval and for 2nd interval.
        idx1 = trialOrder(iiTrial,1);
        idx2 = trialOrder(iiTrial,2);
        
        % Get RGB image for 1st interval and for 2nd interval.
        file1 = fullfile(pathToFolder,imageNames{idx1});
        file2 = fullfile(pathToFolder,imageNames{idx2});
        temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
        temp = load(file2,'RGBImage'); image2 = temp.RGBImage; clear temp;
        
        % Flip images.
        image1 = image1(end:-1:1,:,:);
        image2 = image2(end:-1:1,:,:);
        
        % Create masks.
        image1condition  = imageCondition(idx1);
        image2condition  = imageCondition(idx2);
        image1comparison = imageComparison(idx1);
        image2comparison = imageComparison(idx2);
        maskIdx1 = find(maskCondition==image1condition & maskComparison==image1comparison);
        maskIdx2 = find(maskCondition==image2condition & maskComparison==image2comparison);
        mask1 = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels);
        mask2 = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels);

        % Write the images into the window and disable.
        win.addImage(params.image1Loc, params.image1Size, image1, 'Name', 'image1');
        win.addImage(params.image1Loc, params.image1Size, mask1,  'Name', 'mask1');
        win.addImage(params.image2Loc, params.image2Size, mask2,  'Name', 'mask2');
        win.addImage(params.image2Loc, params.image2Size, image2, 'Name', 'image2');
        win.disableObject('image1');
        win.disableObject('mask1');
        win.disableObject('mask2');
        win.disableObject('image2');
        
        % Enable 1st image and draw.
        win.enableObject('image1');
        win.draw;
        
        % Wait for stimulus duration.
        mglWaitSecs(params.stimDuration);
        win.disableObject('image1');
    
        % Enable 1st mask and draw.
        win.enableObject('mask1');
        win.draw;
        
        % Wait for ISI.
        mglWaitSecs(params.ISI);
        win.disableObject('mask1');
        
        % Enable 2nd mask and draw.
        win.enableObject('mask2');
        win.draw;
        
        % Wait for ISI.
        mglWaitSecs(params.ISI);
        win.disableObject('mask2');
        
        % Enable 2nd image and draw.
        win.enableObject('image2');
        win.draw;
        
        % Store 2nd image start time.
        reactionTimeStart{iiTrial} = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');
        
        % Wait for stimulus duration.
        mglWaitSecs(params.stimDuration);
        win.disableObject('image2');
        win.draw;

        % Wait for key press response.
        FlushEvents;
        key =[];
        while isempty(key)
            % Get user response from keyboard.
            if strcmp(controlSignal, 'keyboard')
                key = mglGetKeyEvent;
                if ~isempty(key)
                    switch key.charCode
                        case {option1Key,option2Key}
                            selectedResponse(iiTrial) = getUserResponse(params,key);
                            % Store observer response time.
                            reactionTimeEnd{iiTrial} = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');
                        case {'q'}
                            fprintf(2,'Do you want to quit? Type Y for Yes, otherwise give your response \n');
                            key2 = [];
                            while isempty(key2)
                                key2 = mglGetKeyEvent;
                                if ~isempty(key2)
                                    switch key2.charCode
                                        case {option1Key,option2Key}
                                            selectedResponse(iiTrial) = getUserResponse(params,key2);
                                            % Store observer response time.
                                            reactionTimeEnd{iiTrial} = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');
                                        case {'y'}
                                            keepLooping = false;
                                        otherwise
                                            key2 = [];
                                    end
                                end
                            end
                        otherwise
                            key = [];
                    end
                end
                % Get user response from gamePad.
            else
                key = gamePad.getKeyEvent();
                if ~isempty(key)
                    switch key.charCode
                        case {option1Key,option2Key}
                            selectedResponse(iiTrial) = getUserResponse(params,key);
                            % Store observer response time.
                            reactionTimeEnd{iiTrial} = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');
                        otherwise
                            key = [];
                    end
                end
                pressedKeyboard = mglGetKeyEvent;
                if ~isempty(pressedKeyboard)
                    switch pressedKeyboard.charCode
                        case {'q'}
                            fprintf(2,'Do you want to quit? Type Y for Yes, otherwise give your response using gamepad \n');
                            key2 = [];
                            keyG = [];
                            FlushEvents;
                            while isempty(key2) && isempty(keyG)
                                key2 = mglGetKeyEvent;
                                keyG = gamePad.getKeyEvent();
                                if ~isempty(key2)
                                    switch key2.charCode
                                        case {'y'}
                                            keepLooping = false;
                                            key = 0; % set key = 0 in order to exit this loop
                                        otherwise
                                            key2 = [];
                                    end
                                elseif ~isempty(keyG)
                                    switch keyG.charCode
                                        case {option1Key,option2Key}
                                            key = keyG;
                                            selectedResponse(iiTrial) = getUserResponse(params,key);
                                            % Store observer response time.
                                            reactionTimeEnd{iiTrial} = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');
                                        otherwise
                                            keyG = [];
                                    end
                                end
                            end
                    end
                end
            end
        end

        % Check if the experiment continues, otherwise quit without saving data.
        if keepLooping
            fprintf('Selected interval: %d\n',selectedResponse(iiTrial));
            % Give feedback if option is on.
            if giveFeedback
                if selectedResponse(iiTrial) == correctResponse(iiTrial)
                    sound(rightSound);
                else
                    sound(wrongSound);
                end
            end
            mglWaitSecs(params.ITI);
        else
            fprintf(2,'Quitting without saving any data.\n');
            saveData = 0;
        end
        
        % Check if 1/6 of experiment is reached.
        if iiTrial == ceil(nTrials/6)
            win.enableObject('oneSixthText');
            win.draw;
            pause(60);
            win.disableObject('oneSixthText');
            win.enableObject('restOver');
            win.draw;
            FlushEvents;
            % Wait for key press.
            if strcmp(controlSignal, 'keyboard')
                key = [];
                while isempty(key)
                    key = mglGetKeyEvent;
                end
            else
                key = [];
                while isempty(key)
                    key = gamePad.getKeyEvent();
                end
            end
            % Turn off text.
            win.disableObject('restOver');
            win.draw;
            mglWaitSecs(params.ITI);
            % Run one practice trial (data not saved).
            PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                          imageComparison,maskCondition,maskComparison,maskPool, ...
                          blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                          gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound);
            % Reset the keyboard queue.
            mglGetKeyEvent;
        end
        
        % Check if 2/6 of experiment is reached.
        if iiTrial == ceil(2*nTrials/6)
            win.enableObject('twoSixthsText');
            win.draw;
            pause(60);
            win.disableObject('twoSixthsText');
            win.enableObject('restOver');
            win.draw;
            FlushEvents;
            % Wait for key press.
            if strcmp(controlSignal, 'keyboard')
                key = [];
                while isempty(key)
                    key = mglGetKeyEvent;
                end
            else
                key = [];
                while isempty(key)
                    key = gamePad.getKeyEvent();
                end
            end
            % Turn off text.
            win.disableObject('restOver');
            win.draw;
            mglWaitSecs(params.ITI);
            % Run one practice trial (data not saved).
            PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                          imageComparison,maskCondition,maskComparison,maskPool, ...
                          blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                          gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound);
            % Reset the keyboard queue.
            mglGetKeyEvent;
        end
        
        % Check if 3/6 of experiment is reached.
        if iiTrial == ceil(3*nTrials/6)
            win.enableObject('threeSixthsText');
            win.draw;
            pause(60);
            win.disableObject('threeSixthsText');
            win.enableObject('restOver');
            win.draw;
            FlushEvents;
            % Wait for key press.
            if strcmp(controlSignal, 'keyboard')
                key = [];
                while isempty(key)
                    key = mglGetKeyEvent;
                end
            else
                key = [];
                while isempty(key)
                    key = gamePad.getKeyEvent();
                end
            end
            % Turn off text.
            win.disableObject('restOver');
            win.draw;
            mglWaitSecs(params.ITI);
            % Run one practice trial (data not saved).
            PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                          imageComparison,maskCondition,maskComparison,maskPool, ...
                          blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                          gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound);
            % Reset the keyboard queue.
            mglGetKeyEvent;
        end
        
        % Check if 4/6 of experiment is reached.
        if iiTrial == ceil(4*nTrials/6)
            win.enableObject('fourSixthsText');
            win.draw;
            pause(60);
            win.disableObject('fourSixthsText');
            win.enableObject('restOver');
            win.draw;
            FlushEvents;
            % Wait for key press.
            if strcmp(controlSignal, 'keyboard')
                key = [];
                while isempty(key)
                    key = mglGetKeyEvent;
                end
            else
                key = [];
                while isempty(key)
                    key = gamePad.getKeyEvent();
                end
            end
            % Turn off text.
            win.disableObject('restOver');
            win.draw;
            mglWaitSecs(params.ITI);
            % Run one practice trial (data not saved).
            PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                          imageComparison,maskCondition,maskComparison,maskPool, ...
                          blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                          gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound);
            % Reset the keyboard queue.
            mglGetKeyEvent;
        end
        
        % Check if 5/6 of experiment is reached.
        if iiTrial == ceil(5*nTrials/6)
            win.enableObject('fiveSixthsText');
            win.draw;
            pause(60);
            win.disableObject('fiveSixthsText');
            win.enableObject('restOver');
            win.draw;
            FlushEvents;
            % Wait for key press.
            if strcmp(controlSignal, 'keyboard')
                key = [];
                while isempty(key)
                    key = mglGetKeyEvent;
                end
            else
                key = [];
                while isempty(key)
                    key = gamePad.getKeyEvent();
                end
            end
            % Turn off text.
            win.disableObject('restOver');
            win.draw;
            mglWaitSecs(params.ITI);
            % Run one practice trial (data not saved).
            PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                          imageComparison,maskCondition,maskComparison,maskPool, ...
                          blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                          gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound);
            % Reset the keyboard queue.
            mglGetKeyEvent;
        end
        
        % Check if end of experiment is reached.
        if iiTrial == nTrials
            keepLooping = false;
        end
    end
end
%% Close up experiment
%
% Show end of experiment text.
win.enableObject('finishText');
win.draw;
pause(10);
win.disableObject('finishText');

% Close task display.
win.close;

% Make sure key capture is off.
ListenChar(0);

% Note end time of experiment now.
endTime = datestr(now,'mm/dd/yyyy HH:MM:SS.FFF');

%% Save experimental parameters and psychophysical responses
%
% Data won't be saved if this option is on.
if isDemo
   saveData = 0;
end
    
if saveData
    % Save data in 'data' struct.
    data = struct;

    % Save inputted experimental parameters.
    data.experimentName = experimentName;
    data.subjectName    = subjectName;
    data.nIterations    = nIterations;
    data.controlSignal  = controlSignal;
    data.option1Key     = option1Key;
    data.option2Key     = option2Key;
    data.giveFeedback   = giveFeedback;
    
    % Save task parameters.
    data.screenDimsCm = params.screenDimsCm;
    data.bgColor      = params.bgColor;
    data.textColor    = params.textColor;
    data.image1Loc    = params.image1Loc;
    data.image2Loc    = params.image2Loc;
    data.image1Size   = params.image1Size;
    data.image2Size   = params.image2Size;
    data.ISI          = params.ISI;
    data.ITI          = params.ITI;
    data.stimDuration = params.stimDuration;
    
    % Save task start time and end time.
    data.startTime = startTime;
    data.endTime   = endTime;
    
    % Save image information.
    data.imageNames        = imageNames;
    data.imageNoiseLevel   = imageNoiseLevel;
    data.imageNoiseAmount1 = imageNoiseAmount1;
    data.imageNoiseAmount2 = imageNoiseAmount2;
    data.imageCondition    = imageCondition;
    data.imageComparison   = imageComparison;
    
    % Save trial information.
    data.trialOrder           = trialOrder;
    data.trialOrderComparison = trialOrderComparison;
    data.trialNoiseLevel      = trialNoiseLevel;
    data.trialNoiseAmount1    = trialNoiseAmount1;
    data.trialNoiseAmount2    = trialNoiseAmount2;
    data.correctResponse      = correctResponse;
    data.selectedResponse     = selectedResponse;
    data.reactionTimeStart    = reactionTimeStart;
    data.reactionTimeEnd      = reactionTimeEnd;
    
    % Create output folder if it doesn't exist.
    if ~exist(pathToOutput, 'dir')
        mkdir(pathToOutput);
    end

    % Set path to specific output folder where data will be saved.
    pathToOutputFolder = fullfile(pathToOutput,sprintf('%s%s','subject',subjectName));
    if ~exist(pathToOutputFolder, 'dir')
        mkdir(pathToOutputFolder);
    end
    
    % Set path to the file to save.
    fileNum = GetNextDataFileNumber(pathToOutputFolder,'.mat');
    fileName = sprintf('%s%s_%d.mat','data',subjectName,fileNum);
    pathToOutputFile = fullfile(pathToOutputFolder,fileName);
    
    % Save data file.
    save(pathToOutputFile,'data');
    fprintf('\nData was saved in:\n%s\n', pathToOutputFile);

    % Set acquisition status.
    acquisitionStatus = 1;
end

end

%% Helper functions

%% Initialize display
function [win, params] = initDisplay(params)

% Create the GLWindow object.
win = GLWindow('SceneDimensions',params.screenDimsCm,'BackgroundColor',params.bgColor);

try
    % Open the display.
    win.open;

    % Add instructions text.
    win.addText('Compared to the 1st banana, is the 2nd banana to the left or right?', ...
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'instructions'); % Identifier for the object

    % Add key option text.
    key1 = params.option1Key;
    key2 = params.option2Key;
    if strcmp(key1,'GP:UpperLeftTrigger');  key1 = 'Upper Left Trigger';  end
    if strcmp(key2,'GP:UpperRightTrigger'); key2 = 'Upper Right Trigger'; end
    if strcmp(key1,'GP:X'); key1 = 'gamepad X'; end
    if strcmp(key2,'GP:A'); key2 = 'gamepad A'; end
    win.addText(['If to left -> ', key1, '     If to right -> ', key2], ... % Text to display
        'Center', [0 5], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ...  % RGB color
        'Name', 'keyOptions'); % Identifier for the object

    % Add start text.
    win.addText('Hit the ''2'' button to start.', ... % Text to display
        'Center', [0 -8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'startText'); % Identifier for the object

    % Add warmup start text.
    win.addText('Today, you will start with some practice trials. Hit any button to continue.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'warmupText'); % Identifier for the object
    
    % Add warmup end text.
    win.addText('Your practice trials are complete. Hit any button to begin the experiment.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'warmupEndText'); % Identifier for the object
    
    % Add text for when experiment is 1/6 over.
    win.addText('One sixth of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'oneSixthText'); % Identifier for the object
    
    % Add text for when rest period is over.
    win.addText('Rest time complete. Hit any button to continue.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'restOver'); % Identifier for the object
    
    % Add text for when experiment is 2/6 over.
    win.addText('One third of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'twoSixthsText'); % Identifier for the object
    
    % Add text for when experiment is 3/6 over.
    win.addText('One half trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'threeSixthsText'); % Identifier for the object
    
    % Add text for when experiment is 4/6 over.
    win.addText('Two thirds of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'fourSixthsText'); % Identifier for the object
    
    % Add text for when experiment is 5/6 over.
    win.addText('Five sixths of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'fiveSixthsText'); % Identifier for the object

    % Add text for when experiment is complete.
    win.addText('Experiment is complete.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'finishText'); % Identifier for the object

    % Turn all objects off for now.
    win.disableAllObjects;
    
catch e
    win.close;
    rethrow(e);
end
end

%% Get user response
function response = getUserResponse(params,key)
response = [];
switch key.charCode
    % Option 1 key response.
    case params.option1Key
        response = 1;
        
    % Option 2 key response.
    case params.option2Key
        response = 2;
end
end

%% Create mask for interstimulus interval
%
% Get 2 quantized images (average intensity calculated per block) from the
% maskPool based on the first 2 inputs (mask pool indices). Per mask block,
% randomly draw an average intensity block from one quantized image or the
% other.

function mask = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels)

% Get number of blocks used to create the quantized images in the maskPool.
nBlocks = size(maskPool,1);

% Calculate the mask size (number of blocks * number of pixels per block).
maskSize = nBlocks * blockPixels;

% Create mask.
mask = zeros(maskSize,maskSize,3);
for ii = 1:nBlocks
    for jj = 1:nBlocks
        for kk = 1:3
            if CoinFlip(1,0.5)
                blockRGB = maskPool(ii,jj,kk,maskIdx1);
            else
                blockRGB = maskPool(ii,jj,kk,maskIdx2);
            end
            mask((ii-1)*blockPixels+1:ii*blockPixels,(jj-1)*blockPixels+1:jj*blockPixels,kk) = blockRGB;
        end
    end
end
end

%% Run a post-break practice trial
%
% Following each break, one practice trial is run. The practice trial is
% selected randomly from the easy trials run at the start of the
% experiment. The data from this practice trial are not saved.

function PracticeTrial(easyTrialOrder,pathToFolder,imageNames,imageCondition, ...
                       imageComparison,maskCondition,maskComparison,maskPool, ...
                       blockPixels,win,controlSignal,option1Key,option2Key,params, ...
                       gamePad,giveFeedback,easyCorrectResponse,rightSound,wrongSound)
                   
% Get image index for 1st interval and for 2nd interval.
iiEasyTrial = randi(size(easyTrialOrder,1));
idx1 = easyTrialOrder(iiEasyTrial,1);
idx2 = easyTrialOrder(iiEasyTrial,2);

% Get RGB image for 1st interval and for 2nd interval.
file1 = fullfile(pathToFolder,imageNames{idx1});
file2 = fullfile(pathToFolder,imageNames{idx2});
temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
temp = load(file2,'RGBImage'); image2 = temp.RGBImage; clear temp;

% Flip images.
image1 = image1(end:-1:1,:,:);
image2 = image2(end:-1:1,:,:);

% Create masks.
image1condition  = imageCondition(idx1);
image2condition  = imageCondition(idx2);
image1comparison = imageComparison(idx1);
image2comparison = imageComparison(idx2);
maskIdx1 = find(maskCondition==image1condition & maskComparison==image1comparison);
maskIdx2 = find(maskCondition==image2condition & maskComparison==image2comparison);
mask1 = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels);
mask2 = MakeBlockMask(maskIdx1,maskIdx2,maskPool,blockPixels);

% Write the images into the window and disable.
win.addImage(params.image1Loc, params.image1Size, image1, 'Name', 'image1');
win.addImage(params.image1Loc, params.image1Size, mask1,  'Name', 'mask1');
win.addImage(params.image2Loc, params.image2Size, mask2,  'Name', 'mask2');
win.addImage(params.image2Loc, params.image2Size, image2, 'Name', 'image2');
win.disableObject('image1');
win.disableObject('mask1');
win.disableObject('mask2');
win.disableObject('image2');

% Enable 1st image and draw.
win.enableObject('image1');
win.draw;

% Wait for stimulus duration.
mglWaitSecs(params.stimDuration);
win.disableObject('image1');

% Enable 1st mask and draw.
win.enableObject('mask1');
win.draw;

% Wait for ISI.
mglWaitSecs(params.ISI);
win.disableObject('mask1');

% Enable 2nd mask and draw.
win.enableObject('mask2');
win.draw;

% Wait for ISI.
mglWaitSecs(params.ISI);
win.disableObject('mask2');

% Enable 2nd image and draw.
win.enableObject('image2');
win.draw;

% Wait for stimulus duration.
mglWaitSecs(params.stimDuration);
win.disableObject('image2');
win.draw;

% Wait for key press response.
FlushEvents;
key =[];
while isempty(key)
    % Get user response from keyboard.
    if strcmp(controlSignal, 'keyboard')
        key = mglGetKeyEvent;
        if ~isempty(key)
            switch key.charCode
                case {option1Key,option2Key}
                    practiceSelectedResponse = getUserResponse(params,key);
                otherwise
                    key = [];
            end
        end
        % Get user response from gamePad.
    else
        key = gamePad.getKeyEvent();
        if ~isempty(key)
            switch key.charCode
                case {option1Key,option2Key}
                    practiceSelectedResponse = getUserResponse(params,key);
                otherwise
                    key = [];
            end
        end
    end
end

fprintf('Selected interval in practice trial: %d\n',practiceSelectedResponse);
% Give feedback if option is on.
if giveFeedback
    if practiceSelectedResponse == easyCorrectResponse(iiEasyTrial)
        sound(rightSound);
    else
        sound(wrongSound);
    end
end
end

%% End
