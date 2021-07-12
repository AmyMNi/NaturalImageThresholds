function acquisitionStatus = runNaturalImageExperiment(varargin)
%runNaturalImageExperiment
%
% Usage:
%   runNaturalImageExperiment('experimentName', 'Experiment000', ...
%                             'subjectName', 'AN');
%
% Description:
%   Run the natural image psychophysics experiment given the specified
%   experiment folder. Save the experimental parameters and psychophysical 
%   responses (struct 'data') in the specified output folder.
%
% Notes:
%   1) Run experiment on calibrated 27-in NEC MultiSync PA271 monitor
%   2) Center observer's eyes horizontal and vertically using chin and forehead rest
%   3) Set monitor distance to 75 cm
%   4) Have observer dark adapt in the dark room for 5 min
%   5) Instruct observers to:
%      a) Keep eyes centered on the screen
%      b) Press firmly on the gamepad buttons
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'subjectName'    : (string)  Name of subject (default: 'AN')
%   'nIterations'    : (scalar)  Number of iterations per image comparison *must be an even number (default: 14)
%   'controlSignal'  : (string)  Input method for user response (options: 'gamePad', 'keyboard') (default: 'gamePad')
%   'option1Key'     : (string)  For gamePad either 'GP:UpperLeftTrigger'  or 'GP:X', for keyboard -> '1' (default: 'GP:UpperLeftTrigger')
%   'option2Key'     : (string)  For gamePad either 'GP:UpperRightTrigger' or 'GP:A', for keyboard -> '2' (default: 'GP:UpperRightTrigger')
%   'giveFeedback'   : (logical) Give feedback if option is on (default: true)
%   'isDemo'         : (logical) Data won't be saved if on (default: false)
%
% History:
%   06/07/21  amn  Adapted from BrainardLab/VirtualWorldPsychophysics
%   07/05/21  amn  Experiment edits based on pilot tests

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.addParameter('subjectName', 'AN', @ischar);
parser.addParameter('nIterations', 14, @isscalar);
parser.addParameter('controlSignal', 'gamePad', @ischar);
parser.addParameter('option1Key', 'GP:UpperLeftTrigger', @ischar);
parser.addParameter('option2Key', 'GP:UpperRightTrigger', @ischar);
parser.addParameter('giveFeedback', true, @islogical);
parser.addParameter('isDemo', false, @islogical);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
subjectName    = parser.Results.subjectName;
nIterations    = parser.Results.nIterations;
controlSignal  = parser.Results.controlSignal;
option1Key     = parser.Results.option1Key;
option2Key     = parser.Results.option2Key;
giveFeedback   = parser.Results.giveFeedback;
isDemo         = parser.Results.isDemo;

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

%% Set paths to folders
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Set path to image folder.
pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBs');

% Set path to output folder.
pathToOutput = fullfile(getpref(projectName,'BaseDir'),experimentName,'PsychophysicalData');
if ~exist(pathToOutput, 'dir')
    mkdir(pathToOutput);
end

%% Set up experiment parameters
%
% Specify feedback sounds.
rightSound = sin(2*pi*(1:1000)/10)/10;
wrongSound = rand(1,1000).*ceil(sin(2*pi*(1:1000)/10))/10;

% Set acquisition status.
acquisitionStatus = 0;

% Set task parameters.
params.screenDimsCm = [59.67 33.57]; %cm
params.fpSize       = [0.1 0.1]; % fixation point size
params.fpColor      = [174 174 128]/255; % fixation point: color of banana (unused)
params.fpColorRed   = [0.6 0.2 0.2]; % fixation point color red
params.bgColor      = [128 128 128]/255; % to match electrophys task
params.textColor    = [0.6 0.2 0.2];
params.image1Loc  = [0 0];
params.image2Loc  = [0 0];
params.image1Size = [10.54 10.54]; % monitor distance=75cm: scene 8 deg vis angle (target 4 deg)
params.image2Size = [10.54 10.54];
params.ISI          = 0.20; % seconds
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
imageNoiseLevel  = nan(nImages,1);
imageNoiseAmount = nan(nImages,1);
imageCondition   = nan(nImages,1);
imageComparison  = nan(nImages,1);

% Break down image file names to get image info.
for ii = 1:nImages
    name = imageNames{ii};
    p    = strfind(name,'_');

    % Label the 'imageNames' indices by their condition, comparison amount,
    % noise level, and noise amount.
    s1 = name(1:p(1));
    s2 = name(p(1):p(2));
    s3 = name(p(2):p(3));
    s4 = name(p(3):end);
    
    imageCondition(ii)   = str2double(regexp(s1,'[+-]?\d*','Match'));
    imageComparison(ii)  = str2double(regexp(s2,'[+-]?\d*','Match'));
    imageNoiseLevel(ii)  = str2double(regexp(s3,'[+-]?\d*','Match'));
    imageNoiseAmount(ii) = str2double(regexp(s4,'[+-]?\d*','Match'));
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

% Preallocate vectors of noise level and noise amount (within the level) per trial.
trialNoiseLevel  = nan(nTrials,1);
trialNoiseAmount = nan(nTrials,2);

% Create random order of runs (each noise level divided into two runs).
allRunsOrdered = [noiseLevels; noiseLevels];
allRuns        = allRunsOrdered(randperm(numel(allRunsOrdered)));

% Create a trial order per run, then combine across runs.
for jj = 1:numel(allRuns)
    noiseLevelthis      = allRuns(jj);
    trialOrderRun       = nan(nTrialsRun,2);
    trialNoiseAmountRun = nan(nTrialsRun,2);
    rrow = 1;
    
    % Include nIterations/2 blocks per run.
    % Complete a block before moving on to the next block.
    % Within each block, create a randomized order of trials.
    for ii = 1:nIterations/2
        trialsPerBlock       = nan(nTrialsBlock,2);
        noiseAmountsPerBlock = nan(nTrialsBlock,2);
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
            trialsPerCondition      = nan(nComparisons,2);
            noiseAmountPerCondition = nan(nComparisons,2);
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
                thisIndices      = [centerthis comparisonthis];
                thisNoiseAmounts = imageNoiseAmount(thisIndices)';
                orderrand        = randperm(2);
                trialsPerCondition(iii,:)      = thisIndices(orderrand);
                noiseAmountPerCondition(iii,:) = thisNoiseAmounts(orderrand);
            end
            
            % Combine ordered trials across all conditions for a single block.
            trialsPerBlock      (brow:brow+nComparisons-1,:) = trialsPerCondition;
            noiseAmountsPerBlock(brow:brow+nComparisons-1,:) = noiseAmountPerCondition;
            brow = brow+nComparisons;
        end
        
        % Randomize trial order within this single block.
        trialsrand = randperm(nTrialsBlock);
        trialsPerBlockRand       = trialsPerBlock      (trialsrand,:);
        noiseAmountsPerBlockRand = noiseAmountsPerBlock(trialsrand,:);
        
        % Combine iterations of blocks for this run.
        trialOrderRun      (rrow:rrow+nTrialsBlock-1,:) = trialsPerBlockRand;
        trialNoiseAmountRun(rrow:rrow+nTrialsBlock-1,:) = noiseAmountsPerBlockRand;
        rrow = rrow+nTrialsBlock;
    end
    
    % Combine trials across runs.
    trialOrder      (srow:srow+nTrialsRun-1,:) = trialOrderRun;
    trialNoiseLevel (srow:srow+nTrialsRun-1,1) = noiseLevelthis;
    trialNoiseAmount(srow:srow+nTrialsRun-1,:) = trialNoiseAmountRun;
    srow = srow+nTrialsRun;
end

%% Create a trial order for the easy trials
nEasyTrials = 4;

% Preallocate matrix of trial order.
easyTrialOrder = nan(nEasyTrials,2);
trow = 1;

% Create an ordered list of image index comparisons for easy comparisons.
for jjj = 1:nConditions
    centerpos = conditions(jjj);
    centerIdx = find(imageComparison==0 & imageCondition==centerpos & ...
                imageNoiseLevel==0);
    compIdx1  = find(imageComparison==comparisons(1) & imageCondition==centerpos & ...
                imageNoiseLevel==0);
    thisIndices = [centerIdx compIdx1];
    easyTrialOrder(trow,:) = thisIndices(randperm(2));
    trow = trow+1;
    compIdx2  = find(imageComparison==comparisons(end) & imageCondition==centerpos & ...
                imageNoiseLevel==0);
    thisIndices = [centerIdx compIdx2];
    easyTrialOrder(trow,:) = thisIndices(randperm(2));
    trow = trow+1;
end

% Randomize easy trial order.
easyTrialOrder = easyTrialOrder(randperm(nEasyTrials),:);
        
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

%% Set up vectors to keep track of trial info
%
% Set up vector for subject response per trial.
selectedResponse     = nan(nTrials,1);
easySelectedResponse = nan(nEasyTrials,1);

% Set up cell array for second stimulus start time per trial
% (to calculate observer reaction time per trial).
reactionTimeStart = cell(nTrials,1);

% Set up cell array for observer response time per trial
% (to calculate observer reaction time per trial).
reactionTimeEnd = cell(nTrials,1);

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

%% Enable fixation and start text
win.enableObject('instructions');
win.enableObject('keyOptions');
win.enableObject('startText');
win.draw;

%% Wait for key press (any) to begin task
key = [];
while isempty(key)
    if strcmp(controlSignal, 'keyboard')
        key = mglGetKeyEvent;
    else
        key = gamePad.getKeyEvent();
    end
end

%% Turn off start text
win.disableObject('instructions');
win.disableObject('keyOptions');
win.disableObject('startText');
win.draw;

%% Run easy trials: per trial, present images and wait for key press response
%
% Run 5 easy trials to acclimate the subject. The data will not be saved.
saveData    = 0;
keepLooping = 1;
iiTrial     = 0;
easyquit    = false;

while keepLooping
    iiTrial = iiTrial + 1; %trial iteration
    
    % Get image index for 1st interval and for 2nd interval.
    idx1 = easyTrialOrder(iiTrial,1);
    idx2 = easyTrialOrder(iiTrial,2);
    
    % Get RGB image for 1st interval and for 2nd interval.
    file1 = fullfile(pathToFolder,imageNames{idx1});
    file2 = fullfile(pathToFolder,imageNames{idx2});
    temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
    temp = load(file2,'RGBImage'); image2 = temp.RGBImage; clear temp;
    
    % Flip images.
    image1 = image1(end:-1:1,:,:);
    image2 = image2(end:-1:1,:,:);

    % Create masks for the ISI.
    mask1 = MakeBlockMask(image1,image2,params.nBlocks);
    mask2 = MakeBlockMask(image1,image2,params.nBlocks);
    mask3 = MakeBlockMask(image1,image2,params.nBlocks);
    
    % Write the images into the window and disable.
    win.addImage(params.image1Loc, params.image1Size, image1, 'Name', 'image1');
    win.addImage(params.image1Loc, params.image1Size, mask1,  'Name', 'mask1');
    win.addImage(params.image1Loc, params.image1Size, mask2,  'Name', 'mask2');
    win.addImage(params.image1Loc, params.image1Size, mask3,  'Name', 'mask3');
    win.addImage(params.image2Loc, params.image2Size, image2, 'Name', 'image2');
    win.disableObject('image1');
    win.disableObject('mask1');
    win.disableObject('mask2');
    win.disableObject('mask3');
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
    
    % Enable 3rd mask and draw.
    win.enableObject('mask3');
    win.draw;
    
    % Wait for ISI.
    mglWaitSecs(params.ISI);
    win.disableObject('mask3');
    
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
                        easySelectedResponse(iiTrial) = getUserResponse(params,key);
                    case {'q'}
                        fprintf(2,'Do you want to quit? Type Y for Yes, otherwise give your response \n');
                        key2 = [];
                        while isempty(key2)
                            key2 = mglGetKeyEvent;
                            if ~isempty(key2)
                                switch key2.charCode
                                    case {option1Key,option2Key}
                                        easySelectedResponse(iiTrial) = getUserResponse(params,key2);
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
                        easySelectedResponse(iiTrial) = getUserResponse(params,key);
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
                                        easySelectedResponse(iiTrial) = getUserResponse(params,key);
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
        fprintf('Selected interval: %d\n',easySelectedResponse(iiTrial));
        % Give feedback if option is on.
        if giveFeedback
            if easySelectedResponse(iiTrial) == easyCorrectResponse(iiTrial)
                sound(rightSound);
            else
                sound(wrongSound);
            end
        end
        mglWaitSecs(params.ITI);        
    else
        easyquit = true;
        fprintf(2,'Quitting without saving any data.\n');
    end

    % Check if end of experiment is reached.
    if iiTrial == nEasyTrials
        keepLooping = false;
    end
end

%% Run experiment: per trial, present images and wait for key press response
if ~easyquit
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
        
        % Create masks for the ISI.
        mask1 = MakeBlockMask(image1,image2,params.nBlocks);
        mask2 = MakeBlockMask(image1,image2,params.nBlocks);
        mask3 = MakeBlockMask(image1,image2,params.nBlocks);
        
        % Write the images into the window and disable.
        win.addImage(params.image1Loc, params.image1Size, image1, 'Name', 'image1');
        win.addImage(params.image1Loc, params.image1Size, mask1,  'Name', 'mask1');
        win.addImage(params.image1Loc, params.image1Size, mask2,  'Name', 'mask2');
        win.addImage(params.image1Loc, params.image1Size, mask3,  'Name', 'mask3');
        win.addImage(params.image2Loc, params.image2Size, image2, 'Name', 'image2');
        win.disableObject('image1');
        win.disableObject('mask1');
        win.disableObject('mask2');
        win.disableObject('mask3');
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
        
        % Enable 3rd mask and draw.
        win.enableObject('mask3');
        win.draw;
        
        % Wait for ISI.
        mglWaitSecs(params.ISI);
        win.disableObject('mask3');

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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: includes both options of 4 runs or 6 runs
        if nNoiseLevels<3
            % Check if one quarter of experiment is reached.
            if iiTrial == ceil(nTrials/4)
                win.enableObject('oneQuarterText');
                win.draw;
                pause(60);
                win.disableObject('oneQuarterText');
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
                % Reset the keyboard queue.
                mglGetKeyEvent;
            end
            
            % Check if two quarters of experiment is reached.
            if iiTrial == ceil(2*nTrials/4)
                win.enableObject('twoQuartersText');
                win.draw;
                pause(60);
                win.disableObject('twoQuartersText');
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
                % Reset the keyboard queue.
                mglGetKeyEvent;
            end
            
            % Check if three quarters of experiment is reached.
            if iiTrial == ceil(3*nTrials/4)
                win.enableObject('threeQuartersText');
                win.draw;
                pause(60);
                win.disableObject('threeQuartersText');
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
                % Reset the keyboard queue.
                mglGetKeyEvent;
            end
        else
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
                % Reset the keyboard queue.
                mglGetKeyEvent;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
    data.fpSize       = params.fpSize;
    data.fpColor      = params.fpColor;
    data.fpColorRed   = params.fpColorRed;
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
    data.imageNames       = imageNames;
    data.imageNoiseLevel  = imageNoiseLevel;
    data.imageNoiseAmount = imageNoiseAmount;
    data.imageCondition   = imageCondition;
    data.imageComparison  = imageComparison;
    
    % Save trial information.
    data.trialOrder           = trialOrder;
    data.trialOrderComparison = trialOrderComparison;
    data.trialNoiseLevel      = trialNoiseLevel;
    data.trialNoiseAmount     = trialNoiseAmount;
    data.correctResponse      = correctResponse;
    data.selectedResponse     = selectedResponse;
    data.reactionTimeStart    = reactionTimeStart;
    data.reactionTimeEnd      = reactionTimeEnd;

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
    
    % Add the fixation point.
    win.addOval([0 0], params.fpSize, params.fpColor, 'Name', 'fp');
    
    % Add the fixation point in red color.
    win.addOval([0 0], params.fpSize, params.fpColorRed, 'Name', 'fpRed');

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
    win.addText('Hit any button to start.', ... % Text to display
        'Center', [0 -8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'startText'); % Identifier for the object
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: includes text for both options of 4 runs or 6 runs
    
    % Add text for when experiment is one quarter over.
    win.addText('One quarter of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'oneQuarterText'); % Identifier for the object

    % Add text for when rest period is over.
    win.addText('Rest time complete. Hit any button to continue.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'restOver'); % Identifier for the object

    % Add text for when experiment is two quarters over.
    win.addText('Half of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'twoQuartersText'); % Identifier for the object
    
    % Add text for when experiment is three quarters over.
    win.addText('Three quarters of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'threeQuartersText'); % Identifier for the object
    
    % Add text for when experiment is 1/6 over.
    win.addText('One sixth of trials complete. Take a minute to stand or stretch.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'oneSixthText'); % Identifier for the object
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
% Take the average intensity per block for each image, then randomly draw
% each block from one image or the other. Thus, the mask will have the same
% basic luminance and color as the images.
function mask = MakeBlockMask(image1,image2,nBlocks)

% The number of blocks must evenly divide the number of image pixels.
nPixels = size(image1,1);
if rem(nPixels,nBlocks)~=0
    error('params.nBlocks must evenly divide the number of image pixels.');
end

% Calculate the number of pixels per block.
blockPixels = nPixels/nBlocks;

% Per mask block, randomly draw an average intensity block from one image or the other.
mask = zeros(size(image1));
for ii = 1:nBlocks
    for jj = 1:nBlocks
        for kk = 1:3
            if CoinFlip(1,0.5)
                theBlock = image1((ii-1)*blockPixels+1:ii*blockPixels,(jj-1)*blockPixels+1:jj*blockPixels,kk);
                blockRGB = mean(theBlock(:));
                mask((ii-1)*blockPixels+1:ii*blockPixels,(jj-1)*blockPixels+1:jj*blockPixels,kk) = blockRGB;
            else
                theBlock = image2((ii-1)*blockPixels+1:ii*blockPixels,(jj-1)*blockPixels+1:jj*blockPixels,kk);
                blockRGB = mean(theBlock(:));
                mask((ii-1)*blockPixels+1:ii*blockPixels,(jj-1)*blockPixels+1:jj*blockPixels,kk) = blockRGB;
            end
        end
    end
end
end

%% End