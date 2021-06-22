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
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'nIterations'    : (scalar)  Number of iterations per image comparison (default: 15)
%   'controlSignal'  : (string)  Input method for user response (options: 'gamePad', 'keyboard') (default: 'gamePad')
%   'option1Key'     : (string)  For gamePad either 'GP:UpperLeftTrigger'  or 'GP:X', for keyboard -> '1' (default: 'GP:UpperLeftTrigger')
%   'option2Key'     : (string)  For gamePad either 'GP:UpperRightTrigger' or 'GP:A', for keyboard -> '2' (default: 'GP:UpperRightTrigger')
%   'giveFeedback'   : (logical) Give feedback if option is on (default: true)
%   'isDemo'         : (logical) Data won't be saved if on (default: false)
%   'subjectName'    : (string)  Name of subject (default: 'AN')
%
% History:
%   06/07/21  amn  Adapted from BrainardLab/VirtualWorldPsychophysics

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.addParameter('nIterations', 15, @isscalar);
parser.addParameter('controlSignal', 'gamePad', @ischar);
parser.addParameter('option1Key', 'GP:UpperLeftTrigger', @ischar);
parser.addParameter('option2Key', 'GP:UpperRightTrigger', @ischar);
parser.addParameter('giveFeedback', true, @islogical);
parser.addParameter('isDemo', false, @islogical);
parser.addParameter('subjectName', 'AN', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
nIterations    = parser.Results.nIterations;
controlSignal  = parser.Results.controlSignal;
option1Key     = parser.Results.option1Key;
option2Key     = parser.Results.option2Key;
giveFeedback   = parser.Results.giveFeedback;
isDemo         = parser.Results.isDemo;
subjectName    = parser.Results.subjectName;

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
params.fpColor      = [34 70 34]/255; % fixation point color
params.fpColorRed   = [0.6 0.2 0.2]; % fixation point color red
params.bgColor      = [128 128 128]/255; % to match electrophys task
params.textColor    = [0.6 0.2 0.2];
params.image1Loc  = [0 0];
params.image2Loc  = [0 0];
params.image1Size = [10.54 10.54]; % monitor distance=75cm: scene 8 deg vis angle (target 4 deg)
params.image2Size = [10.54 10.54];
params.ISI = 1; % seconds
params.ITI = 0.25; % seconds
params.stimDuration = 0.25; % seconds
params.option1Key = option1Key;
params.option2Key = option2Key;

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
%               11 comparisons per condition * 2 conditions = 22 total comparisons per block.
%
% BLOCK       : Per block, trials will be interspersed from 2 conditions. 
%
% ITERATION   : 'nIterations' (default 15) of each block.
%               Run a complete block (22 total comparisons, pseudorandom order) 
%               before moving on to next block.
%               22 total comparisons per block * 15 iterations of each block = 330 trials.
%
% NOISE LEVEL : 330 trials make up 1 noise level.
%               Noise is in a task-irrelevant feature/object.
%               Each noise level is a pool with a Gaussian distribution,
%               mean of 0, and a set standard deviation.
%               Per trial, take a random draw from the pool.
%               15 iterations of a comparison will be made up of 15 random draws.
%
% SESSION     : Per session, 3 noise levels (0, 1, and 2).
%               Order of noise levels will be randomly assigned per session.
%               330 trials per noise level * 3 noise levels = 990 trials.
%
% RUN         : Each session will be divided into 4 runs.
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
    p1   = strfind(name,'center');
    p2   = strfind(name,'_');
    p2a  = p2(1);
    p3   = strfind(name,'comp');
    p2b  = p2(2);
    p4   = strfind(name,'noise');
    p2c  = p2(3);
    p5   = strfind(name,'.mat');
    
    % Label the 'imageNames' indices by their condition, comparison amount,
    % noise level, and noise amount.
    imageCondition(ii)   = sscanf(name(p1+6:p2a-3),'%f');
    imageComparison(ii)  = sscanf(name(p3+4:p2b-3),'%f');
    imageNoiseLevel(ii)  = sscanf(name(p4+5:p2c-1),'%f');
    imageNoiseAmount(ii) = sscanf(name(p2c+1:p5-4),'%f');
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
nTrialsNoise = nTrialsBlock*nIterations;
nTrials      = nTrialsNoise*nNoiseLevels;
trialOrder   = nan(nTrials,2);
srow         = 1;

% Preallocate vectors of noise level and noise amount (within the level) per trial.
trialNoiseLevel  = nan(nTrials,1);
trialNoiseAmount = nan(nTrials,1);

% Create a trial order per noise level, then combine across noise levels.
for jj = 1:nNoiseLevels
    trialOrderNoise       = nan(nTrialsNoise,2);
    trialNoiseAmountNoise = nan(nTrialsNoise,1);
    nrow = 1;
    
    % Per iteration of a block, create a randomized order of trials within the block.
    for ii = 1:nIterations
        trialsPerBlock       = nan(nTrialsBlock,2);
        noiseAmountsPerBlock = nan(nTrialsBlock,1);
        brow = 1;
        
        % First create an ordered list of image index comparisons for trials of a single block.
        for jjj = 1:nConditions
            
            % Find 'imageNames' index for the image with the center position for this condition.
            % This reference image will always be from Noise Level 0.
            centerpos = conditions(jjj);
            centerIdx = find(imageNoiseLevel==0 & imageCondition==centerpos & imageComparison==0);
            
            % Create ordered list of image index comparisons and noise amounts
            % for trials for this condition.
            trialsPerCondition      = nan(nComparisons,2);
            noiseAmountPerCondition = nan(nComparisons,1);
            for iii = 1:nComparisons
                % For each comparison amount, create pool of images (of the various
                % noise amounts for this noise level, including the noise amount of
                % 0 from Noise Level 0) for this condition.
                imagePool = find(imageComparison==comparisons(iii) & ...
                    (imageNoiseLevel==noiseLevels(jj) | imageNoiseLevel==0) & ...
                    imageCondition==centerpos);
                
                % Randomly draw the comparison image index from the above pool.
                thisImage = imagePool(randi(numel(imagePool)));
                
                % Randomize whether the center position is shown in the 1st/2nd interval.
                thisIndices = [centerIdx thisImage];
                trialsPerCondition(iii,:)    = thisIndices(randperm(2));
                noiseAmountPerCondition(iii) = imageNoiseAmount(thisImage);
            end
            
            % Combine ordered trials across all conditions for a single block.
            trialsPerBlock      (brow:brow+nComparisons-1,:) = trialsPerCondition;
            noiseAmountsPerBlock(brow:brow+nComparisons-1,1) = noiseAmountPerCondition;
            brow = brow+nComparisons;
        end
        
        % Randomize trial order within this single block.
        trialsrand = randperm(nTrialsBlock);
        trialsPerBlockRand       = trialsPerBlock      (trialsrand,:);
        noiseAmountsPerBlockRand = noiseAmountsPerBlock(trialsrand,1);
        
        % Combine iterations of blocks for this noise level.
        trialOrderNoise      (nrow:nrow+nTrialsBlock-1,:) = trialsPerBlockRand;
        trialNoiseAmountNoise(nrow:nrow+nTrialsBlock-1,1) = noiseAmountsPerBlockRand;
        nrow = nrow+nTrialsBlock;
    end
    
    % Combine trials across noise levels.
    trialOrder      (srow:srow+nTrialsNoise-1,:) = trialOrderNoise;
    trialNoiseLevel (srow:srow+nTrialsNoise-1,1) = noiseLevels(jj);
    trialNoiseAmount(srow:srow+nTrialsNoise-1,1) = trialNoiseAmountNoise;
    srow = srow+nTrialsNoise;
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

%% Set up vectors to keep track of trial info
%
% Set up vector for subject response per trial.
selectedResponse = nan(nTrials,1);

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
win.enableObject('fp');
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

% Wait the duration of the intertrial interval.
mglWaitSecs(params.ITI);

%% Per trial, present images and wait for key press response
%
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
    
    % Write the images into the window and disable.
    win.addImage(params.image1Loc, params.image1Size, image1, 'Name', 'image1');
    win.addImage(params.image2Loc, params.image2Size, image2, 'Name', 'image2');
    win.disableObject('image1');
    win.disableObject('image2');
    
    % Enable 1st image and draw.
    win.enableObject('image1');
    win.draw;
    
    % Wait for stimulus duration.
    mglWaitSecs(params.stimDuration);
    win.disableObject('image1');
    win.draw;
    
    % Wait for ISI and show 2nd image.
    mglWaitSecs(params.ISI);
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
    
    % Check if one quarter of experiment is reached.
    if iiTrial == ceil(nTrials/4)
        win.enableObject('oneQuarterText');
        win.disableObject('fp');        
        win.enableObject('fpRed');
        win.draw;
        pause(60);
        win.disableObject('oneQuarterText');
        win.enableObject('restOver');
        win.disableObject('fpRed');
        win.enableObject('fp');
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
        % Reset the keyboard queue.
        mglGetKeyEvent;
        % Wait the duration of the intertrial interval.
        mglWaitSecs(params.ITI);
    end
    
    % Check if two quarters of experiment is reached.
    if iiTrial == ceil(2*nTrials/4)
        win.enableObject('twoQuartersText');
        win.disableObject('fp');        
        win.enableObject('fpRed');
        win.draw;
        pause(60);
        win.disableObject('twoQuartersText');
        win.enableObject('restOver');
        win.disableObject('fpRed');
        win.enableObject('fp');
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
        % Reset the keyboard queue.
        mglGetKeyEvent;
        % Wait the duration of the intertrial interval.
        mglWaitSecs(params.ITI);
    end
    
    % Check if three quarters of experiment is reached.
    if iiTrial == ceil(3*nTrials/4)
        win.enableObject('threeQuartersText');
        win.disableObject('fp');        
        win.enableObject('fpRed');
        win.draw;
        pause(60);
        win.disableObject('threeQuartersText');
        win.enableObject('restOver');
        win.disableObject('fpRed');
        win.enableObject('fp');
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
        % Reset the keyboard queue.
        mglGetKeyEvent;
        % Wait the duration of the intertrial interval.
        mglWaitSecs(params.ITI);
    end
        
    % Check if end of experiment is reached.
    if iiTrial == nTrials
        keepLooping = false;
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
    
    % Add text for when experiment is one quarter over.
    win.addText('One quarter of trials complete. Take a 1 minute break.', ... % Text to display
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
    win.addText('Half of trials complete. Take a 1 minute break.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'twoQuartersText'); % Identifier for the object
    
    % Add text for when experiment is three quarters over.
    win.addText('Three quarters of trials complete. Take a 1 minute break.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'threeQuartersText'); % Identifier for the object

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

%% End