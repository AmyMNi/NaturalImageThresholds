function acquisitionStatus = runNaturalImageExperiment(varargin)
%runNaturalImageExperiment
%
% Usage:
%   runNaturalImageExperiment();
%
% Description:
%   Run the natural image psychophysics experiment given the specified
%   experiment folder. Save the experimental parameters and psychophysical 
%   responses (struct 'data') in the specified output folder.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'nPresentations' : (scalar)  Number of presentations per image comparison (default: 5)
%   'controlSignal'  : (string)  Input method for user response (options: 'keyboard', 'gamePad') (default: 'keyboard')
%   'option1Key'     : (string)  For keyboard -> '1', for gamePad either 'GP:UpperLeftTrigger'  or 'GP:X' (default: '1')
%   'option2Key'     : (string)  For keyboard -> '2', For gamePad either 'GP:UpperRightTrigger' or 'GP:A' (default: '2')
%   'giveFeedback'   : (logical) Give feedback if option is on (default: true)
%   'isDemo'         : (logical) Data won't be saved if on (default: false)
%   'subjectName'    : (string)  Name of subject (default: 'testSubject')
%
% History:
%   06/07/21  amn  Adapted from BrainardLab/VirtualWorldPsychophysics

%% Parse the inputs
parser = inputParser();
parser.addParameter('experimentName', 'Experiment000', @ischar);
parser.addParameter('nPresentations', 5, @isscalar);
parser.addParameter('controlSignal', 'keyboard', @ischar);
parser.addParameter('option1Key', '1', @ischar);
parser.addParameter('option2Key', '2', @ischar);
parser.addParameter('giveFeedback', true, @islogical);
parser.addParameter('isDemo', false, @islogical);
parser.addParameter('subjectName', 'testSubject', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;
nPresentations = parser.Results.nPresentations;
controlSignal  = parser.Results.controlSignal;
option1Key     = parser.Results.option1Key;
option2Key     = parser.Results.option2Key;
giveFeedback   = parser.Results.giveFeedback;
isDemo         = parser.Results.isDemo;
subjectName    = parser.Results.subjectName;

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

% Get the image file names. Images will be called by their index here.
imageNames = {fileInfo(:).name}';

% Get the number of images in the folder.
nImages = numel(imageNames);

%% Create a random trial order for this experiment
%
% Make a list of all image comparison pairs.
% Each image is compared to each other image.
%   col 1: image index shown in 1st interval
%   col 2: image index shown in 2nd interval
comparisons = nchoosek(1:nImages,2);

% Also include images shown in the opposite order to the above.
comparisonsFlip = flip(comparisons,2);

% Finally, each image is also compared to itself.
comparisonI = (1:nImages)';
comparisonI = [comparisonI comparisonI];

% Combine all image comparison pairs into one list.
comparisons = [comparisons; comparisonsFlip; comparisonI];

% For the specified number of presentations per image comparison (blocks),
% create a random trial order for the above image comparison pairs.
numcomp = size(comparisons,1);
trialOrder = nan(numcomp*nPresentations,2);
for ii = 1:nPresentations
    indexorder = randperm(numcomp);
    trialOrderthis = comparisons(indexorder,:);
    trialOrder((ii-1)*numcomp+1:ii*numcomp,:) = trialOrderthis;
end

% As is, one full set of image comparison pairs (a block) will be presented
% before a next block is presented. This evenly distributes the tested 
% pairs across time. However, this may also allow the subject to keep track 
% of the remaining comparisons in a block.
% To avoid this, randomly shuffle the trial order within the full
% experiment trial list (all blocks).
numTrials = size(trialOrder,1);
indexorder = randperm(numTrials);
trialOrder = trialOrder(indexorder,:);

%% Calculate the correct response for each trial
%
% The correct response is 1 if the 2nd target is to the left of the 1st target.
% The correct response is 2 if the 2nd target is to the right of the 1st target.
% The images are indexed in left-right order, with image index 1 the most left.
trialDiff = diff(trialOrder,1,2);
correctResponse = nan(numTrials,1);
correctResponse(trialDiff <0) = 1; % 2nd target is to the left
correctResponse(trialDiff >0) = 2; % 2nd target is to the right
correctResponse(trialDiff==0) = randi(2); % no difference: random assignment

%% Set up vector to keep track of subject response for each trial
selectedResponse = nan(numTrials,1);

%% Begin task

% Note start time of experiment now.
startTime = datestr(now);

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
                    case {'q'}
                        fprintf('Do you want to quit? Type Y for Yes, otherwise give your response \n');
                        key2 = [];
                        while isempty(key2)
                            key2 = mglGetKeyEvent;
                            if ~isempty(key2)
                                switch key2.charCode
                                    case {option1Key,option2Key}
                                        selectedResponse(iiTrial) = getUserResponse(params,key2);
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
                    otherwise
                        key = [];
                end
            end
            pressedKeyboard = mglGetKeyEvent;
            if ~isempty(pressedKeyboard)
                switch pressedKeyboard.charCode
                    case {'q'}
                        fprintf('Do you want to quit? Type Y for Yes, otherwise give your response using gamepad \n');
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
        fprintf('Quitting without saving any data.\n');
        saveData = 0;
    end
    
    % Check if one third of experiment is reached.
    if iiTrial == ceil(numTrials/3)
        win.enableObject('oneThirdText');
        win.disableObject('fp');        
        win.enableObject('fpRed');
        win.draw;
        pause(60);
        win.disableObject('oneThirdText');
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
    end
    
    % Check if two thirds of experiment is reached.
    if iiTrial == ceil(2*numTrials/3)
        win.enableObject('twoThirdText');
        win.disableObject('fp');        
        win.enableObject('fpRed');
        win.draw;
        pause(60);
        win.disableObject('twoThirdText');
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
    end
        
    % Check if end of experiment is reached.
    if iiTrial == numTrials
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
endTime = datestr(now);

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
    data.nPresentations = nPresentations;
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
    data.imageNames = imageNames;
    
    % Save image indices for each trial.
    data.trialOrder = trialOrder;

    % Save correct response for each trial.
    data.correctResponse = correctResponse;
    
    % Save selected response for each trial.
    data.selectedResponse = selectedResponse;
    
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
    
    % Add text for when experiment is one third over.
    win.addText('One third of trials complete. Take a 1 minute break.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'oneThirdText'); % Identifier for the object

    % Add text for when rest period is over.
    win.addText('Rest time complete. Hit any button to continue.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'restOver'); % Identifier for the object

    % Add text for when experiment is two thirds over.
    win.addText('Two thirds of trials complete. Take a 1 minute break.', ... % Text to display
        'Center', [0 8], ... % Where to center the text (x,y)
        'FontSize', 75, ... % Font size
        'Color', params.textColor, ... % RGB color
        'Name', 'twoThirdText'); % Identifier for the object

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