function data = combineElectrophysData(dataNames,dataFolder)
%combineElectrophysData
%
% Usage:
%	data = combineElectrophysData(dataNames,dataFolder)
%
% Description:
%   Combine the data sets described by the input variable 
%   and output the data in a struct.
%
% Inputs:
%   'dataNames'  : (cell) Cell array of data file names
%   'dataFolder' : (char) Data folder containing data files
%
% Output:
%   'data' : (struct) Struct of combined data for analysis
%
% History:
%   10/18/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addRequired('dataNames',@(x)(iscell(x)));
parser.addRequired('dataFolder',@(x)(ischar(x)));
parser.parse(dataNames,dataFolder);

dataNames  = parser.Results.dataNames;
dataFolder = parser.Results.dataFolder;

%% Create struct to hold data
data = struct;

%% Electrophys data & image info organization

% Electrophys data variables:
 
% params : (struct 1 x num stimuli)
%       set: for the Natural Image Thresholds experiment, set = 7
%       num: number of the stimulus (to match to 'imgInfo' table with image info)
%       x: position in pixels of the stimulus
%       y: position in pixels of the stimulus
%       s: size in pixels of the stimulus
%       good: 0 = exclude (broke fixation), 1 = include (trial was completed)
%       lims_trial: not relevant to this experiment
%       lims_stim: not relevant to this experiment
%       eyes: position of eye in x and y coordinates in microvolts
%           column 1: x position
%           column 2: y position
%       spikes: spike trains for a -100 from stim onset : +100 from stim offset
%           column 1: spike times in milliseconds, indexed to 0 from the stimulus onset
%           column 2: channel spike was on

% eid : (256 total electrodes x 2 columns)
%       column 1: array number (1 = V4, 2 = V1/V2, 3 = 7A)
%       column 2: electrode number
   
% resp : (num stimuli x 256 total electrodes) stimulus response (Hz)
%       analysis window: stim onset + 50 ms to stim offset + 100 ms
%       has a value of NaN when params.good = 0
    
% resp_base : (num stimuli x 256 total electrodes) baseline response (Hz)
%       analysis window: stim onset - 100 ms to stim onset


% Image info variable:

% imgInfo : (table num images x 4 columns)
%       column 1: imgNumber
%       column 2: bananaPosition
%       column 3: backgroundRotation
%       column 4: backgroundDepth

%% Combine data from each data file
for ii = 1:numel(dataNames)

    % Set path to data file.
    nameThis = dataNames{ii};
    pathToDataFile = fullfile(dataFolder,nameThis);
    
    % Load specified data file.
    temp = load(pathToDataFile,'params');    params    = temp.params;    clear temp;
    temp = load(pathToDataFile,'eid');       eid       = temp.eid;       clear temp;
    temp = load(pathToDataFile,'resp');      resp      = temp.resp;      clear temp;
    temp = load(pathToDataFile,'resp_base'); resp_base = temp.resp_base; clear temp;
    
    % Get image set name.
    p_   = strfind(nameThis,'_');
    pend = strfind(nameThis,'.');
    imageSet = str2double(regexp(nameThis(p_:pend),'[+-]?\d*','Match'));
    
    % Set path to image info file.
    pathToDataInfo = fullfile(dataFolder,sprintf('imgInfo_set%d.mat',imageSet));
    
    % Load specified image info file.
    temp = load(pathToDataInfo,'imgInfo'); imgInfo = temp.imgInfo; clear temp;
    



    
    
    
end

%% Check that all of the stimuli were from the Natural Image Thresholds experiment
%
% Make sure that all params.set values = 7 (the Natural Image Thresholds images).
% If not, print a warning.
experimentSet = [params.set];
if any(experimentSet~=7)
    warning('The stimuli are not all from the Natural Image Thresholds project');
end

%% Exclude stimuli presented during incomplete trials (fixation broken)
%
% Exclude stimuli for which params.good = 0.
excludeStim = ~[params.good];
params   (excludeStim)   = [];
resp     (excludeStim,:) = [];
resp_base(excludeStim,:) = [];

% Number of included stimuli.
numStim = numel(params);

% Save to analysis output struct.
dataAnalysis.numStim = numStim;

%% Get electrode numbers for each brain area
%
% Get electrode numbers for the V1/V2 array (V1/V2: array number 2).
electrodeV1 = eid(eid(:,1)==2,2);

% Get electrode numbers for the V4 array (V4: array number 1).
electrodeV4 = eid(eid(:,1)==1,2);

% Save to analysis output struct.
dataAnalysis.electrodeV1 = electrodeV1;
dataAnalysis.electrodeV4 = electrodeV4;

%% Get responses for each brain area
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: responses currently based on same analysis period
%       for both V1/V2 and V4 (50-350 ms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get responses for V1/V2.
V1resp      = resp     (:,ismember(eid(:,2),electrodeV1));
V1resp_base = resp_base(:,ismember(eid(:,2),electrodeV1));

% Get responses for V4.
V4resp      = resp     (:,ismember(eid(:,2),electrodeV4));
V4resp_base = resp_base(:,ismember(eid(:,2),electrodeV4));

% Save to analysis output struct.
dataAnalysis.V1resp = V1resp;
dataAnalysis.V1resp_base = V1resp_base;
dataAnalysis.V4resp = V4resp;
dataAnalysis.V4resp_base = V4resp_base;

%% Determine electrodes to include per brain area, based on stimulus-evoked firing rate
%
% Base electrode inclusion on the following criteria:
%   stimulus-evoked firing rate statsig greater than baseline firing rate,
%   minimum stimulus-evoked firing rate.
frMin = 10; %Hz
alphaLevel = .05;
electrodeV1included = false(1,size(V1resp,2));
electrodeV4included = false(1,size(V4resp,2));

% Calculate included electrodes: V1.
for ii = 1:size(V1resp,2)
    % Two-sample one-tailed t-test.
    h = ttest2(V1resp(:,ii), V1resp_base(:,ii), 'Tail', 'right', 'Alpha', alphaLevel);
    if h==1 && nanmean(V1resp(:,ii)) > frMin
        electrodeV1included(ii) = true;
    end
end

% Calculate included electrodes: V4.
for ii = 1:size(V4resp,2)
    % Two-sample one-tailed t-test.
    h = ttest2(V4resp(:,ii), V4resp_base(:,ii), 'Tail', 'right', 'Alpha', alphaLevel);
    if h==1 && nanmean(V4resp(:,ii)) > frMin
        electrodeV4included(ii) = true;
    end
end

electrodeV1included = electrodeV1(electrodeV1included');
electrodeV4included = electrodeV4(electrodeV4included');

% Save to analysis output struct.
dataAnalysis.electrodeV1included = electrodeV1included;
dataAnalysis.electrodeV4included = electrodeV4included;

%% Get responses for each brain area, based on included electrodes
%         
% Get included responses for V1/V2.
V1respInc      = resp     (:,ismember(eid(:,2),electrodeV1included));
V1resp_baseInc = resp_base(:,ismember(eid(:,2),electrodeV1included));

% Get included responses for V4.
V4respInc      = resp     (:,ismember(eid(:,2),electrodeV4included));
V4resp_baseInc = resp_base(:,ismember(eid(:,2),electrodeV4included));

% Save to analysis output struct.
dataAnalysis.V1respInc = V1respInc;
dataAnalysis.V1resp_baseInc = V1resp_baseInc;
dataAnalysis.V4respInc = V4respInc;
dataAnalysis.V4resp_baseInc = V4resp_baseInc;

%% Per stimulus, get image number and info
%
% Per stimulus, get image number to match to imgInfo table with image info.
imageNum = [params.num]';

% Convert imgInfo table to matrix.
% columns:  imgNumber  bananaPosition  backgroundRotation  backgroundDepth
imgInfo = table2array(imgInfo);

% Per stimulus, get central object (banana) position.
[~,indImageNum] = ismember(imageNum, imgInfo(:,1));
imagePosition = imgInfo(indImageNum,2);

% Per stimulus, get background object (branches & leaves) rotation.
imageRotation = imgInfo(indImageNum,3);

% Per stimulus, get background object (branches & leaves) depth.
imageDepth = imgInfo(indImageNum,4);

% Save to analysis output struct.
dataAnalysis.imageNum = imageNum;
dataAnalysis.imagePosition = imagePosition;
dataAnalysis.imageRotation = imageRotation;
dataAnalysis.imageDepth = imageDepth;

%% Get unique positions, rotations, and depths
%
% Get unique central object positions and background object rotations and depths.
positions = unique(imagePosition);
rotations = unique(imageRotation);
depths    = unique(imageDepth);

% Save to analysis output struct.
dataAnalysis.positions = positions;
dataAnalysis.rotations = rotations;
dataAnalysis.depths = depths;

%% Save data struct
if saveData
    save(pathToOutput,'data');
    fprintf('\ndata was saved in:\n%s\n', pathToOutput);
end

end

%% Helper functions

%%

%% End