function dataAnalysis = analyzeElectrophysData(varargin)
%analyzeElectrophysData
%
% Usage:
%   dataAnalysis = analyzeElectrophysData('dataDate','211013','imageSet',1)
%
% Description:
%   Analyze electrophysiological data from a single Sdata collection session.
%   Save the results in the specified output folder.
%
% Optional parameters/values:
%   'dataDate'    : (string)  Electrophys data file name (default: '211013')
%   'imageSet'    : (scalar)  Image set presented (default: 1)
%   'plotFigures' : (logical) Plot figures if option is on (default: true)
%   'saveData'    : (logical) Save data if option is on (default: true)
%
% History:
%   10/13/21  amn  Wrote it.

%% Record of electrophys data file name & corresponding image set presented

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ('dataDate','211013','imageSet',1) : sRGB image set to test position/rotation/depth ranges



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse the inputs
parser = inputParser();
parser.addParameter('dataDate', '211013', @ischar);
parser.addParameter('imageSet', 1, @isscalar);
parser.addParameter('plotFigures', true, @islogical);
parser.addParameter('saveData', true, @islogical);
parser.parse(varargin{:});

dataDate    = parser.Results.dataDate;
imageSet    = parser.Results.imageSet;
plotFigures = parser.Results.plotFigures;
saveData    = parser.Results.saveData;

%% Set paths to input and output files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO:
% Data on Raptor. Example: ram/data/main/ivory/211013/ivory_map_211013-133559_dense.mat
% Save data file in folder below. Example: 211013data.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set path to data file.
dataFolder = fullfile('/Users','amy','Desktop','Brainard','Natural Image Thresholds','Electrophys Data');
pathToDataFile = fullfile(dataFolder,sprintf('%sdata.mat',dataDate));

% Set path to image info file.
pathToDataInfo = fullfile(dataFolder,sprintf('imgInfo_set%d.mat',imageSet));

% Set path to output file.
pathToOutput = fullfile(dataFolder,sprintf('%sanalysis.mat',dataDate));

%% Load data
%
% Load specified data file.
temp = load(pathToDataFile,'params');    params    = temp.params;    clear temp;
temp = load(pathToDataFile,'eid');       eid       = temp.eid;       clear temp;
temp = load(pathToDataFile,'resp');      resp      = temp.resp;      clear temp;
temp = load(pathToDataFile,'resp_base'); resp_base = temp.resp_base; clear temp;

% Load specified image info file.
temp = load(pathToDataInfo,'imgInfo');   imgInfo   = temp.imgInfo;   clear temp;

%% Electrophys data organization

% Data variables:
 
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

%% Get electrode numbers for each brain area
%
% Get electrode numbers for the V1/V2 array (V1/V2: array number 2).
electrodeV1 = eid(eid(:,1)==2,2);

% Get electrode numbers for the V4 array (V4: array number 1).
electrodeV4 = eid(eid(:,1)==1,2);

%% Get responses for each brain area
%
% Get responses for V1/V2.
V1resp      = resp     (:,ismember(eid(:,2),electrodeV1));
V1resp_base = resp_base(:,ismember(eid(:,2),electrodeV1));

% Get responses for V4.
V4resp      = resp     (:,ismember(eid(:,2),electrodeV4));
V4resp_base = resp_base(:,ismember(eid(:,2),electrodeV4));

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

%% Get responses for each brain area, based on included electrodes
%
% Get included responses for V1/V2.
V1respInc      = resp     (:,ismember(eid(:,2),electrodeV1included));
V1resp_baseInc = resp_base(:,ismember(eid(:,2),electrodeV1included));

% Get included responses for V4.
V4respInc      = resp     (:,ismember(eid(:,2),electrodeV4included));
V4resp_baseInc = resp_base(:,ismember(eid(:,2),electrodeV4included));

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

%% Get unique positions, rotations, and depths
%
% Get all presented central object positions and background object rotations and depths.
positions = unique(imagePosition);
rotations = unique(imageRotation);
depths    = unique(imageDepth);

%% Calculate specific decoder performance for central object (banana) position
%
% Calculate for each combination of positions, averaged across all rotation/depth combos.














%% Save data analysis results

if saveData
    % Save dataAnalysis struct.
    save(pathToOutput,'dataAnalysis');
    fprintf('\nData was saved in:\n%s\n', pathToOutput);
end
end
%% End