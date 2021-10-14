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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO:
% Data on Raptor. Example: ram/data/main/ivory/211013/ivory_map_211013-133559_dense.mat
% Save data file in folder below. Example: 211013data.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Exclude stimuli presented during incomplete trials (fixation broken)
%
% Exclude stimuli for which params.good = 0.
excludeStim = ~[params.good];
params   (excludeStim)   = [];
resp     (excludeStim,:) = [];
resp_base(excludeStim,:) = [];

%% Analyze responses from specific brain areas
%
% Analyze only V4 array responses.
v4resp = resp(:,eid(:,1)==1);
stimNum = [params.num];







%% 

% EXAMPLE RASTER PLOT FOR EACH CHANNEL:
% plot time on x-axis, plot channel number on y-axis
plot(params(1).spikes(:,1), params(1).spikes(:,2), 'k.')




%% Save data analysis results

if saveData
    % Save dataAnalysis struct.
    save(pathToOutput,'dataAnalysis');
    fprintf('\nData was saved in:\n%s\n', pathToOutput);
end
end
%% End