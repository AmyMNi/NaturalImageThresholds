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

%% Create struct to hold all analysis results
dataAnalysis = struct;

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

% Save to analysis output struct.
dataAnalysis.pathToDataFile = pathToDataFile;
dataAnalysis.pathToDataInfo = pathToDataInfo;

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

%% Calculate specific decoder performance for central object position
%
% Calculate decoder performance for each combination of positions.
[specificPositionV1,specificPositionV4,specificPositionInfo] = helperSpecificDecoder ...
    (positions,rotations,depths,imagePosition,imageRotation,imageDepth,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.specificPositionV1 = specificPositionV1;
dataAnalysis.specificPositionV4 = specificPositionV4;
dataAnalysis.specificPositionInfo = specificPositionInfo;


%% Plot the mean performance for specific decoders of position, per size of discriminated position difference
%
% Separate the specific decoders into groups based on the size of their 
% discriminated position difference, then calculate the mean decoder 
% performance per group.

% Get the unique values of the difference between the two positions
% discriminated by the decoders.
uniquePositionDiff = unique(specificPositionInfo(:,3));

% Calculate mean decoder performance per group.
specificPositionV1mean = nan(numel(uniquePositionDiff),1);
specificPositionV4mean = nan(numel(uniquePositionDiff),1);
for ii = 1:numel(uniquePositionDiff)
    positionDiffThis = uniquePositionDiff(ii);
    
    % Get performances of decoders that discriminated this position difference.
    V1this = specificPositionV1(specificPositionInfo(:,3)==positionDiffThis);
    V4this = specificPositionV4(specificPositionInfo(:,3)==positionDiffThis);
    
    % Calculate mean decoder performance.
    specificPositionV1mean(ii,1) = nanmean(V1this);
    specificPositionV4mean(ii,1) = nanmean(V4this);
end

% Plot the size of discriminated position difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    plot(uniquePositionDiff, specificPositionV1mean, '.-m');
    plot(uniquePositionDiff, specificPositionV4mean, '.-b');
    % Plot parameters.
    title('Specific decoder performance for position');
    legend('V1','V4');
    xlabel('Difference between discriminated positions');
    ylabel('Proportion correct');
    axis([min(uniquePositionDiff)-2 max(uniquePositionDiff)+2 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniquePositionDiff);
    set(gca,'XTickLabel',uniquePositionDiff);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniquePositionDiff = uniquePositionDiff;
dataAnalysis.specificPositionV1mean = specificPositionV1mean;
dataAnalysis.specificPositionV4mean = specificPositionV4mean;

%% Save data analysis results
if saveData
    save(pathToOutput,'dataAnalysis');
    fprintf('\nData was saved in:\n%s\n', pathToOutput);
end

end

%% Helper functions

%% Helper function: Calculate specific decoder performance without noise
%
% Description:
%   Calculate specific decoder performance for discriminating two values of a feature. 
%   No task-irrelevant noise is included in the stimuli.
%
% Inputs:
%   X      : (num values x 1) values of the feature to be discriminated
%   Y      : (num values x 1) values of a first task-irrelevant feature
%   Z      : (num values x 1) values of a second task-irrelevant feature
%   imageX : (num stimuli x 1) Per stimulus, value of the feature to be discriminated
%   imageY : (num stimuli x 1) Per stimulus, value of the first task-irrelevant feature
%   imageZ : (num stimuli x 1) Per stimulus, value of the second task-irrelevant feature 
%   V1resp : (num stimuli x num V1 neurons) Neuronal responses for V1
%   V4resp : (num stimuli x num V4 neurons) Neuronal responses for V4
%
% Outputs:
%   specificDecoderV1 : (num decoders tested x 1) for V1, proportion correct per decoder tested
%   specificDecoderV4 : (num decoders tested x 1) for V4, proportion correct per decoder tested

function [decoderV1,decoderV4,info] = helperSpecificDecoder(X,Y,Z,imageX,imageY,imageZ,V1resp,V4resp)
    
% Calculate the number of combinations of the feature to be discriminated.
numCombos = nchoosek(numel(X),2);

% Set up decoder performance vectors.
decoderV1 = nan(numCombos*numel(Y)*numel(Z),1);
decoderV4 = nan(numCombos*numel(Y)*numel(Z),1);

% Set up vector of the difference between the two values discriminated, per decoder.
info = nan(numCombos*numel(Y)*numel(Z),1);

% Calculate proportion correct per decoder.
row = 0;
for yy = 1:numel(Y)
    % Get a value of the first task-irrelevant feature.
    Ythis = Y(yy);
    
    for zz = 1:numel(Z)
        % Get a value of the second task-irrelevant feature.
        Zthis = Z(zz);
        
        for ii = 1:numel(X)
            % Get Value #1 of the task-relevant feature.
            Xthis1 = X(ii);
                
            for jj = ii+1:numel(X)
                % Get Value #2 of the task-relevant feature.
                Xthis2 = X(jj);
                
                % Record the difference between the two values discriminated.
                row = row+1;
                info(row) = abs(Xthis1-Xthis2);

                % Get the neuronal responses for Value #1 of the task-relevant
                % feature, for this value of the first task-irrelevant feature and
                % this value of the second task-irrelevant feature. 
                stimIdx1 = imageX==Xthis1 & imageY==Ythis & imageZ==Zthis;  
                V1respThis1 = V1resp(stimIdx1,:);
                V4respThis1 = V4resp(stimIdx1,:);
                
                % Same as above but for Value #2 of the task-relevant feature.
                stimIdx2 = imageX==Xthis2 & imageY==Ythis & imageZ==Zthis;  
                V1respThis2 = V1resp(stimIdx2,:);
                V4respThis2 = V4resp(stimIdx2,:);
                
                % Calculate decoder proportion correct on discriminating
                % Value #1 from Value #2.
                [~,V1pc] = calcSpecificDecoder([V1respThis1;V1respThis2],[zeros(size(V1respThis1,1),1);ones(size(V1respThis2,1),1)]);
                [~,V4pc] = calcSpecificDecoder([V4respThis1;V4respThis2],[zeros(size(V4respThis1,1),1);ones(size(V4respThis2,1),1)]);
                decoderV1(row) = V1pc;
                decoderV4(row) = V4pc;
            end
        end
    end
end
end

%% End