function dataAnalysis = analyzeElectrophysData(varargin)
%analyzeElectrophysData
%
% Usage:
%   dataAnalysis = analyzeElectrophysData('dataName','211013','imageSet',1)
%
% Description:
%   Analyze electrophysiological data from a single Sdata collection session.
%   Save the results in the specified output folder.
%
% Optional parameters/values:
%   'dataName'    : (string)  Electrophys data file name (default: '211013-133559data')
%   'imageSet'    : (scalar)  Image set presented (default: 1)
%   'plotFigures' : (logical) Plot figures if option is on (default: true)
%   'saveData'    : (logical) Save data if option is on (default: true)
%
% History:
%   10/13/21  amn  Wrote it.

%% Record of electrophys data file name & corresponding image set presented

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% ('dataName','211013-133559data','imageSet',1) :sRGB image set, 2 iterations/stim
% ('dataName','211014-153026data','imageSet',2) : RGB image set, 5 iterations/stim
% ('dataName','211015-155409data','imageSet',2) : RGB image set, 4 iterations/stim
% ('dataName','211015-160308data','imageSet',3) : RGB image set, 4 iterations/stim





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse the inputs
parser = inputParser();
parser.addParameter('dataName', '211013-133559data', @ischar);
parser.addParameter('imageSet', 1, @isscalar);
parser.addParameter('plotFigures', true, @islogical);
parser.addParameter('saveData', true, @islogical);
parser.parse(varargin{:});

dataName    = parser.Results.dataName;
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
% Data on Raptor. Example: Ram/data/main/ivory/211013/ivory_map_211013-133559_dense.mat
% Save data file in folder below. Example: 211013-133559data.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set path to data file.
dataFolder = fullfile('/Users','amy','Desktop','Brainard','Natural Image Thresholds','Electrophys Data');
pathToDataFile = fullfile(dataFolder,sprintf('%s.mat',dataName));

% Set path to image info file.
pathToDataInfo = fullfile(dataFolder,sprintf('imgInfo_set%d.mat',imageSet));

% Set path to output file.
pathToOutput = fullfile(dataFolder,sprintf('%sAnalysis.mat',dataName));

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

%% Calculate specific decoder performance for central object POSITION
%
% Calculate decoder performance for each combination of positions.
[diffBetweenPositionValues,specificPositionV1,specificPositionV4] = helperSpecificDecoder ...
    (positions,rotations,depths,imagePosition,imageRotation,imageDepth,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.diffBetweenPositionValues = diffBetweenPositionValues;
dataAnalysis.specificPositionV1 = specificPositionV1;
dataAnalysis.specificPositionV4 = specificPositionV4;

%% Plot the mean performance for specific decoders of POSITION, per size of discriminated POSITION difference
%
% Calculate the mean decoder performance per size difference in the discriminated values.
[uniquePositionDiff,specificPositionV1mean,specificPositionV4mean,specificPositionV1sem,specificPositionV4sem] = helperDecoderMean...
    (diffBetweenPositionValues,specificPositionV1,specificPositionV4);

% Plot the size of discriminated position difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    errorbar(uniquePositionDiff, specificPositionV1mean, specificPositionV1sem, '.-m');
    errorbar(uniquePositionDiff, specificPositionV4mean, specificPositionV4sem, '.-b');
    % Plot parameters.
    title('Specific decoder performance for POSITION');
    legend('V1','V4');
    xlabel('Difference between discriminated positions');
    ylabel('Proportion correct');
    axis([min(uniquePositionDiff)-8 max(uniquePositionDiff)+8 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniquePositionDiff);
    set(gca,'XTickLabel',uniquePositionDiff);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniquePositionDiff = uniquePositionDiff;
dataAnalysis.specificPositionV1mean = specificPositionV1mean;
dataAnalysis.specificPositionV4mean = specificPositionV4mean;
dataAnalysis.specificPositionV1sem = specificPositionV1sem;
dataAnalysis.specificPositionV4sem = specificPositionV4sem;

%% Calculate specific decoder performance for background object ROTATION
%
% Calculate decoder performance for each combination of rotations.
[diffBetweenRotationValues,specificRotationV1,specificRotationV4] = helperSpecificDecoder ...
    (rotations,positions,depths,imageRotation,imagePosition,imageDepth,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.diffBetweenRotationValues = diffBetweenRotationValues;
dataAnalysis.specificRotationV1 = specificRotationV1;
dataAnalysis.specificRotationV4 = specificRotationV4;

%% Plot the mean performance for specific decoders of ROTATION, per size of discriminated ROTATION difference
%
% Calculate the mean decoder performance per size difference in the discriminated values.
[uniqueRotationDiff,specificRotationV1mean,specificRotationV4mean,specificRotationV1sem,specificRotationV4sem] = helperDecoderMean...
    (diffBetweenRotationValues,specificRotationV1,specificRotationV4);

% Plot the size of discriminated rotation difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    errorbar(uniqueRotationDiff, specificRotationV1mean, specificRotationV1sem, '.-m');
    errorbar(uniqueRotationDiff, specificRotationV4mean, specificRotationV4sem, '.-b');
    % Plot parameters.
    title('Specific decoder performance for ROTATION');
    legend('V1','V4');
    xlabel('Difference between discriminated rotations');
    ylabel('Proportion correct');
    axis([min(uniqueRotationDiff)-6 max(uniqueRotationDiff)+6 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniqueRotationDiff);
    set(gca,'XTickLabel',uniqueRotationDiff);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniqueRotationDiff = uniqueRotationDiff;
dataAnalysis.specificRotationV1mean = specificRotationV1mean;
dataAnalysis.specificRotationV4mean = specificRotationV4mean;
dataAnalysis.specificRotationV1sem = specificRotationV1sem;
dataAnalysis.specificRotationV4sem = specificRotationV4sem;

%% Calculate specific decoder performance for background object DEPTH
%
% Calculate decoder performance for each combination of depths.
[diffBetweenDepthValues,specificDepthV1,specificDepthV4] = helperSpecificDecoder ...
    (depths,positions,rotations,imageDepth,imagePosition,imageRotation,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.diffBetweenDepthValues = diffBetweenDepthValues;
dataAnalysis.specificDepthV1 = specificDepthV1;
dataAnalysis.specificDepthV4 = specificDepthV4;

%% Plot the mean performance for specific decoders of DEPTH, per size of discriminated DEPTH difference
%
% Calculate the mean decoder performance per size difference in the discriminated values.
[uniqueDepthDiff,specificDepthV1mean,specificDepthV4mean,specificDepthV1sem,specificDepthV4sem] = helperDecoderMean...
    (diffBetweenDepthValues,specificDepthV1,specificDepthV4);

% Plot the size of discriminated depth difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    errorbar(uniqueDepthDiff, specificDepthV1mean, specificDepthV1sem, '.-m');
    errorbar(uniqueDepthDiff, specificDepthV4mean, specificDepthV4sem, '.-b');
    % Plot parameters.
    title('Specific decoder performance for DEPTH');
    legend('V1','V4');
    xlabel('Difference between discriminated depths');
    ylabel('Proportion correct');
    axis([min(uniqueDepthDiff)-63 max(uniqueDepthDiff)+63 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniqueDepthDiff);
    set(gca,'XTickLabel',uniqueDepthDiff);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniqueDepthDiff = uniqueDepthDiff;
dataAnalysis.specificDepthV1mean = specificDepthV1mean;
dataAnalysis.specificDepthV4mean = specificDepthV4mean;
dataAnalysis.specificDepthV1sem = specificDepthV1sem;
dataAnalysis.specificDepthV4sem = specificDepthV4sem;

%% Calculate specific decoder performance for central object POSITION with background object ROTATION included as noise
%
% Calculate decoder performance for each combination of positions,
% for each value of depth.
[diffBetweenPositionValues_NoiseRotation,specificPositionV1_NoiseRotation,specificPositionV4_NoiseRotation] = helperSpecificDecoderNoise1 ...
    (positions,depths,imagePosition,imageDepth,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.diffBetweenPositionValues_NoiseRotation = diffBetweenPositionValues_NoiseRotation;
dataAnalysis.specificPositionV1_NoiseRotation = specificPositionV1_NoiseRotation;
dataAnalysis.specificPositionV4_NoiseRotation = specificPositionV4_NoiseRotation;

%% Plot the mean performance for specific decoders of POSITION, per size of discriminated POSITION difference, with ROTATION noise
%
% Calculate the mean decoder performance per size difference in the discriminated values.
[uniquePositionDiff_NoiseRotation,specificPositionV1mean_NoiseRotation,specificPositionV4mean_NoiseRotation,...
 specificPositionV1sem_NoiseRotation,specificPositionV4sem_NoiseRotation] = helperDecoderMean...
    (diffBetweenPositionValues_NoiseRotation,specificPositionV1_NoiseRotation,specificPositionV4_NoiseRotation);

% Plot the size of discriminated position difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    errorbar(uniquePositionDiff_NoiseRotation, specificPositionV1mean_NoiseRotation, specificPositionV1sem_NoiseRotation, '.-m');
    errorbar(uniquePositionDiff_NoiseRotation, specificPositionV4mean_NoiseRotation, specificPositionV4sem_NoiseRotation, '.-b');
    % Plot parameters.
    title('Specific decoder performance for POSITION with ROTATION noise');
    legend('V1','V4');
    xlabel('Difference between discriminated positions');
    ylabel('Proportion correct');
    axis([min(uniquePositionDiff_NoiseRotation)-3 max(uniquePositionDiff_NoiseRotation)+3 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniquePositionDiff_NoiseRotation);
    set(gca,'XTickLabel',uniquePositionDiff_NoiseRotation);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniquePositionDiff_NoiseRotation = uniquePositionDiff_NoiseRotation;
dataAnalysis.specificPositionV1mean_NoiseRotation = specificPositionV1mean_NoiseRotation;
dataAnalysis.specificPositionV4mean_NoiseRotation = specificPositionV4mean_NoiseRotation;
dataAnalysis.specificPositionV1sem_NoiseRotation = specificPositionV1sem_NoiseRotation;
dataAnalysis.specificPositionV4sem_NoiseRotation = specificPositionV4sem_NoiseRotation;

%% Calculate specific decoder performance for central object POSITION with background object DEPTH included as noise
%
% Calculate decoder performance for each combination of positions,
% for each value of rotation.
[diffBetweenPositionValues_NoiseDepth,specificPositionV1_NoiseDepth,specificPositionV4_NoiseDepth] = helperSpecificDecoderNoise1 ...
    (positions,rotations,imagePosition,imageRotation,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.diffBetweenPositionValues_NoiseDepth = diffBetweenPositionValues_NoiseDepth;
dataAnalysis.specificPositionV1_NoiseDepth = specificPositionV1_NoiseDepth;
dataAnalysis.specificPositionV4_NoiseDepth = specificPositionV4_NoiseDepth;

%% Plot the mean performance for specific decoders of POSITION, per size of discriminated POSITION difference, with DEPTH noise
%
% Calculate the mean decoder performance per size difference in the discriminated values.
[uniquePositionDiff_NoiseDepth,specificPositionV1mean_NoiseDepth,specificPositionV4mean_NoiseDepth,...
 specificPositionV1sem_NoiseDepth,specificPositionV4sem_NoiseDepth] = helperDecoderMean...
    (diffBetweenPositionValues_NoiseDepth,specificPositionV1_NoiseDepth,specificPositionV4_NoiseDepth);

% Plot the size of discriminated position difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    errorbar(uniquePositionDiff_NoiseDepth, specificPositionV1mean_NoiseDepth, specificPositionV1sem_NoiseDepth, '.-m');
    errorbar(uniquePositionDiff_NoiseDepth, specificPositionV4mean_NoiseDepth, specificPositionV4sem_NoiseDepth, '.-b');
    % Plot parameters.
    title('Specific decoder performance for POSITION with DEPTH noise');
    legend('V1','V4');
    xlabel('Difference between discriminated positions');
    ylabel('Proportion correct');
    axis([min(uniquePositionDiff_NoiseDepth)-3 max(uniquePositionDiff_NoiseDepth)+3 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniquePositionDiff_NoiseDepth);
    set(gca,'XTickLabel',uniquePositionDiff_NoiseDepth);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniquePositionDiff_NoiseDepth = uniquePositionDiff_NoiseDepth;
dataAnalysis.specificPositionV1mean_NoiseDepth = specificPositionV1mean_NoiseDepth;
dataAnalysis.specificPositionV4mean_NoiseDepth = specificPositionV4mean_NoiseDepth;
dataAnalysis.specificPositionV1sem_NoiseDepth = specificPositionV1sem_NoiseDepth;
dataAnalysis.specificPositionV4sem_NoiseDepth = specificPositionV4sem_NoiseDepth;

%% Calculate specific decoder performance for central object POSITION with background object ROTATION and NOISE included as noise
%
% Calculate decoder performance for each combination of positions.
[diffBetweenPositionValues_NoiseRotationDepth,specificPositionV1_NoiseRotationDepth,specificPositionV4_NoiseRotationDepth] = helperSpecificDecoderNoise2 ...
    (positions,imagePosition,V1respInc,V4respInc);

% Save to analysis output struct.
dataAnalysis.diffBetweenPositionValues_NoiseRotationDepth = diffBetweenPositionValues_NoiseRotationDepth;
dataAnalysis.specificPositionV1_NoiseRotationDepth = specificPositionV1_NoiseRotationDepth;
dataAnalysis.specificPositionV4_NoiseRotationDepth = specificPositionV4_NoiseRotationDepth;

%% Plot the mean performance for specific decoders of POSITION, per size of discriminated POSITION difference, with ROTATION and DEPTH noise
%
% Calculate the mean decoder performance per size difference in the discriminated values.
[uniquePositionDiff_NoiseRotationDepth,specificPositionV1mean_NoiseRotationDepth,specificPositionV4mean_NoiseRotationDepth,...
 specificPositionV1sem_NoiseRotationDepth,specificPositionV4sem_NoiseRotationDepth] = helperDecoderMean...
    (diffBetweenPositionValues_NoiseRotationDepth,specificPositionV1_NoiseRotationDepth,specificPositionV4_NoiseRotationDepth);

% Plot the size of discriminated position difference on the x-axis
% and the mean decoder proportion correct on the y-axis.
if plotFigures
    figure; hold on; axis square;
    errorbar(uniquePositionDiff_NoiseRotationDepth, specificPositionV1mean_NoiseRotationDepth, specificPositionV1sem_NoiseRotationDepth, '.-m');
    errorbar(uniquePositionDiff_NoiseRotationDepth, specificPositionV4mean_NoiseRotationDepth, specificPositionV4sem_NoiseRotationDepth, '.-b');
    % Plot parameters.
    title('Specific decoder performance for POSITION with ROTATION & DEPTH noise');
    legend('V1','V4');
    xlabel('Difference between discriminated positions');
    ylabel('Proportion correct');
    axis([min(uniquePositionDiff_NoiseRotationDepth)-3 max(uniquePositionDiff_NoiseRotationDepth)+3 0 1]);
    set(gca,'tickdir','out');
    set(gca,'XTick',uniquePositionDiff_NoiseRotationDepth);
    set(gca,'XTickLabel',uniquePositionDiff_NoiseRotationDepth);
    box off; hold off;
end

% Save to analysis output struct.
dataAnalysis.uniquePositionDiff_NoiseRotationDepth = uniquePositionDiff_NoiseRotationDepth;
dataAnalysis.specificPositionV1mean_NoiseRotationDepth = specificPositionV1mean_NoiseRotationDepth;
dataAnalysis.specificPositionV4mean_NoiseRotationDepth = specificPositionV4mean_NoiseRotationDepth;
dataAnalysis.specificPositionV1sem_NoiseRotationDepth = specificPositionV1sem_NoiseRotationDepth;
dataAnalysis.specificPositionV4sem_NoiseRotationDepth = specificPositionV4sem_NoiseRotationDepth;

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
%   imageX : (num stimuli x 1) per stimulus, value of the feature to be discriminated
%   imageY : (num stimuli x 1) per stimulus, value of the first task-irrelevant feature
%   imageZ : (num stimuli x 1) per stimulus, value of the second task-irrelevant feature 
%   V1resp : (num stimuli x num V1 neurons) neuronal responses for V1
%   V4resp : (num stimuli x num V4 neurons) neuronal responses for V4
%
% Outputs:
%   diffBetweenValues : (num decoders tested x 1) difference between the two values discriminated per decoder
%   decoderV1         : (num decoders tested x 1) for V1, proportion correct per decoder
%   decoderV4         : (num decoders tested x 1) for V4, proportion correct per decoder

function [diffBetweenValues,decoderV1,decoderV4] = helperSpecificDecoder(X,Y,Z,imageX,imageY,imageZ,V1resp,V4resp)
    
% Calculate the number of combinations of the feature to be discriminated.
numCombos = nchoosek(numel(X),2);

% Set up vectors of decoder performance 
% and of the difference between the two values discriminated, per decoder.
decoderV1 = nan(numCombos*numel(Y)*numel(Z),1);
decoderV4 = nan(numCombos*numel(Y)*numel(Z),1);
diffBetweenValues = nan(numCombos*numel(Y)*numel(Z),1);

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
                diffBetweenValues(row) = abs(Xthis1-Xthis2);

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
                
                % Calculate decoder proportion correct on discriminating Value #1 from Value #2.
                [~,V1pc] = calcSpecificDecoder([V1respThis1;V1respThis2],[zeros(size(V1respThis1,1),1);ones(size(V1respThis2,1),1)]);
                [~,V4pc] = calcSpecificDecoder([V4respThis1;V4respThis2],[zeros(size(V4respThis1,1),1);ones(size(V4respThis2,1),1)]);
                decoderV1(row) = V1pc;
                decoderV4(row) = V4pc;
            end
        end
    end
end
end

%% Helper function: Calculate specific decoder performance with 1 task-irrelevant feature included as noise
%
% Description:
%   Calculate specific decoder performance for discriminating two values of a feature. 
%   One task-irrelevant feature is included as noise in the stimuli,
%   while the other task-irrelevant feature is NOT included as noise, and 
%   is held constant, and separate decoders are calculated for each value 
%   of this non-included feature.
%
% Inputs:
%   X        : (num values x 1) values of the feature to be discriminated
%   Y        : (num values x 1) values of the task-irrelevant feature that is NOT included as noise
%   imageX   : (num stimuli x 1) per stimulus, value of the feature to be discriminated
%   imageY   : (num stimuli x 1) per stimulus, value of the task-irrelevant feature that is NOT included as noise
%   V1resp   : (num stimuli x num V1 neurons) neuronal responses for V1
%   V4resp   : (num stimuli x num V4 neurons) neuronal responses for V4
%
% Outputs:
%   diffBetweenValues : (num decoders tested x 1) difference between the two values discriminated per decoder
%   decoderV1         : (num decoders tested x 1) for V1, proportion correct per decoder
%   decoderV4         : (num decoders tested x 1) for V4, proportion correct per decoder

function [diffBetweenValues,decoderV1,decoderV4] = helperSpecificDecoderNoise1(X,Y,imageX,imageY,V1resp,V4resp)
    
% Calculate the number of combinations of the feature to be discriminated.
numCombos = nchoosek(numel(X),2);

% Set up vectors of decoder performance 
% and of the difference between the two values discriminated, per decoder.
decoderV1         = nan(numCombos*numel(Y),1);
decoderV4         = nan(numCombos*numel(Y),1);
diffBetweenValues = nan(numCombos*numel(Y),1);

% Calculate proportion correct per decoder.
row = 0;
for yy = 1:numel(Y)
    % Get a value of the task-irrelevant feature that will NOT be included as noise,
    % to calculate decoder performance at this value of this task-irrelevant feature.
    Ythis = Y(yy);
    
    for ii = 1:numel(X)
        % Get Value #1 of the task-relevant feature.
        Xthis1 = X(ii);
        
        for jj = ii+1:numel(X)
            % Get Value #2 of the task-relevant feature.
            Xthis2 = X(jj);
            
            % Record the difference between the two values discriminated.
            row = row+1;
            diffBetweenValues(row) = abs(Xthis1-Xthis2);
            
            % Get the neuronal responses for Value #1 of the task-relevant
            % feature, for this value of task-irrelevant feature that is
            % NOT included as noise. 
            stimIdx1 = imageX==Xthis1 & imageY==Ythis;
            V1respThis1 = V1resp(stimIdx1,:);
            V4respThis1 = V4resp(stimIdx1,:);
            
            % Same as above but for Value #2 of the task-relevant feature.
            stimIdx2 = imageX==Xthis2 & imageY==Ythis;
            V1respThis2 = V1resp(stimIdx2,:);
            V4respThis2 = V4resp(stimIdx2,:);
            
            % Calculate decoder proportion correct on discriminating Value #1 from Value #2.
            [~,V1pc] = calcSpecificDecoder([V1respThis1;V1respThis2],[zeros(size(V1respThis1,1),1);ones(size(V1respThis2,1),1)]);
            [~,V4pc] = calcSpecificDecoder([V4respThis1;V4respThis2],[zeros(size(V4respThis1,1),1);ones(size(V4respThis2,1),1)]);
            decoderV1(row) = V1pc;
            decoderV4(row) = V4pc;
        end
    end
end
end

%% Helper function: Calculate specific decoder performance with 2 task-irrelevant features included as noise
%
% Description:
%   Calculate specific decoder performance for discriminating two values of a feature. 
%   Two task-irrelevant features are included as noise in the stimuli.
%
% Inputs:
%   X        : (num values x 1) values of the feature to be discriminated
%   imageX   : (num stimuli x 1) per stimulus, value of the feature to be discriminated
%   V1resp   : (num stimuli x num V1 neurons) neuronal responses for V1
%   V4resp   : (num stimuli x num V4 neurons) neuronal responses for V4
%
% Outputs:
%   diffBetweenValues : (num decoders tested x 1) difference between the two values discriminated per decoder
%   decoderV1         : (num decoders tested x 1) for V1, proportion correct per decoder
%   decoderV4         : (num decoders tested x 1) for V4, proportion correct per decoder

function [diffBetweenValues,decoderV1,decoderV4] = helperSpecificDecoderNoise2(X,imageX,V1resp,V4resp)
    
% Calculate the number of combinations of the feature to be discriminated.
numCombos = nchoosek(numel(X),2);

% Set up vectors of decoder performance 
% and of the difference between the two values discriminated, per decoder.
decoderV1         = nan(numCombos,1);
decoderV4         = nan(numCombos,1);
diffBetweenValues = nan(numCombos,1);

% Calculate proportion correct per decoder.
row = 0;
for ii = 1:numel(X)
    % Get Value #1 of the task-relevant feature.
    Xthis1 = X(ii);
    
    for jj = ii+1:numel(X)
        % Get Value #2 of the task-relevant feature.
        Xthis2 = X(jj);
        
        % Record the difference between the two values discriminated.
        row = row+1;
        diffBetweenValues(row) = abs(Xthis1-Xthis2);
        
        % Get the neuronal responses for Value #1 of the task-relevant feature.
        stimIdx1 = imageX==Xthis1;
        V1respThis1 = V1resp(stimIdx1,:);
        V4respThis1 = V4resp(stimIdx1,:);
        
        % Same as above but for Value #2 of the task-relevant feature.
        stimIdx2 = imageX==Xthis2;
        V1respThis2 = V1resp(stimIdx2,:);
        V4respThis2 = V4resp(stimIdx2,:);
        
        % Calculate decoder proportion correct on discriminating Value #1 from Value #2.
        [~,V1pc] = calcSpecificDecoder([V1respThis1;V1respThis2],[zeros(size(V1respThis1,1),1);ones(size(V1respThis2,1),1)]);
        [~,V4pc] = calcSpecificDecoder([V4respThis1;V4respThis2],[zeros(size(V4respThis1,1),1);ones(size(V4respThis2,1),1)]);
        decoderV1(row) = V1pc;
        decoderV4(row) = V4pc;
    end
end
end

%% Helper function: Calculate mean decoder performance per size difference in discriminated values
%
% Description:
%	Separate the decoders into groups based each decoder's difference in 
%   size between the two discriminated values. Then calculate the mean
%   decoder performance per group.
%
% Inputs:
%   diffBetweenValues : (num decoders tested x 1) difference between the two values discriminated per decoder
%   decoderV1         : (num decoders tested x 1) for V1, proportion correct per decoder
%   decoderV4         : (num decoders tested x 1) for V4, proportion correct per decoder
%
% Outputs:
%   valueDiff     : (num x 1) unique values of the size difference in the values discriminated by the decoder 
%   decoderV1mean : (num x 1) for V1, mean decoder performance per size difference
%   decoderV4mean : (num x 1) for V4, mean decoder performance per size difference
%   decoderV1sem  : (num x 1) for V1, standard error of the mean (SEM) for mean decoder performance
%   decoderV4sem  : (num x 1) for V4, standard error of the mean (SEM) for mean decoder performance

function [valueDiff,decoderV1mean,decoderV4mean,decoderV1sem,decoderV4sem] = helperDecoderMean(diffBetweenValues,decoderV1,decoderV4)
% Get the unique values of the size difference in the values discriminated by the decoder.
valueDiff = unique(diffBetweenValues);

% Calculate the mean and sem of the decoder performance per group.
decoderV1mean = nan(numel(valueDiff),1);
decoderV4mean = nan(numel(valueDiff),1);
decoderV1sem  = nan(numel(valueDiff),1);
decoderV4sem  = nan(numel(valueDiff),1);
for ii = 1:numel(valueDiff)
    diffThis = valueDiff(ii);
    
    % Get performances of decoders that discriminated this size difference.
    V1this = decoderV1(diffBetweenValues==diffThis);
    V4this = decoderV4(diffBetweenValues==diffThis);
    
    % Calculate mean decoder performance for this size difference.
    decoderV1mean(ii,1) = nanmean(V1this);
    decoderV4mean(ii,1) = nanmean(V4this);
    
    % Calculate standard error of the mean (SEM).
    decoderV1sem(ii,1) = nanstd(V1this) / sqrt(sum(~isnan(V1this)));
    decoderV4sem(ii,1) = nanstd(V4this) / sqrt(sum(~isnan(V4this)));
end
end

%% End