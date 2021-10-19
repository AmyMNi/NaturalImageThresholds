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

%% Electrophys data & image info organization

% Electrophys data variables:
 
% params : (struct 1 x num stimuli)
%       set: for the Natural Image Thresholds experiment, set = 7
%       num: number of the stimulus (to match to imgInfo table with image info)
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
%       column 1: imgNumber (image number to match to params.num)
%       column 2: bananaPosition (central object position in iset3d mm compared to scene center)
%       column 3: backgroundRotation (background objects rotation in degrees compared to original scene)
%       column 4: backgroundDepth (background objects depth in mm compared to original scene)

%% Create struct to hold data
data = struct;

%% Combine data from each data file
V1resp      = [];
V1resp_base = [];
V4resp      = [];
V4resp_base = [];
imagePosition = [];
imageRotation = [];
imageDepth    = [];

% Get data from one data set at a time.
for dd = 1:numel(dataNames)

    % Set path to data file.
    nameThis = dataNames{dd};
    pathToDataFile = fullfile(dataFolder,nameThis);
    
    % Load specified data file.
    clear params;    temp = load(pathToDataFile,'params');    params    = temp.params;    clear temp;
    clear eid;       temp = load(pathToDataFile,'eid');       eid       = temp.eid;       clear temp;
    clear resp;      temp = load(pathToDataFile,'resp');      resp      = temp.resp;      clear temp;
    clear resp_base; temp = load(pathToDataFile,'resp_base'); resp_base = temp.resp_base; clear temp;
    
    % Get image set name.
    p_   = strfind(nameThis,'_');
    pend = strfind(nameThis,'.');
    imageSet = str2double(regexp(nameThis(p_:pend),'[+-]?\d*','Match'));
    
    % Set path to image info file.
    pathToDataInfo = fullfile(dataFolder,sprintf('imgInfo_set%d.mat',imageSet));
    
    % Load specified image info file.
    clear imgInfo; temp = load(pathToDataInfo,'imgInfo'); imgInfo = temp.imgInfo; clear temp;
    
    % Check that all of the stimuli were from the Natural Image Thresholds experiment.
    if any([params.set] ~= 7)
        fprintf(2,'Warning: not all stimuli from %s were from the Natural Image Thresholds project\n',nameThis);
    end
    
    % Check that the eid matrix from the last data set matches this one.
    if exist('eidLast','var')==1
        fprintf(2,'Warning: the electrode order from %s does''t match the electrode order from the prior data set\n',nameThis);
    end

    % Exclude stimuli presented during incomplete trials (fixation broken).
    excludeStim = ~[params.good];
    params   (excludeStim)   = [];
    resp     (excludeStim,:) = [];
    resp_base(excludeStim,:) = [];
    
    % Get electrode numbers for the V1/V2 array (V1/V2: array number 2).
    electrodeV1 = eid(eid(:,1)==2,2);
    % Get electrode numbers for the V4 array (V4: array number 1).
    electrodeV4 = eid(eid(:,1)==1,2);
     
    % Get responses for each brain area.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: based on the same time period for both V1/V2 and V4 (50-350 ms).
    % NOTE: for different time periods, see params.spikes described above.
    % Get responses for V1/V2.
    V1respThis      = resp     (:,ismember(eid(:,2),electrodeV1));
    V1resp_baseThis = resp_base(:,ismember(eid(:,2),electrodeV1));
    % Get responses for V4.
    V4respThis      = resp     (:,ismember(eid(:,2),electrodeV4));
    V4resp_baseThis = resp_base(:,ismember(eid(:,2),electrodeV4));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine included electrodes for this data set, based on the following criteria:
    %   stimulus-evoked firing rate statsig greater than baseline firing rate, and
    %   minimum stimulus-evoked firing rate.
    % The electrodes included in the combined data will be determined below.
    frMin = 10; %Hz
    alphaLevel = .05;
    % Calculate included electrodes: V1.
    electrodeV1includedThis = false(1,size(V1respThis,2));
    for ii = 1:size(V1respThis,2)
        % Two-sample one-tailed t-test.
        h = ttest2(V1respThis(:,ii),V1resp_baseThis(:,ii),'Tail','right','Alpha',alphaLevel);
        if h==1 && nanmean(V1respThis(:,ii)) > frMin
            electrodeV1includedThis(ii) = true;
        end
    end
    % Calculate included electrodes: V4.
    electrodeV4includedThis = false(1,size(V4respThis,2));
    for ii = 1:size(V4respThis,2)
        % Two-sample one-tailed t-test.
        h = ttest2(V4respThis(:,ii),V4resp_baseThis(:,ii),'Tail','right','Alpha',alphaLevel);
        if h==1 && nanmean(V4respThis(:,ii)) > frMin
            electrodeV4includedThis(ii) = true;
        end
    end
    
    % Determine electrodes to includes across all combined data so far.
    if exist('electrodeV1included','var')==1 && exist('electrodeV4included','var')==1
        electrodeV1included = electrodeV1included & electrodeV1includedThis;
        electrodeV4included = electrodeV4included & electrodeV4includedThis;
    else
        electrodeV1included = electrodeV1includedThis;
        electrodeV4included = electrodeV4includedThis;
    end
    
    % Convert imgInfo table to matrix.
    % columns:  imgNumber  bananaPosition  backgroundRotation  backgroundDepth
    imgInfo = table2array(imgInfo);
    
    % Per stimulus, get image number to match to imgInfo table with image info.
    imageNum = [params.num]';
    [~,indImageNum] = ismember(imageNum, imgInfo(:,1));
    
    % Per stimulus, get central object (banana) position.
    imagePositionThis = imgInfo(indImageNum,2);
    
    % Per stimulus, get background object (branches & leaves) rotation.
    imageRotationThis = imgInfo(indImageNum,3);
    
    % Per stimulus, get background object (branches & leaves) depth.
    imageDepthThis = imgInfo(indImageNum,4);
    
    % Save eid matrix to compare to the next data set's eid matrix.
    eidLast = eid;
    
    % Add this data set to the combined data.
    V1resp      = [V1resp; V1respThis];
    V1resp_base = [V1resp_base; V1resp_baseThis];
    V4resp      = [V4resp; V4respThis];
    V4resp_base = [V4resp_base; V4resp_baseThis];
    imagePosition = [imagePosition; imagePositionThis];
    imageRotation = [imageRotation; imageRotationThis];
    imageDepth    = [imageDepth; imageDepthThis];
end

%% Get responses for electrodes included across combined data
% These are the electrodes that were included in each of the data sets.
V1respInc      = V1resp(:,electrodeV1included);
V1resp_baseInc = V1resp_base(:,electrodeV1included);
V4respInc      = V4resp(:,electrodeV4included);
V4resp_baseInc = V4resp_base(:,electrodeV4included);

%% Set up data struct output
%
% Save the image info per stimulus.
data.imagePosition = imagePosition;
data.imageRotation = imageRotation;
data.imageDepth = imageDepth;

% Save the electrode number order for all of the electrodes.
data.electrodeV1 = electrodeV1;
data.electrodeV4 = electrodeV4;

% Save the neuronal response matrices with all electrodes included.
data.V1resp = V1resp;
data.V1resp_base = V1resp_base;
data.V4resp = V4resp;
data.V4resp_base = V4resp_base;

% Save the electrode number order for the included electrodes.
data.electrodeV1included = electrodeV1(electrodeV1included);
data.electrodeV4included = electrodeV4(electrodeV4included);

% Save the neuronal response matrices with only the included electrodes.
data.V1respInc = V1respInc;
data.V1resp_baseInc = V1resp_baseInc;
data.V4respInc = V4respInc;
data.V4resp_baseInc = V4resp_baseInc;

end
%% End