function dataAnalysis = analyzeElectrophysData(varargin)
%analyzeElectrophysData
%
% Usage:
%   data = analyzeElectrophysData()
%
% Description:
%   Analyze electrophysiological data from a single data collection session.
%   Save the results (struct 'data') in the specified output folder.
%
% Optional parameters/values:
%   'plotFigures'    : (logical) Plot figures if option is on (default: true)
%   'saveData'       : (logical) Save data if option is on (default: true)
%
% History:
%   10/13/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addParameter('plotFigures', true, @islogical);
parser.addParameter('saveData', true, @islogical);
parser.parse(varargin{:});

plotFigures    = parser.Results.plotFigures;
saveData       = parser.Results.saveData;

%% Set paths to input and output files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set path to data file.
pathToInput  = fullfile('/Users','amy','Desktop','Brainard','Natural Image Thresholds','Electrophys Data','XXXXX.mat');

% Set path to output file.
pathToOutput = fullfile('/Users','amy','Desktop','Brainard','Natural Image Thresholds','Electrophys Data','XXXXXAnalysis.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Load data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: update with data file variable

% Load specified data file.
temp = load(pathToInput,'data'); data = temp.data; clear temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save data analysis results

if saveData
    % Save dataAnalysis struct.
    save(pathToOutput,'dataAnalysis');
    fprintf('\nData was saved in:\n%s\n', pathToOutput);
end
end
%% End