%% Combine electrophysiological data sets and run analyzeElectrophysData.m
%
% History:
%   10/18/21  amn  Wrote it.

%% TO DO ITEMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: stim numbers (PCs)

% TO DO: general decoder

% TO DO: electrode inclusion criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Record of electrophys data file names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: update record below
% Raptor data example: Ram/data/main/ivory/211013/ivory_map_211013-133559_dense.mat
% Save file as name example: 211013-133559_set1.mat
% Save in folder: /Users/amy/Desktop/Brainard/Natural Image Thresholds/Electrophys Data/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 211013-133559_set1.mat :sRGB image set, 2 iterations/stim

% 211014-153026_set2.mat : RGB image set, 5 iterations/stim

% 211015-155409_set2.mat : RGB image set, 4 iterations/stim
% 211015-160308_set3.mat : RGB image set, 4 iterations/stim

% 211016-150335_set3.mat : RGB image set, 4 iterations/stim
% 211016-152640_set3.mat : RGB image set, 4 iterations/stim
% 211016-161735_set3.mat : RGB image set, 3 iterations/stim

% 211017-155234_set3.mat : RGB image set, 8 iterations/stim
% 211017-164400_set3.mat : RGB image set, 3 iterations/stim

% 211018-145436_set3.mat : RGB image set, 1 iterations/stim
% 211018-153731_set3.mat : RGB image set, 1 iterations/stim
% 211018-155308_set3.mat : RGB image set, 8 iterations/stim

%% Select electrophysiological data sets to analyze

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: select data sets below.
dataNames = { ...
            '211013-133559_set1.mat', ...
            '211014-153026_set2.mat', ...
            '211015-155409_set2.mat', ...
            '211015-160308_set3.mat', ...
            '211016-150335_set3.mat', ...
            '211016-152640_set3.mat', ...
            '211016-161735_set3.mat', ...
            '211017-155234_set3.mat', ...
            '211017-164400_set3.mat', ...
            '211018-145436_set3.mat', ...
            '211018-153731_set3.mat', ...
            '211018-155308_set3.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set path to data folder.
dataFolder = fullfile('/Users','amy','Desktop','Brainard','Natural Image Thresholds','Electrophys Data');

%% Run combineElectrophysData.m to combine data sets selected above
data = combineElectrophysData(dataNames,dataFolder);

%% Run analyzeElectrophysData.m to analyze data
dataAnalysis = analyzeElectrophysData(data,'plotFigures',true);

%% Save data
%{
% Set file name.
fileName = 'TEST';

% Save data struct.
pathToDataFile = fullfile(dataFolder,sprintf('%sdata.mat',fileName));
save(pathToDataFile,'data');
fprintf('\ndata was saved in:\n%s\n', pathToDataFile);
    
% Save data analysis struct.
pathToDataAnalysisFile = fullfile(dataFolder,sprintf('%sdataAnalysis.mat',fileName));
save(pathToDataAnalysisFile,'dataAnalysis');
fprintf('\ndataAnalysis was saved in:\n%s\n', pathToDataAnalysisFile);
%}

%% End