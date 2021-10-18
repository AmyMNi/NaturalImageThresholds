function dataAnalysis = analyzeElectrophysData(data,varargin)
%analyzeElectrophysData
%
% Usage:
%   dataAnalysis = analyzeElectrophysData(data)
%
% Description:
%   Analyze electrophysiological data from the inputted data struct.
%   Save the results in the output struct, in the specified output folder.
%
% Input:
%   'data' : (struct) Struct of data to be analyzed, created by combineElectrophysData.m
%
% Optional parameters/values:
%   'plotFigures' : (logical) Plot figures if option is on (default: true)
%   'saveData'    : (logical) Save data if option is on (default: true)
%
% Output:
%   'dataAnalysis' : (struct) data struct with additional analysis results
%
% History:
%   10/13/21  amn  Wrote it.
%   10/18/21  amn  Updated to work with runAnalyzeElectrophysData.m and
%                  combineElectrophysData.m to analyze multiple data sets

%% Parse the inputs
parser = inputParser();
parser.addRequired('data',@(x)(isstruct(x)));
parser.addParameter('plotFigures', true, @islogical);
parser.addParameter('saveData', true, @islogical);
parser.parse(data,varargin{:});

data        = parser.Results.data;
plotFigures = parser.Results.plotFigures;
saveData    = parser.Results.saveData;

%% Set up new dataAnalysis struct
dataAnalysis = data;

%% Get output file path
pathToOutput = dataAnalysis.pathToOutput;

%% Get data variables

positions     = dataAnalysis.positions;
rotations     = dataAnalysis.rotations;
depths        = dataAnalysis.depths;
imagePosition = dataAnalysis.imagePosition;
imageRotation = dataAnalysis.imageRotation;
imageDepth    = dataAnalysis.imageDepth;
V1respInc     = dataAnalysis.V1respInc;
V4respInc     = dataAnalysis.V4respInc;

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
    fprintf('\ndataAnalysis was saved in:\n%s\n', pathToOutput);
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