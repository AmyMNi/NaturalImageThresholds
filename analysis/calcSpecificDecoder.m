function [predicted,pcor,bad] = calcSpecificDecoder(data,group)
%calcSpecificDecoder
%
% Usage:
%   [~,pcor] = calcSpecificDecoder(data,group)
%
% Description:
%   Calculate the performance of a linear specific decoder, trained and
%   tasked (with leave-one-out cross validation) with discriminating the
%   neuronal population's responses to Stimulus #1 from the neuronal
%   population's responses to Stimulus #2.
%
% Inputs:
%   data  : (num stimuli x num neurons) neuronal responses to each stimulus
%   group : (num stimuli x 1) stimulus group (Stimulus #1 = 0, Stimulus #2  = 1)
%
% Outputs:
%   predicted : (num stimuli x 1) Stimulus group predicted by the decoder
%   pcor      : (scalar) Proportion of stimuli correctly predicted by the decoder
%   bad       : (num stimuli x 1) (logical) Stimuli not included due to nonfinite response values
%
% History:
%   10/14/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addRequired('data',@(data)(ismatrix(data)));
parser.addRequired('group',@(group)(isvector(group)));
parser.parse(data,group);

data  = parser.Results.data;
group = parser.Results.group;

%% Set up data.
s=sum(data,2);
bad=~isfinite(s);
numStim=length(s);
data(bad,:)=nan;
predicted=nan(numStim,1);

%% Perform dimension reduction on the neuronal population responses
%
% Because the Matlab function 'classify' used to determine decoder
% performance below requires that there are more training stimuli than
% dimensions of the neuronal population response, perform PCA on the
% neuronal population responses and analyze a number of dimensions that is
% less than the number of stimuli that will be used to train the decoder.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: this method maximizes the number of dimensions analyzed for this
%       decoder. May want to, instead, use the same number of dimensions
%       across all analzyed decoders.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of training stimuli that will be used.
numStimTrain = numStim - 1;

% Calculate the number of PCs (principal components) to analyze.
numPC = numStimTrain - 2;

% The number of PCs cannot exceed the number of neurons.
if numPC > size(data,2)
    numPC = size(data,2);
end

% Give warning if the number of dimensions is less than 2.
if numPC < 2
    fprintf(2,'Warning: number of population dimensions to analyze = %d\n',numPC);
end

% Perform dimensionality reduction on the neuronal population responses.
[y,~,~,psi] = pca2(data',numPC);

% Project responses onto the above PCs.
dataPC = calcProjPCA(data',y,psi)';

%% Per leave-one-out stimulus, calculate stimulus group predicted by the decoder.
for ii=1:numStim
    if isfinite(group(ii))
        this=1:numStim;
        this(bad)=nan;
        this(ii)=nan;
        this=this(isfinite(this));
        class = nan;
        try
            % Matlab function inputs: (sample,training,group)
            class = classify(dataPC(ii,:),dataPC(this,:),group(this));
        catch
            warning('classify function failed')
        end
        predicted(ii)=class; 
    end
end

%% Calculate proportion of stimuli correctly predicted by the decoder
correct = sum(predicted==group);
total = numStim - sum(bad) - sum(isnan(predicted));
pcor = correct/total;
if pcor==0
    pcor = nan;
end

end