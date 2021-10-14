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
%   pcor      : (scalar)
%   bad       :
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

%% Calculate specific decoder performance
%
% First set up data.
s=sum(data,2);
bad=~isfinite(s);
numstim=length(s);
data(bad,:)=nan;
predicted=nan(numstim,1);







for i=1:numstim,
    if isfinite(group(i))
        this=1:numstim;
        this(bad)=nan;
        this(i)=nan;
        this=this(isfinite(this));
        class = nan;
        try
            class = classify(data(i,:),data(this,:),group(this));
        end
        predicted(i)=class; 
    end
end

pcor=sum(predicted==group)/(numstim-sum(bad));

if pcor<.5,
    pcor=1-pcor;
end

end