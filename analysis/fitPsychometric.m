function [xx,FittedCurve,threshold,pse] = fitPsychometric(x,NumPos,OutOfNum)
%fitPsychometric
%
% Usage:
%   [xx,FittedCurve,threshold] = fitPsychometric(x,NumPos,OutOfNum);
%
% Inputs:
%   x        : (vector) stimulus levels tested
%   NumPos   : (vector) for each stimulus level tested, the number of trials a positive response was given
%   OutOfNum : (vector) for each stimulus level tested, the total number of trials
%
% History:
%   07/05/21  amn  Wrote it.

%% Calculate x-axis values to plot
%
% Increase the number of stimulus level values to plot for a smooth fit.
xx = linspace(min(x),max(x),1000);

%% Calculate psychometric function fit using Palamedes Toolbox

% Psychometric function form (alternative: PAL_Weibull).
PF = @PAL_CumulativeNormal;

% 'paramsFree' is a boolean vector that determines what parameters get
% searched over (1: free parameter, 0: fixed parameter).
paramsFree = [1 1 1 1];  

% Set up starting points:
%   1st (mean of the cumulative normal): set to mean of stimulus levels tested
%   2nd (standard deviation): set to a fraction of the range of stimulus levels tested
searchGrid = [mean(x) 1/(max(x)-min(x)) 0 0];

% Set up lapse limits.
lapseLimits = [0 0.05];

% Set up standard options for Palamedes search.
options = PAL_minimize('options');

% Fit with Palemedes Toolbox.
[paramsValues] = PAL_PFML_Fit(x,NumPos,OutOfNum,searchGrid,paramsFree,PF, ...
    'lapseLimits',lapseLimits,'guessLimits',[],'searchOptions',options,'gammaEQlambda',true);

% Get fitted curve values.
FittedCurve = PF(paramsValues,xx);

% Calculate threshold: the difference between the stimulus levels at
% performances equal to 0.7602 and 0.5.
pse = PF(paramsValues,0.5,'inverse');
threshold = PF(paramsValues,0.7602,'inverse')-pse;