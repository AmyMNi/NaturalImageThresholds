function x = calcTukey(x,varargin)
%calcTukey
%
% Usage:
%   x = calcTukey(x);
%
% Description:
%   Exclude outliers with the Tukey Method:
%   1) Convert input into a single vector if necessary.
%   2) Replace outliers with NaN.
%   3) Convert vector back into original matrix size if necessary.
%
% Input:
%   x : (matrix) values to check for outliers
%
% Optional parameter/value:
%   'multiplier' : (scalar) Inclusion multiplier for the interquartile range (default: 1.5)
%
% History:
%   07/08/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addRequired('x',@(x)(ismatrix(x)));
parser.addParameter('multiplier',1.5,@isscalar);
parser.parse(x,varargin{:});

multiplier = parser.Results.multiplier;

%% Replace outliers with NaN
x(isinf(x)) = nan;
sizex = size(x);
numelx = numel(x);
x = reshape(x,numelx,1);
xiqr = iqr(x);
xTukey = xiqr*multiplier;
xquar = quantile(x,3);
xq1 = xquar(1);
xq3 = xquar(3);
x(x>(xq3+xTukey) | x<(xq1-xTukey)) = nan; %outlier exclusion
x = reshape(x,sizex);
end