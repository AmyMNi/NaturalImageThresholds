function [estimates,model] = fitWeibullMax(x,y)
%fitWeibull
%
% Usage:
%   fitWeibull(x,y);
%
% Description:
%   Estimates alpha and beta parameters of a Weibull function.
%   Gamma parameter (max) is contrained to 1 (generally, based on gamma
%   parameter fit >1 using fitWeibull.m).
%   p(1): threshold
%   p(2): slope
%
% Inputs:
%   x : (vector) x-axis values (values of tested parameter)
%   y : (vector) y-axis values (performance per parameter value)
%
% History:
%   06/22/21  amn  Wrote it.

%% Set fit parameters
options    = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'FunValCheck','off');
startpoint = [median(x), 3];

%% Fit Weibull function
model      = @myfun;
estimates  = fminsearch(model,startpoint,options);

    function [sse,FittedCurve,p] = myfun(p)
        FittedCurve = 1-exp(-power(x/p(1),p(2)));
        ErrorVector = FittedCurve-y;
        sse = sum(ErrorVector.^2);
    end

%% If estimates equal startpoints, the fit failed
if isequal(estimates,startpoint)
    estimates = nan(1,2);
end

end