function [estimates,model] = fitWeibull(x,y)
%fitWeibull
%
% Usage:
%   fitWeibull(x,y);
%
% Description:
%   Estimates alpha, beta, and gamma parameters of a Weibull function.
%   p(1): threshold
%   p(2): slope
%   p(3): max
%
% Inputs:
%   x : (vector) x-axis values (values of tested parameter)
%   y : (vector) y-axis values (performance per parameter value)
%
% History:
%   06/18/21  amn  Wrote it.

%% Set fit parameters
options    = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'FunValCheck','off');
startpoint = [median(x), 3, .95];

%% Fit Weibull function
model      = @myfun;
estimates  = fminsearch(model,startpoint,options);

    function [sse,FittedCurve,p] = myfun(p)
        FittedCurve = p(3)-exp(-power(x/p(1),p(2)));
        ErrorVector = FittedCurve-y;
        sse = sum(ErrorVector.^2);
    end

%% If estimates equal startpoints, the fit failed
if isequal(estimates,startpoint)
    estimates = nan(1,3);
end

end