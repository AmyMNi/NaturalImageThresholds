function newvec = calcProjPCA(newdata,y,oldmean)
%calcProjPCA
%
% Usage:
%   newvec = calcProjPCA(newdata,y,oldmean)
%
% Description:
%   Project neuronal responses onto the input eigenvectors.
%
% Inputs:
%   newdata : (num neurons x num stimuli) neuronal responses to each stimulus
%   y       : (num neurons x numvecs) of numvecs first eigenvectors
%   oldmean : (num neurons x 1) the mean response per neuron
%
% Output:
%   newvec : (numvecs x num stimuli) new vectors
%
% History:
%   10/14/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addRequired('newdata',@(x)(ismatrix(x)));
parser.addRequired('y',@(x)(ismatrix(x)));
parser.addRequired('oldmean',@(x)(isvector(x)));
parser.parse(newdata,y,oldmean);

newdata = parser.Results.newdata;
y       = parser.Results.y;
oldmean = parser.Results.oldmean;

%% Project the neuronal responses
s=size(newdata);
oldmean=repmat(oldmean,1,s(2));
newvec=((newdata-oldmean)'*y)';
end