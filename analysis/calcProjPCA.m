function newvec = calcProjPCA(newdata,y,oldmean)
%calcProjPCA
%
% Usage:
%   newvec = calcProjPCA(newdata,y,oldmean)
%
% Description: 
%   Project neuronal responses onto the inputted principal components.
%
% Inputs:
%   newdata : (num neurons x num stimuli) neuronal responses to each stimulus
%   y       : (num neurons x numvecs) of numvecs first eigenvectors of the
%             correlation matrix
%   oldmean : (num neurons x 1) the mean response per neuron
%
% Output:
%   newvec : (numvecs x num stimuli) new vectors
%
% History:
%   10/14/21  amn  Wrote it.

%% Parse the inputs
parser = inputParser();
parser.addRequired('newdata',@(newdata)(ismatrix(newdata)));
parser.addRequired('y',@(y)(ismatrix(y)));
parser.addRequired('oldmean',@(oldmean)(isvector(oldmean)));
parser.parse(newdata,y,oldmean);

newdata = parser.Results.newdata;
y       = parser.Results.y;
oldmean = parser.Results.oldmean;

%% Project the neuronal responses
s=size(newdata);
oldmean=repmat(oldmean,1,s(2));
newvec=((newdata-oldmean)'*y)';
end