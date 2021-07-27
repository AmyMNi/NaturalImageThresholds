function runMaskPool(varargin)
%runMaskPool
%
% Usage:
%   runMaskPool()
%
% Description:
%   Create a mask pool from all of the images such that all of the masks 
%   are drawn from the same distribution. For each image, take the average
%   intensity (for each RGB channel) per image block. During the 
%   experiment, for each mask, each mask block will be randomly drawn from 
%   one of the quantized images in the mask pool. Thus, the masks will have 
%   the same basic luminance and color as the images.
%
% Optional parameter/value:
%   'experimentName' : (string) Name of experiment folder (default: 'Experiment100')
% 
% History:
%   07/27/21  amn  Wrote it.

%% Parse the input
parser = inputParser();
parser.addParameter('experimentName', 'Experiment100', @ischar);
parser.parse(varargin{:});

experimentName = parser.Results.experimentName;

%% Set paths to folder
%
% Specify project name.
projectName = 'NaturalImageThresholds';

% Get calibration file (set for the local experiment machine by the project 
% local hook file) to determine which image folder to use.
calFile = getpref(projectName,'CalDataFile');

% Set path to folder.
if strcmp(calFile,'NaturalImageThresholdsCal_Amy')
    pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBsAmy');
else
    pathToFolder = fullfile(getpref(projectName,'BaseDir'),experimentName,'ImageRGBs');
end

%% Get names of all scene files in the input folder
%
% List .mat files in the image folder.
fileInfo = dir([pathToFolder '/*.mat']);

% Get the number of images in the folder.
nImages = numel(fileInfo);

%% Calculate the number of pixels per block
%
% Define the number of blocks/image for the mask.
nBlocks = 16;

% Load one image to get image size.
file1 = fullfile(pathToFolder,fileInfo(1).name);
temp = load(file1,'RGBImage'); image1 = temp.RGBImage; clear temp;
imageSize = size(image1,1);

% The number of blocks must evenly divide the number of image pixels.
nPixels = imageSize;
if rem(nPixels,nBlocks)~=0
    error('params.nBlocks must evenly divide the number of image pixels.');
end

% Calculate the number of pixels per block.
blockPixels = nPixels/nBlocks;

%% For each image, take the average intensity (for each RGB channel) per block
%
% Set up a matrix (nBlocks x nBlocks x RGB channels x number of images) for the mask pool.
maskPool = nan(nBlocks,nBlocks,3,nImages);

% Analyze each image.
for iii = 1:nImages
    
    % Load the the image.
    fileToLoad = fullfile(pathToFolder,fileInfo(iii).name);
    temp = load(fileToLoad,'RGBImage'); image1 = temp.RGBImage; clear temp;
    
    % Flip image.
    image1 = image1(end:-1:1,:,:);
    
    % Take the average intensity (for each RGB channel) per block.
    for ii = 1:nBlocks
        for jj = 1:nBlocks
            for kk = 1:3
                theBlock = image1((ii-1)*blockPixels+1:ii*blockPixels,(jj-1)*blockPixels+1:jj*blockPixels,kk);
                blockRGB = mean(theBlock(:));
                maskPool(ii,jj,kk,iii) = blockRGB;
            end
        end
    end
end

%% Save the mask pool

fileToSave = fullfile(pathToFolder,'maskPool.mat');
save(fileToSave,'maskPool','nBlocks','blockPixels');
fprintf('maskPool saved as %s\n', fileToSave);

%% End