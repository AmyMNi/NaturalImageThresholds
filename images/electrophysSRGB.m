function sRGBImage = electrophysSRGB(scene,varargin)
% electrophysSRGB
%
% Usage:
%   electrophysSRGB(scene,varargin)
%
% Description:
%   Read in a rendered recipe and convert to an sRGB image (not gamma 
%   corrected) for the Natural Image Thresholds electrophysiology experiment.
%   Adapted from t_renderISET3dHyperspectral:
%   https://github.com/isetbio/ISET3DProjects/tree/main/renderISET3dHyperspectral
%
% Inputs:
%   scene : (struct) ISET3d scene info
%
% Optional parameters/values:
%   'showSRGB' : (logical) Whether to display the sRGBImage (default: false)
%
% History:
%   06/05/21  dhb  Wrote it.
%   10/12/21  amn  Edits to t_renderISET3dHyperspectral.

%% Parse the inputs
parser = inputParser();
parser.addRequired('scene',@(x)(isstruct(x)));
parser.addParameter('showSRGB',false,@islogical);
parser.parse(scene,varargin{:});

showSRGB = parser.Results.showSRGB;

%% Get the hyperspectral image data out of the scene
%
% The image is rows x cols x nWls - the third dimension
% is the spectral radiance at each pixel in photons/sec-m2-sr-nm.
wls = sceneGet(scene,'wave');
S = MakeItS(wls);
radianceImageQuanta = sceneGet(scene,'photons');

%% Get standard LMS cone spectral sensitivities in quantal units
%
% This gets us the standard CIE fundamentals, for a 2 degree field.
% We could adjust for observer age, if we wanted.  32 years old is
% the standard default.
coneParams = DefaultConeParams('cie_asano');
coneParams.ageYears = 32;
coneParams.fieldSizeDegrees = 2;
[~,T_energy,T_quanta] = ComputeObserverFundamentals(coneParams,S);

%% Convert image to cal format and get LMS
%
% LMS coordinates in units of isomerizations/cone-sec (foveal cone
% geometric parameters used to estimate cone quantal capture).
%
% The multiplication by S(2) handles the wavelength spacing in 
% the matrix multiplication approximation of the integral over
% wavelength.  The convention in ISET code is that units of radiance
% are per nm.  Our calibration routines, below, use a convention of
% units in per wavelength band, so once we are entirely in PTB land
% we don't multiply by the delta wavelength factor.
[radianceQuantaCalFormat,nX,nY] = ImageToCalFormat(radianceImageQuanta);
LMSExcitationsCalFormat = T_quanta*radianceQuantaCalFormat*S(2);

%% Check on energy/quanta conversion
%
% Convert radiance to Watts/m2-sr-nm
radianceEnergyCalFormat = QuantaToEnergy(S,radianceQuantaCalFormat);
LMSExcitationsCalFormatChk = T_energy*radianceEnergyCalFormat*S(2);
if (max(abs(LMSExcitationsCalFormatChk(:) - LMSExcitationsCalFormat(:))) > 1e-12*max(abs(LMSExcitationsCalFormat(:))))
    error('Energy/quanta conversion glitch somewhere');
end

%% Render sRGB versions of the rendered image
%
% sRGB is sort of a generic monitor standard.  Having sRGB
% versions of images is useful for talks and papers, where 
% you don't really know what device will be used to show it
% and thus using a standard is as good a guess as anything.
%
% Notice that the image comes out looking achromatic here. 
% sRGB is closer to most modern LCD monitors than the monitor
% described by that old calibration file.
%
% sRGB is based on the CIE XYZ color matching functions, so
% first step is to go from Radiance to XYZ.  And first step
% to do that is get the color matching functions.  The magic
% 683 makes the units of Y cd/m2 when we specificy radiance
% as Watts/sr-m2-nm and take the wavelength delta properly
% into account in the summation over wavelength
load T_xyz1931
T_xyz = 683*SplineCmf(S_xyz1931,T_xyz1931,S);
XYZCalFormat = T_xyz*radianceEnergyCalFormat*S(2);

%% Convert XYZ to sRGB primary values
%
% Same general issues with scaling as above confront us here
sRGBPrimaryCalFormat = XYZToSRGBPrimary(XYZCalFormat);
maxPrimaryValue = max(sRGBPrimaryCalFormat(:));
headroomFactor = 1;
if (maxPrimaryValue > headroomFactor)
    fprintf('Warning: Maximum sRGB primary intensity of %0.2g exceeds desired max of %0.2g\n',maxPrimaryValue,headroomFactor);
end
sRGBPrimaryCalFormatScaled = headroomFactor*sRGBPrimaryCalFormat/maxPrimaryValue;

if (any(sRGBPrimaryCalFormatScaled(:) < 0))
    fprintf('Warning: some sRGB primary values in rendered image are negative.\n');
    fprintf('\tWorth looking into how many and by how much\n');
    fprintf('\tThis routine simply truncates such values to 0\n');
end
sRGBPrimaryCalFormatScaled(sRGBPrimaryCalFormatScaled < 0) = 0;

% Convert a calibration format image back to a real image.
sRGBImage = CalFormatToImage(sRGBPrimaryCalFormatScaled,nX,nY);
if showSRGB
    figure; imshow(sRGBImage);
    title('sRGB rendering');
end

%% End