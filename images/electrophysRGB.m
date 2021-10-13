function RGBImage = electrophysRGB(scene,cal,varargin)
% electrophysRGB
%
% Usage:
%   electrophysRGB(scene,cal)
%
% Description:
%   Read in a rendered recipe, extract the hyperspectral image data, 
%   compute LMS cone excitations, and use PTB routines to go from there to
%   a (nongamma-corrected) metameric RGB image for the Natural Image
%   Thresholds electrophysiology experiment.
%   Adapted from t_renderISET3dHyperspectral:
%   https://github.com/isetbio/ISET3DProjects/tree/main/renderISET3dHyperspectral
%
% Inputs:
%   scene : (struct) ISET3d scene info
%   cal   : (struct) Calibration file for the electrophysiology experimental machine
%
% Optional parameters/values:
%   'showRGB'  : (logical) Whether to display the RGBImage  (default: false)
%
% History:
%   06/05/21  dhb  Wrote it.
%   10/13/21  amn  Edits to t_renderISET3dHyperspectral.

%% Parse the inputs
parser = inputParser();
parser.addRequired('scene',@(x)(isstruct(x)));
parser.addRequired('cal',@(x)(isstruct(x)));
parser.addParameter('showRGB',false,@islogical);
parser.parse(scene,cal,varargin{:});

showRGB  = parser.Results.showRGB;

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

%% Initialize the sensor color space for use in calibration
cal = SetSensorColorSpace(cal,T_energy,S);

%% Go from LMS to device primary space
rgbCalFormat = SensorToPrimary(cal,LMSExcitationsCalFormat);

%% Scale into gamut
%
% Nothing in our rendering pipeline guarantees that the maximum
% intensity of the image is within the gamut of the monitor.  We
% could address this by scaling the illumination intensity of
% the light source to bring the maximum primary RGB value down
% lower than 1, or we can scale at this stage. 
%
% When we do the experiment, we have to be careful to scale all
% of the images the same way, so it may be cleaner to scale the
% intensity of the light source in the rendering, and then throw
% an error at this stage if the intensity is out of gamut.
headroomFactor = 1;
maxPrimaryValue = max(rgbCalFormat(:));
if (maxPrimaryValue > headroomFactor)
    fprintf('Warning: Maximum primary intensity of %0.2g exceeds desired max of %0.2g\n',maxPrimaryValue,headroomFactor);
end
rgbCalFormatScaled = headroomFactor*rgbCalFormat/maxPrimaryValue;

% It's also possible to get negative rgb values.  This happens
% if the saturation of one of the rendered pixels exceeds
% what the monitor gamut can display.  Not much to do here 
% except truncate to positive, and perhaps let the user know.
if (any(rgbCalFormatScaled(:) < 0))
    fprintf('Warning: some primary values in rendered rgb image are negative.\n');
    fprintf('\tWorth looking into how many and by how much\n');
    fprintf('\tThis routine simply truncates such values to 0\n');
end
rgbCalFormatScaled(rgbCalFormatScaled < 0) = 0;

%% Convert back to an image
RGBImage = CalFormatToImage(rgbCalFormatScaled,nX,nY);

%% Display
if showRGB
    figure; imshow(RGBImage);
    title('Calibrated RGB rendering');
end

%% End