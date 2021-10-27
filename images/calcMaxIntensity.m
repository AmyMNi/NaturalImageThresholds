function maxIntensity = calcMaxIntensity(scene,cal)
% calcMaxIntensity
%
% Usage:
%   maxIntensity = calcMaxIntensity(scene,cal)
%
% Description:
%   Calculate the maximum intensity within the inputted scene.
%   Adapted from t_renderISET3dHyperspectral:
%   https://github.com/isetbio/ISET3DProjects/tree/main/renderISET3dHyperspectral
%
% Inputs:
%   scene : (struct) ISET3d scene info
%   cal   : (struct) Calibration file for the electrophysiology experimental machine
%
% Output:
%   maxIntensity : (scalar) Maximum intensity within the inputted scene
%
% History:
%   06/05/21  dhb  Wrote it.
%   10/27/21  amn  Edits to t_renderISET3dHyperspectral.

%% Parse the inputs
parser = inputParser();
parser.addRequired('scene',@(x)(isstruct(x)));
parser.addRequired('cal',@(x)(isstruct(x)));
parser.parse(scene,cal);

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

%% Calculate the maximum intensity in the image
maxIntensity = max(rgbCalFormat(:));

%% End