function [excitations, theOI, L, M, S] = computeResponse(display, fovDeg, psf, theMosaic, image)
% REPLACE BY ISETPipelineToolbox/ConeResponse (class)
% Compute cone mosaic response, given mosaic parameters and RGB image.
%
% Syntax: 
%   [excitations, theOI] = computeResponse(presentationDisplay, fovDeg, psf, theMosaic, image)
%
% Description: 
%   This funciton computes the mean cone excitations and the optical image
%   for a given mosaic and RGB image. 
%
% Inputs: 
%   display     - Display struct, see "displayCreate"
%   fovDeg      - Foveal visual degree
%   psf         - Point spread function (optical image structure), see "oiCreate"
%   theMosaic   - Cone mosaic object, see "coneMosaicHex"
%   image       - a RGB image
%
% Outputs:
%   excitations - Array of mean cone excitations
%   theOI       - Optical image (of the RGB image)
%   L           - L cone mean excitations
%   M           - M cone mean excitations
%   S           - S cone mean excitations


% Load picture
meanLuminanceCdPerM2 = 100;
realizedStimulusScene = sceneFromFile(image, 'rgb', ...
    meanLuminanceCdPerM2, display);

% Set the angular scene width
realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', fovDeg);

% Human optics
theOI = oiCompute(psf, realizedStimulusScene);

% Cone excitations
nTrialsNum = 2;
emPath = zeros(nTrialsNum, 1, 2);

% Compute mosaic excitation responses
% with/without the eye movement path
excitations = theMosaic.compute(theOI, 'emPath', emPath);

sizeExci = size(excitations);
meanExci = reshape(excitations(1, :, :), [sizeExci(2), sizeExci(3)]);

% LMS cone excitation
L = meanExci(theMosaic.pattern == 2);
M = meanExci(theMosaic.pattern == 3);
S = meanExci(theMosaic.pattern == 4);

end

