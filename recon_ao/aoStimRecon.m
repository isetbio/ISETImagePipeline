%% aoStimRecon
%
% Descriptoin:
%    Script to see how well we can measure unique yellow percepts under AO
%    conditions.

% History:
%   07/29/22  lz    Wrote this sometime in the past
%   07/29/22  dhb, chr  Starting to dig into this

%% Clear
clear; close all;

%% Point at directory with data files for this subproject
%
% This will allow us to load in project specific precomputed information.
aoReconDir = getpref('ISETImagePipeline','aoReconDir');

%% Setup / Simulation parameters
%
% Note these are only for book-keeping, if these are changed, we need to
% build a new render matrix using the code below
nPixels = 64;
fieldSizeMinutes = 30;
fieldSizeDegs = fieldSizeMinutes/60;
eccXDegs = 2.0;
eccYDegs = 0.0; 

%% Load in some important variables
% 
% The file monoDispRender.mat contains (at least)
%    monoDisplay: special monochromatic display with narrow primaries
%    theConeMosaic: Wrapper object that contains a pre-generated,
%       half-degree cone mosaic at eccX = 2.0 and eccY = 0.0
%    render: Linear transformation between input pixel and cone response,
%       for the pre-generated cone mosaic. Used as the likelihood function.
load(fullfile(aoReconDir,'monoDispRender.mat'));

% If we were clever, we'd add code here to make sure what we loaded matched
% variable values set above.

%% Code for building new render matrix
%
% Run this code segement if you would like to
% rebuild a new mosaic and render matrix
buildNew = false;
if buildNew
    theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3.0);
    theConeMosaic.Display = monoDisplay;
    render = theConeMosaic.forwardRender([nPixels nPixels 3], ...
        'validation', false);
    render = double(render);
end

%% Change the optics to diffraction limited optics with no LCA
% 
% This allows us to simulate the conditions we have in AO experiment,
% and reconstruct from those conditions
pupilDiameterMm = 7.0;
theConeMosaic.PSF = ConeResponse.psfDiffLmt(pupilDiameterMm);

%% Need new render if we want to reconstruct with respect to the AO stimulus
reconstructWrtAO = true;
if (reconstructWrtAO)
    render = theConeMosaic.forwardRender([nPixels nPixels 3], ...
        'validation', false);
    render = double(render);
end

%% Show cone mosaic
theConeMosaic.visualizeMosaic();

%% Generate an image stimulus
% stimulus in the size of retinal degree
% should not exceed 'fieldSizeDegs'
stimSizeDegs = 0.1;
size = stimSizeDegs / fieldSizeDegs;
idxLB = round(nPixels * (0.5 - size / 2));
idxUB = round(nPixels * (0.5 + size / 2));
idxRange = idxLB : idxUB;

% Image stimulus with a gray background and yellow color
testImage = ones(nPixels, nPixels, 3) * 0.20;
testImage(idxRange, idxRange, 1) = 0.8;
testImage(idxRange, idxRange, 2) = 0.65;

% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimScene, ~, linear] = sceneFromFile(testImage, 'rgb', ...
                meanLuminanceCdPerM2, theConeMosaic.Display);
visualizeScene(stimScene, 'displayRadianceMaps', false);

%% Compute cone response
coneExcitations = theConeMosaic.compute(testImage);
theConeMosaic.visualizeOI()

% Visualization of the cone response note that we are using
% 'activationRange', [0 max(coneExcitations)] to avoid confusions due to
% small stimulus
figureHandle = figure(); axesHandle = [];
theConeMosaic.Mosaic.visualize(...
        'figureHandle', figureHandle, ...
        'axesHandle', axesHandle, ...
        'activation', theConeMosaic.LastResponse, ...
        'activationRange', [0 max(coneExcitations)], ...
        'plotTitle',  'Cone Response');

%% Run reconstruction
% with special prior built for the display
prior = load(fullfile(aoReconDir,'sparsePrior.mat'));

% Construct onstruct image estimator
regPara = 0.001; stride = 2;
estimator = PoissonSparseEstimator(render, inv(prior.regBasis), ...
            prior.mu', regPara, stride, [nPixels nPixels 3]);

% Estimate image
%
% Note: could scale up/down the cone excitation vector if the
% reconstruction is going out of range
scaleFactor = 1;
reconImage = estimator.runEstimate(coneExcitations * scaleFactor, ...
    'maxIter', 500, 'display', 'iter', 'gpu', false);

%% Show reconstruction
meanLuminanceCdPerM2 = [];
[stimScene, ~, linear] = sceneFromFile(gammaCorrection(reconImage, monoDisplay), 'rgb', ...
                                        meanLuminanceCdPerM2, theConeMosaic.Display);
visualizeScene(stimScene, 'displayRadianceMaps', false);
