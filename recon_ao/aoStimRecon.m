%% Setup / Simulation parameters
% Note these are only for book-keeping,
% if these are changed, we need to build
% a new render matrix using the code below
nPixels = 64;
fieldSizeMinutes = 30;
fieldSizeDegs = fieldSizeMinutes/60;

eccX = 2.0; % In degs
eccY = 0.0; % In degs

%% Point at directory with data files for this subproject
aoReconDir = getpref('ISETImagePipeline','aoReconDir');

%% Load in some important variables
% monoDisplay: special monochromatic display with narrow primaries
% theConeMosaic: Wrapper object that contains a pre-generated,
% half-degree cone mosaic at eccX = 2.0 and eccY = 0.0
% render: Linear transformation between input pixel and cone response,
% for the pre-generated cone mosaic. Used as the likelihood function.
load(fullfile(aoReconDir,'monoDispRender.mat'));

%% Code for building new render matrix
% Run this code segement if you would like to
% rebuild a new mosaic and render matrix
buildNew = false;

if buildNew
theConeMosaic = ConeResponseCmosaic(eccX, eccY, ...
    'fovealDegree', fieldSizeDegs, 'pupilSize', 3.0);

theConeMosaic.Display = monoDisplay;
render = theConeMosaic.forwardRender([nPixels nPixels 3], ...
                                        validation, 'false')
render = double(render);
end

%% Show cone mosaic
theConeMosaic.visualizeMosaic();

%% Change the optics to diffraction limited optics with no LCA
% to simulate the condition we have in AO experiment
pupilDiameterMm = 7.0;
theConeMosaic.PSF = ConeResponse.psfDiffLmt(pupilDiameterMm);

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

% Visualization of the cone response
% note that we are using
% 'activationRange', [0 max(coneExcitations)]
% to avoid confusions due to small stimulus
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

% Run image reconstruction
% construct image estimator
regPara = 0.001; stride = 2;
estimator = PoissonSparseEstimator(render, inv(prior.regBasis), ...
            prior.mu', regPara, stride, [nPixels nPixels 3]);

% Run image reconstruction to estimate images
% Note: could scale up/down the cone excitation vector
% if the reconstruction is going out of range
reconImage = estimator.runEstimate(coneExcitations * 0.5, ...
    'maxIter', 500, 'display', 'iter', 'gpu', false);

%% Show reconstruction
meanLuminanceCdPerM2 = [];
[stimScene, ~, linear] = sceneFromFile(gammaCorrection(reconImage, monoDisplay), 'rgb', ...
                                        meanLuminanceCdPerM2, theConeMosaic.Display);
visualizeScene(stimScene, 'displayRadianceMaps', false);
