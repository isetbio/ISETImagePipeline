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
forwardDefocusDiopters = 0;
displayGammaBits = 16;
displayGammaGamma = 2;
AOForwardRender = false;
displayName = 'mono';
if (AOForwardRender)
    forwardPupilDiamMM = 7;
else
    forwardPupilDiamMM = 3;
end

% Set forward render name
if (AOForwardRender)
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%d.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters);
else
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%d.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters);
end

% Build and save control
buildNewForward = true;

%% Load in some important variables
% 
% The file monoDisplayRender_FieldSize_EccX_EccY_Pix_PupilDiam_AOFlag_DefocusDiopters.mat contains
%    theDisplay: Special monochromatic display with narrow-band primaries
%    theConeMosaic: Wrapper object that contains a pre-generated,
%       half-degree cone mosaic at eccX = 2.0 and eccY = 0.0
%    renderMatrix: Linear transformation between input pixel and cone response,
%       for the pre-generated cone mosaic. Used as the likelihood function.
% The name has numbers filed in for the various parameters, to help keep
% versions for different parameters separate.

% Load cached forward file if it exists.  Force build if not.
% Check that cached parameters match current parameters
if (~buildNewForward)
    if (exist(fullfile(aoReconDir,forwardRenderStructureName),'file'))
        % Read and check that loaded structure is as expected
        forwardRenderStructure = load(fullfile(aoReconDir,forwardRenderStructureName));
        if (forwardRenderStructure.eccX ~= eccXDegs || forwardRenderStructure.eccY ~= eccY)
            error('Precomputed forward rendering matrix not computed for current eccentricity');
        end
        if (forwardRenderStructure.fieldSizeDegs ~= fieldSizeDegs)
            error('Precomputed forward rendering matrix not computed for current field size');
        end
        if (forwardRenderStructure.nPixels ~= nPixels)
            error('Precomputed forward rendering nPixels not equal to current value');
        end
        if (any(forwardRenderStructure.imageSize ~= imageSize))
            error('Precomputed forward rendering imageSize not equal to current value')
        end
        if (forwardRenderStructure.pupilDiamMM ~= forwardPupilDiamMM)
            error('Precompued forward pupil size not equal to current value');
        end
        if (forwardRenderStructure.AORender ~= forwardAORender)
            error('Precompued forward AO state not equal to current value');
        end
        if (forwardRenderStructure.defocusDiopters ~= forwardDefocusDiopters)
            error('Precompued forward defocus diopters not equal to current value');
        end
    end
else
    % Force build when desired file does not exist
    buildNewForward = true;
end

%% Build new render matrix if desired/needed
%
% Run this code segement if you would like to rebuild a new mosaic and
% render matrix.  It also gets run if there is no cached file corresponding
% to the desired parameters.
if (buildNewForward)
    % Get display
    theDisplayLoad = load(fullfile(aoReconDir,[displayName 'Display.mat']));
    eval(['theDisplay = theDisplayLoad.' displayName 'Display;']);
    gammaInput = linspace(0,1,2^displayGammaBits);
    gammaOutput = gammaInput.^displayGammaGamma;
    theDisplgammaTable = gammaOutput(:,[1 1 1]);
    clear theDisplayLoad;

    % Create and setup cone mosaic
    %
    % For AO, we create a dummy object with 3 mm pupil and then adjust
    % pupil and make the OI diffraction limited with no LCA.  The field
    % name PSF is not optimal, because it is actually OI. We need the dummy
    % 3 mm pupil because the code that pulls out the initial Polens optics
    % checks that the desired pupil size is smaller than the 4 mm at which
    % those data were measured.
    if (AOForwardRender)
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', false); 
        theConeMosaic.PSF = ConeResponse.psfDiffLmt(forwardPupilDiamMM);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', forwardPupilDiamMM, 'useRandomSeed', false); 
    end
    theConeMosaic.Display = theDisplay;

    % Generate render matrix
    forwardRenderMatrix = theConeMosaic.forwardRender([nPixels nPixels 3], ...
        true, true, 'useDoublePrecision', true);
    forwardRenderMatrix = double(forwardRenderMatrix);

    % Push new info back into structure and save
    forwardRenderStructure.theDisplay = theDisplay;
    forwardRenderStructure.renderMatrix = forwardRenderMatrix;
    forwardRenderStructure.theConeMosaic = theConeMosaic;
    forwardRenderStructure.fieldSizeDegs = fieldSizeDegs;
    forwardRenderStructure.eccX = eccXDegs;
    forwardRenderStructure.eccY = eccYDegs;
    forwardRenderStructure.nPixels = nPixels;
    forwardRenderStructure.pupilDiamMM = forwardPupilDiamMM;
    forwardRenderStructure.AORender = AOForwardRender;
    forwardRenderStructure.defocusDiopters = forwardDefocusDiopters;
    save(fullfile(aoReconDir,forwardRenderStructureName));
end

% Set forward variables from loaded/built structure
forwardRenderMatrix = forwardRenderStructure.renderMatrix;
forwardConeMosaic = forwardRenderStructure.theConeMosiac;
forwardOI = forwardConeMosaic.PSF;
theDisplay = forwardRenderStructure.theDisplay;
clear forwardRenderStructure;

% Get the optical image structures
OIAO = ConeResponse.psfDiffLmt;

%% Need new render if we want to reconstruct with respect to the AO stimulus
reconstructWrtAO = true;
buildNewAO = false;
if (reconstructWrtAO)
    if (buildNewAO)
        forwardRenderMatrix = theConeMosaic.forwardRender([nPixels nPixels 3], ...
            'validation', false);
        forwardRenderMatrix = double(forwardRenderMatrix);
    else
        load(fullfile(aoReconDir,'aoOpticsRenderMatrix.mat'));
    end
end

%% Show cone mosaic
theConeMosaic.visualizeMosaic();

%% Generate an image stimulus
% stimulus in the size of retinal degree
% should not exceed 'fieldSizeDegs'
stimSizeDegs = 0.4;
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
stimScene = sceneSet(stimScene, 'fov', stimSizeDegs);
visualizeScene(stimScene, 'displayRadianceMaps', false);

% Compute and visual retinal images
forwardOI = oiCompute(stimScene,forwardOI);
visualizeOpticalImage(forwardOI);
theResponseRegular = theConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered');
coneExcitationsCheckRegular = theResponseRegular(:);

OIAO = oiCompute(stimScene,OIAO);
visualizeOpticalImage(OIAO);
theResponseAO = theConeMosaic.Mosaic.compute(OIAO, 'opticalImagePositionDegs', 'mosaic-centered');
coneExcitationsCheckAO = theResponseAO(:);

for ii = 1:101
    temp = OIAO.data.photons(:,:,ii);
    if (any(max(temp(:)) >= 1e12))
        fprintf('Plane %d has non-zero photons, max %g, mean %g, center %g\n',ii,max(temp(:)),mean(temp(:)),temp(40,40));
    end
end

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
estimator = PoissonSparseEstimator(forwardRenderMatrix, inv(prior.regBasis), ...
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
[stimScene, ~, linear] = sceneFromFile(gammaCorrection(reconImage, theDisplay), 'rgb', ...
                                        meanLuminanceCdPerM2, theConeMosaic.Display);
visualizeScene(stimScene, 'displayRadianceMaps', false);
