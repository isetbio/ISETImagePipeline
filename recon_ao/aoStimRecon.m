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
% Spatial parameters.  Common to forward and recon models
nPixels = 100;
fieldSizeMinutes = 30;
fieldSizeDegs = fieldSizeMinutes/60;
eccXDegs = 2.0;
eccYDegs = 0.0;

% Display parameters.
%
% Display options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
%
% Also specify gamma table parameters
displayName = 'conventional';
displayGammaBits = 16;
displayGammaGamma = 2;
switch (displayName)
    case 'conventional'
        displayFieldName = 'CRT12BitDisplay';
        overwriteDisplayGamma = false;
    case 'mono'
        overwriteDisplayGamma = true;
        displayFieldName = 'monoDisplay';
    otherwise
        error('Unknown display specified');
end

% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional
%                           display.
sparsePriorName = 'conventionalSparsePrior.mat';

% Forward rendering parameters
%
% Use AO in forward rendering?
% This determins pupil diameter which typically differs in AO
forwardAORender = true;
if (forwardAORender)
    forwardPupilDiamMM = 7;
else
    forwardPupilDiamMM = 3;
end

% Residual defocus for forward rendering
forwardDefocusDiopters = 0;

% Force build and save
buildNewForward = false;

% Recon rendering parameters
useForwardRenderingForRecon = true;

%% Set forward render filenname
if (forwardAORender)
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%d.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters);
else
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%d.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters);
end

%% Build render matrix if desired/needed
%
% Run this code segement if you would like to rebuild a new mosaic and
% render matrix.  This also gets run if there is no cached file corresponding
% to the desired parameters. Once built, this file can be loaded from cache
% for quicker running.
if (buildNewForward || ~exist(fullfile(aoReconDir,forwardRenderStructureName),'file'))
    % Get display
    theDisplayLoad = load(fullfile(aoReconDir,[displayName 'Display.mat']));
    eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
    if (overwriteDisplayGamma)
        gammaInput = linspace(0,1,2^displayGammaBits)';
        gammaOutput = gammaInput.^displayGammaGamma;
        theDisplay.gamma = gammaOutput(:,[1 1 1]);
    end
    clear theDisplayLoad;

    % Create and setup cone mosaic
    %
    % For AO, we create a dummy object with 3 mm pupil and then adjust
    % pupil and make the OI diffraction limited with no LCA.  The field
    % name PSF is not optimal, because it is actually OI. We need the dummy
    % 3 mm pupil because the code that pulls out the initial Polens optics
    % checks that the desired pupil size is smaller than the 4 mm at which
    % those data were measured.
    if (forwardAORender)
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
    forwardRenderStructure.AORender = forwardAORender;
    forwardRenderStructure.defocusDiopters = forwardDefocusDiopters;
    save(fullfile(aoReconDir,forwardRenderStructureName),'forwardRenderStructure');

else
    % If not building, load cached file.  If we're doing this, it exists
    % as checked above. After laod, check that cached parameters match current parameters
    %
    % Read and check that loaded structure is as expected
    clear forwardRenderStructure;
    load(fullfile(aoReconDir,forwardRenderStructureName),'forwardRenderStructure');
    if (forwardRenderStructure.eccX ~= eccXDegs || forwardRenderStructure.eccY ~= eccYDegs)
        error('Precomputed forward rendering matrix not computed for current eccentricity');
    end
    if (forwardRenderStructure.fieldSizeDegs ~= fieldSizeDegs)
        error('Precomputed forward rendering matrix not computed for current field size');
    end
    if (forwardRenderStructure.nPixels ~= nPixels)
        error('Precomputed forward rendering nPixels not equal to current value');
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

% Set forward variables from loaded/built structure
forwardRenderMatrix = forwardRenderStructure.renderMatrix;
forwardConeMosaic = forwardRenderStructure.theConeMosaic;
forwardOI = forwardConeMosaic.PSF;

% Set display variable
theDisplay = forwardRenderStructure.theDisplay;

% Clear forward render structure
clear forwardRenderStructure;

%% Reconstruction rendering parameters
%
% These are how the reconstruction algorithm thinks the image
% was formed.  Need not be the same as the forward rendering.
if (useForwardRenderingForRecon)
    reconRenderMatrix = forwardRenderMatrix;
    reconConeMosaic = forwardConeMosaic;
    reconOI = forwardOI;
    reconPupilDiamMM = forwardPupilDiamMM;
    reconAORender = forwardAORender;
    reconDefocusDiopters = forwardDefocusDiopters;
else
    error('Need to implement separate recon rendering setup');
end

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

%% Show forward cone mosaic
forwardConeMosaic.visualizeMosaic();

%% Generate an image stimulus
% stimulus in the size of retinal degree
% should not exceed 'fieldSizeDegs'
stimSizeDegs = 0.4;
stimSizeFraction = stimSizeDegs / fieldSizeDegs;
idxLB = round(nPixels * (0.5 - stimSizeFraction / 2));
idxUB = round(nPixels * (0.5 + stimSizeFraction / 2));
idxRange = idxLB:idxUB;

% Image stimulus with a gray background and yellow color
stimulusImageRGB = ones(nPixels, nPixels, 3) * 0.20;
stimulusImageRGB(idxRange, idxRange, 1) = 0.8;
stimulusImageRGB(idxRange, idxRange, 2) = 0.65;

% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, testImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false);

%% Compute forward retinal image and excitations.
%
% We'll reconstruct from these.
forwardOI = oiCompute(stimulusScene,forwardOI);
visualizeOpticalImage(forwardOI);
forwardExcitationsToStimulus = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));

% Check forward exciations calculation another way.  Also shows another way
% to visualize the retinal image, but this only is done when the check
% fails.
forwardExcitationsToStimulusCheck = forwardConeMosaic.compute(stimulusImageRGB);
if (max(abs(forwardExcitationsToStimulusCheck-forwardExcitationsToStimulus)) ~= 0)
    forwardConeMosaic.visualizeOI()
    error('Two ways of doing the same thing do not agree');
end
temp = squeeze(forwardConeMosaic.LastResponse);
if (max(abs(temp-forwardExcitationsToStimulus)) ~= 0)
    error('Last excitations in object not as we expect');
end

% This code lets us figure out which wavelengths had non-zero photon counts
%
% for ii = 1:size(forwardOI.data.photons,3)
%     temp = forwardOI.data.photons(:,:,ii);
%     if (any(max(temp(:)) >= 1e12))
%         fprintf('Plane %d has non-zero photons, max %g, mean %g, center pixel value %g\n',ii,max(temp(:)),mean(temp(:)),temp(round(nPixels/2),round(nPixels/2)));
%     end
% end

% Visualization of the cone response note that we are using
% 'activationRange', [0 max(coneExcitations)] to avoid confusions due to
% small stimulus
figureHandle = figure(); axesHandle = [];
forwardConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', axesHandle, ...
    'activation', forwardConeMosaic.LastResponse, ...
    'activationRange', [0 max(forwardExcitationsToStimulus)], ...
    'plotTitle',  'Cone Response');

%% Run reconstruction
% with special prior built for the display
prior = load(fullfile(aoReconDir,sparsePriorName));

% Construct onstruct image estimator
regPara = 0.001; stride = 2;
estimator = PoissonSparseEstimator(reconRenderMatrix, inv(prior.regBasis), ...
    prior.mu', regPara, stride, [nPixels nPixels 3]);

% Estimate
%
% Scale excitations to take into account difference in forward and recon
% pupil sizes.  This helps keep things in range.
scaleFactor = (forwardPupilDiamMM/reconPupilDiamMM)^2;
reconImage = estimator.runEstimate(forwardExcitationsToStimulus * scaleFactor, ...
    'maxIter', 500, 'display', 'iter', 'gpu', false);

% Show reconstruction
meanLuminanceCdPerM2 = [];
[reconScene, ~, reconImageLinear] = sceneFromFile(gammaCorrection(reconImage, theDisplay), 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
reconScene = sceneSet(reconScene, 'fov', fieldSizeDegs);
visualizeScene(reconScene, 'displayRadianceMaps', false);

% Compute forward excitations from reconstruction
% And compare with stimulus exciations
forwardOI = oiCompute(reconScene,forwardOI);
forwardExcitationsToRecon = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));
figure; clf; hold on;
plot(forwardExcitationsToStimulus,forwardExcitationsToRecon,'ro','MarkerFaceColor','r','MarkerSize',10);
axis('square');

