%% aoStimRecon
%
% Descriptoin:
%    Script to see how well we can measure unique yellow percepts under AO
%    conditions.

% History:
%   07/29/22  lz    Wrote this sometime in the past
%   07/29/22  dhb, chr  Starting to dig into this
%   08/13/22  dhb   Lots of bookkeeping, cleaning, etc.

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
displayName = 'mono';
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

% Stimulus parameters
stimSizeDegs = 0.4;
stimBgVal = 0.2;
% stimRVal = 0.8;
% stimGVal = 0.65;
% stimBVal = 0.2;
stimRVal = 0.5;
stimGVal = 0.5;
stimBVal = 0.5;

% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional
%                           display.
sparsePriorStr = 'conventional';
sparsePriorName = [sparsePriorStr 'SparsePrior.mat'];

% Forward rendering parameters
%
% Use AO in forward rendering?
% This determins pupil diameter which typically differs in AO
forwardAORender = true;
if (forwardAORender)
    forwardPupilDiamMM = 7;
    forwardAOStr = 'AO';
else
    forwardPupilDiamMM = 3;
    forwardAOStr = 'NOAO';
end

% Residual defocus for forward rendering
forwardDefocusDiopters = 0.05;

% Force build and save
buildNewForward = false;

% Recon rendering parameters
useForwardRenderingForRecon = true;
reconstructfromRenderMatrix = true;

%% Set forward render filenname
if (forwardAORender)
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters);
else
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters);
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
        % We build a normal optics structure, and then overwrite the PSF
        % property.  Allow specified defocus.
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', false);
        theConeMosaic.PSF = ConeResponse.psfDiffLmt(forwardPupilDiamMM,'defocusDiopters', forwardDefocusDiopters);
    else
        % Build normal optics structure. 
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', forwardPupilDiamMM, 'useRandomSeed', false, ...
            'defocusDiopters',forwardDefocusDiopters);
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
    reconAOStr = forwardAOStr;
else
    error('Need to implement separate recon rendering setup');
end

%% Setup output directory
outputName = sprintf('%s_%s_%d_%0.2f_%0.2f_%s_%s_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f', ...
    forwardAOStr,reconAOStr,nPixels,forwardDefocusDiopters,reconDefocusDiopters,displayName, sparsePriorStr, ...
    stimSizeDegs,stimBgVal,stimRVal,stimGVal,stimBVal);
outputDir = fullfile(aoReconDir,outputName);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Need new render if we want to reconstruct with respect to the AO stimulus
% reconstructWrtAO = true;
% buildNewAO = false;
% if (reconstructWrtAO)
%     if (buildNewAO)
%         forwardRenderMatrix = theConeMosaic.forwardRender([nPixels nPixels 3], ...
%             'validation', false);
%         forwardRenderMatrix = double(forwardRenderMatrix);
%     else
%         load(fullfile(aoReconDir,'aoOpticsRenderMatrix.mat'));
%     end
% end

%% Show forward cone mosaic
forwardConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(outputDir,'Mosaic.jpg'),'jpg');

%% Generate an image stimulus
% stimulus in the size of retinal degree
% should not exceed 'fieldSizeDegs'
stimSizeFraction = stimSizeDegs / fieldSizeDegs;
idxLB = round(nPixels * (0.5 - stimSizeFraction / 2));
idxUB = round(nPixels * (0.5 + stimSizeFraction / 2));
idxRange = idxLB:idxUB;

% Image stimulus with a gray background and yellow color
stimulusImageRGB = ones(nPixels, nPixels, 3) * stimBgVal;
stimulusImageRGB(idxRange, idxRange, 1) = stimRVal;
stimulusImageRGB(idxRange, idxRange, 2) = stimGVal;
stimulusImageRGB(idxRange, idxRange, 3) = stimBVal;

% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false);
saveas(gcf,fullfile(outputDir,'Stimulus.jpg'),'jpg');

%% Compute forward retinal image and excitations using ISETBio
%
% We'll reconstruct from these.
forwardOI = oiCompute(stimulusScene,forwardOI);
visualizeOpticalImage(forwardOI);
forwardExcitationsToStimulusISETBio = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));

% Check forward exciations calculation another way.  Also shows another way
% to visualize the retinal image, but this only is done when the check
% fails.
forwardExcitationsToStimulusCheck = forwardConeMosaic.compute(stimulusImageRGB);
if (max(abs(forwardExcitationsToStimulusCheck-forwardExcitationsToStimulusISETBio)) ~= 0)
    forwardConeMosaic.visualizeOI()
    error('Two ways of doing the same thing do not agree');
end
temp = squeeze(forwardConeMosaic.LastResponse);
if (max(abs(temp-forwardExcitationsToStimulusISETBio)) ~= 0)
    error('Last excitations in object not as we expect');
end

%% Compute excitations using forward render matrix.
%
% We would like these to agree with direct rendering, but
% for reasons we are slowly coming to understand (quantization
% is a non-linear process), they don't always.
forwardExcitationsToStimulusRenderMatrix = forwardRenderMatrix*stimulusImageLinear(:);
figure; clf; hold on;
plot(forwardExcitationsToStimulusISETBio,forwardExcitationsToStimulusRenderMatrix,'ro','MarkerFaceColor','r','MarkerSize',10);
axis('square');
maxVal = max([forwardExcitationsToStimulusISETBio; forwardExcitationsToStimulusRenderMatrix]);
plot([0 maxVal],[0 maxVal],'k');
xlim([0 maxVal]); ylim([0 maxVal]);
xlabel('Excitations to stimulus ISETBio');
ylabel('Excitations to stimulus render matrix');
title('Exciations ISETBio and render matrix');
saveas(gcf,fullfile(outputDir,'ISETBioVsRenderMatrixExciations.jpg'),'jpg');

%% Choose which excitations to reconstruct form
if (reconstructfromRenderMatrix)
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusRenderMatrix;
else
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusISETBio;
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
if (reconstructfromRenderMatrix)
    titleStr = 'Excitations using render matrix';
else
    titleStr = 'Excitations using ISETBio';
end
figureHandle = figure(); axesHandle = [];
forwardConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', axesHandle, ...
    'activation', reshape(forwardExcitationsToStimulusUse,1,1,length(forwardExcitationsToStimulusUse)), ...
    'activationRange', [0 max(forwardExcitationsToStimulusUse)], ...
    'plotTitle',  titleStr);
saveas(gcf,fullfile(outputDir,'MosaicExcitations.jpg'),'jpg');

%% Run reconstruction
%
% Load prior
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
reconImage = estimator.runEstimate(forwardExcitationsToStimulusUse * scaleFactor, ...
    'maxIter', 500, 'display', 'iter', 'gpu', false);

% Show reconstruction
meanLuminanceCdPerM2 = [];
[reconScene, ~, reconImageLinear] = sceneFromFile(gammaCorrection(reconImage, theDisplay), 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
reconScene = sceneSet(reconScene, 'fov', fieldSizeDegs);
visualizeScene(reconScene, 'displayRadianceMaps', false);
saveas(gcf,fullfile(outputDir,'Recon.jpg'),'jpg');

% Compute forward excitations from reconstruction
% And compare with stimulus exciations
forwardOI = oiCompute(reconScene,forwardOI);
forwardExcitationsToRecon = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));
figure; clf; hold on;
plot(forwardExcitationsToStimulusUse,forwardExcitationsToRecon,'ro','MarkerFaceColor','r','MarkerSize',10);
axis('square');
maxVal = max([forwardExcitationsToStimulusUse; forwardExcitationsToRecon]);
plot([0 maxVal],[0 maxVal],'k');
xlim([0 maxVal]); ylim([0 maxVal]);
xlabel('Excitations to stimulus');
ylabel('Excitations to reconstruction');
if (reconstructfromRenderMatrix)
    title('Reconstruction from forward render matrix');
else
    title('Reconstruction from forward ISETBio');
end
saveas(gcf,fullfile(outputDir,'StimulusVsReconExcitations.jpg'),'jpg');

%% Evaluate prior and likelihood of stimulus and reconstruction
[stimNegLogPrior,~,stimNegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, stimulusImageLinear(:));
[reconNegLogPrior,~,reconNegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, reconImageLinear(:));
fprintf('Stimulus: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
    -stimNegLogPrior,-stimNegLogLikely,-(stimNegLogPrior+stimNegLogLikely));
fprintf('Recon: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
    -reconNegLogPrior,-reconNegLogLikely,-(reconNegLogPrior+reconNegLogLikely));
fprintf('Each of the folowing should be *higher* for a valid reconstruction');
if (-stimNegLogPrior > -reconNegLogPrior)
    fprintf('\tReconstruction prior *lower* than stimulus\n');
else
    fprintf('\tReconstruction prior *higher* than stimulus\n');
end
if (-stimNegLogLikely > -reconNegLogLikely)
    fprintf('\tStimulus likelihood *higher* than reconstruction\n');
else
    fprintf('\tStimulus likelihood *lower* than reconstruction\n');
end
if (-(stimNegLogPrior+stimNegLogLikely) > -(reconNegLogPrior+reconNegLogLikely))
    fprintf('\tReconstruction neg objective *lower* than stimulus\n');
else
    fprintf('\ttReconstruction neg objective *higher* than stimulus\n');
end

%% Save workspace
close all;
save(fullfile(outputDir,'xRunOutput.mat'));
