function aoStimRecon_chr(displayName,sparsePriorStr,...
    forwardAORender, reconAORender, ...
    forwardDefocusDiopters, reconDefocusDiopters, ...
    stimSizeDegs,stimBgVal,stimRVal,stimGVal,stimBVal,...
    regPara,stride)
% Synopsis:
%    Driver to run AO recon simulations.
%
% Description:
%    Script to see how well we can measure unique yellow percepts under AO
%    conditions.  This has many parameters.  See aoStimReconRunMany for
%    those not fixed here.
%
%    This script organizes and saves ts output in the directory hierarchy
%    set up by the local hook file.
%
% See also: aoStimReconRunMany

% History:
%   07/29/22  lz    Wrote this sometime in the past
%   07/29/22  dhb, chr  Starting to dig into this
%   08/13/22  dhb   Lots of bookkeeping, cleaning, etc.
%   08/14/22  dhb   Made it a callable function.
%   08/17/22  chr   Incorporated separate forward and recon callings 

%% Close existing figures
% close all;

%% Point at directory with data files for this subproject
%
% This will allow us to load in project specific precomputed information.
% Also records initials of version editors, otherwise set to 'main'
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
versEditor = 'TESTING2';

%% Setup / Simulation parameters
%
% Spatial parameters.  Common to forward and recon models
nPixels = 58;
fieldSizeMinutes = 30;
fieldSizeDegs = fieldSizeMinutes/60;
eccXDegs = 2.0;
eccYDegs = 0.0;

% Gamma table parameters
displayGammaBits = 12;
displayGammaGamma = 2;
switch (displayName)
    case 'conventional'
        displayFieldName = 'CRT12BitDisplay';
        overwriteDisplayGamma = false;
    case 'mono'
        displayFieldName = 'monoDisplay';
        overwriteDisplayGamma = true;
    otherwise
        error('Unknown display specified');
end

% Sparse prior name
sparsePriorName = [sparsePriorStr 'SparsePrior.mat'];

% Determine forward pupil diameter, allowing it to differ in AO case
if (forwardAORender)
    forwardPupilDiamMM = 7;
    forwardAOStr = ['AO' num2str(forwardPupilDiamMM)];
else
    forwardPupilDiamMM = 3;
    forwardAOStr = ['NOAO' num2str(forwardPupilDiamMM)];
end

if (reconAORender)
    reconPupilDiamMM = 7;
    reconAOStr = ['AO' num2str(reconPupilDiamMM)];
else
    reconPupilDiamMM = 3;
    reconAOStr = ['NOAO' num2str(reconPupilDiamMM)];
end



% Force build and save
buildNewForward = false;
buildNewRecon = false;



% Recon rendering parameters
% 
% DHB: See "DHB:" comment below. I think we no
% longer need these two "use" variables.
% useForwardRenderingForRecon = true;
% if (useForwardRenderingForRecon)
%     forwardRandSeed = true;
% end
% useReconRenderingForRecon = true;
% if (useReconRenderingForRecon)
%     reconRandSeed = true;
% end

reconstructfromRenderMatrix = true;
forwardRandSeed = false;
reconRandSeed = true;

if (forwardRandSeed)
    forwardSeed = 'rand';
else
    forwardSeed = 'noRand';
end
if (reconRandSeed)
    reconSeed = 'rand';
else
    reconSeed = 'noRand';
end

%% Set render filennames
if (forwardAORender)
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters, forwardSeed);
else
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters, forwardSeed);
end
if (reconAORender)
    reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,reconPupilDiamMM,reconDefocusDiopters, reconSeed);
else
    reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,reconPupilDiamMM,reconDefocusDiopters, reconSeed);
end

%% Build render matrices/files or load from existing cache
%
% Have to do a workaround if loading render structures (line 141) since as of now
% all of them use the name 'forwardRenderStructure'. Instead of having to
% rebuild, make the call for the forwardRender and then manually assign
% that to the reconRenderStructure (line 142). Then in the second if statement replace
% the forwardRenderStructure with what it should actaully be. 
if (buildNewRecon || ~exist(fullfile(aoReconDir,reconRenderStructureName),'file'))
reconRenderStructure = buildRenderStruct_chr(aoReconDir, eccXDegs, eccYDegs, ...
    fieldSizeDegs, nPixels, reconPupilDiamMM, reconAORender, reconDefocusDiopters, ...
    overwriteDisplayGamma, displayName, displayFieldName, displayGammaBits, ...
    displayGammaGamma, reconRandSeed);
save(fullfile(aoReconDir,reconRenderStructureName),'reconRenderStructure');
else
    clear reconRenderStructure;
    load(fullfile(aoReconDir,reconRenderStructureName),'forwardRenderStructure');
    reconRenderStructure = forwardRenderStructure;
    grabRenderStruct_chr(reconRenderStructure, eccXDegs, eccYDegs, fieldSizeDegs, ... 
        nPixels, reconPupilDiamMM, reconAORender, reconDefocusDiopters)
end 

if (buildNewForward || ~exist(fullfile(aoReconDir,forwardRenderStructureName),'file'))
forwardRenderStructure = buildRenderStruct_chr(aoReconDir, eccXDegs, eccYDegs, ...
    fieldSizeDegs, nPixels, forwardPupilDiamMM, forwardAORender, forwardDefocusDiopters, ...
    overwriteDisplayGamma, displayName, displayFieldName, displayGammaBits, ...
    displayGammaGamma, forwardRandSeed);
save(fullfile(aoReconDir,forwardRenderStructureName),'forwardRenderStructure');
else
    clear forwardRenderStructure;
    load(fullfile(aoReconDir,forwardRenderStructureName),'forwardRenderStructure');
    grabRenderStruct_chr(forwardRenderStructure, eccXDegs, eccYDegs, fieldSizeDegs, ... 
        nPixels, forwardPupilDiamMM, forwardAORender, forwardDefocusDiopters)
end

% Set forward variables from loaded/built structure
forwardRenderMatrix = forwardRenderStructure.renderMatrix;
forwardConeMosaic = forwardRenderStructure.theConeMosaic;
forwardOI = forwardConeMosaic.PSF;

% Set recon variables from loaded/built structure
reconRenderMatrix = reconRenderStructure.renderMatrix;
reconConeMosaic = reconRenderStructure.theConeMosaic;
reconOI = reconConeMosaic.PSF;

% Set display variable
theForwardDisplay = forwardRenderStructure.theDisplay;
theReconDisplay = reconRenderStructure.theDisplay;

% Clear forward render structure
clear forwardRenderStructure;
clear reconRenderStructure;

%% Reconstruction rendering parameters
%
% These are how the reconstruction algorithm thinks the image
% was formed.  Need not be the same as the forward rendering.
%
% DHB: Now that you have the ability to create and load 
% separate forward and recon files, we don't really need
% a conditinal here.  Rather, just proceed as above where
% you set the forward and recon variables appropriately and
% get rid of this conditional altogether.
% if (useForwardRenderingForRecon)
%     reconRenderMatrix = forwardRenderMatrix;
%     reconConeMosaic = forwardConeMosaic;
%     reconOI = forwardOI;
%     reconPupilDiamMM = forwardPupilDiamMM;
%     reconAORender = forwardAORender;
%     reconDefocusDiopters = forwardDefocusDiopters;
%     reconAOStr = forwardAOStr;
% elseif (useReconRenderingForRecon)
% else % ----------------------------------
%     error('Need to implement separate recon rendering setup');
% end


%% Setup output directories
outputMainName = sprintf('%s_%s_%0.2f_%0.2f_%d_%d_%0.1f_%s_%s_%s', ...
    forwardAOStr,reconAOStr,forwardDefocusDiopters,reconDefocusDiopters,nPixels,fieldSizeMinutes,60*stimSizeDegs,displayName,sparsePriorStr, versEditor);
outputSubName = sprintf('%0.4f_%d_%0.2f_%0.2f_%0.2f_%0.2f',regPara,stride,stimBgVal,stimRVal,stimGVal,stimBVal);
outputDir = fullfile(aoReconDir,outputMainName,outputSubName);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Show forward cone mosaic
forwardConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(outputDir,'forwardMosaic.jpg'),'jpg');

reconConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(outputDir,'reconMosaic.jpg'),'jpg');

%% Generate an image stimulus
%
% Stimulus size in retinal degrees should not exceed 'fieldSizeDegs'
stimSizeFraction = stimSizeDegs / fieldSizeDegs;
if (stimSizeFraction > 1)
    error('Stimulus size too big given field size');
end
idxLB = round(nPixels * (0.5 - stimSizeFraction / 2));
idxUB = round(nPixels * (0.5 + stimSizeFraction / 2));
idxRange = idxLB:idxUB;

% Set image pixels
stimulusImageRGB = ones(nPixels, nPixels, 3) * stimBgVal;
stimulusImageRGB(idxRange, idxRange, 1) = stimRVal;
stimulusImageRGB(idxRange, idxRange, 2) = stimGVal;
stimulusImageRGB(idxRange, idxRange, 3) = stimBVal;

% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(outputDir,'Stimulus.jpg'),'jpg');

%% Compute forward retinal image and excitations using ISETBio
%
% We may or may not reconstruct from these
forwardOI = oiCompute(stimulusScene,forwardOI);
visualizeOpticalImage(forwardOI, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(outputDir,'forwardStimulusRetinalImage.jpg'),'jpg');
forwardExcitationsToStimulusISETBio = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));

reconOI = oiCompute(stimulusScene,forwardOI);
visualizeOpticalImage(reconOI, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(outputDir,'reconStimulusRetinalImage.jpg'),'jpg');
reconExcitationsToStimulusISETBio = squeeze(reconConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));

% Check forward exciations calculation another way.  Also shows another way
% to visualize the retinal image, but this only is done when the check
% fails.
forwardExcitationsToStimulusCheck = forwardConeMosaic.compute(stimulusImageRGB);
reconExcitationsToStimulusCheck = reconConeMosaic.compute(stimulusImageRGB);


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
reconExcitationsToStimulusRenderMatrix = reconRenderMatrix*stimulusImageLinear(:);



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
    reconExcitationsToStimulusUse = reconExcitationsToStimulusRenderMatrix;
else
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusISETBio;
    reconExcitationsToStimulusUse = reconExcitationsToStimulusISETBio;
end

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
saveas(gcf,fullfile(outputDir,'forwardMosaicExcitations.jpg'),'jpg');

figureHandle = figure(); axesHandle = [];
reconConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', axesHandle, ...
    'activation', reshape(reconExcitationsToStimulusUse,1,1,length(reconExcitationsToStimulusUse)), ...
    'activationRange', [0 max(reconExcitationsToStimulusUse)], ...
    'plotTitle',  titleStr);
saveas(gcf,fullfile(outputDir,'reconMosaicExcitations.jpg'),'jpg');

%% Run reconstruction
%
% Load prior
prior = load(fullfile(aoReconDir,sparsePriorName));

% Construct onstruct image estimator
estimator = PoissonSparseEstimator(reconRenderMatrix, inv(prior.regBasis), ...
    prior.mu', regPara, stride, [nPixels nPixels 3]);

% Estimate
%
% Scale excitations to take into account difference in forward and recon
% pupil sizes.  This helps keep things in range.
%
% DHB: I think I might have scale factor backwards, now that I look at it.
% If we use a large forward pupil and reconstruct with a small one, we want
% to scale the cone exciations passed to the recon to be smaller.  So it
% should be (reconPupilDiamMM)^2/(fowardPupilDiamMM)^2.  Leaving it for
% now, but as soon as the two pupil sizes differ, we'll want to check
% carefully.
meanLuminanceCdPerM2 = [];
scaleFactor = (forwardPupilDiamMM/reconPupilDiamMM)^2;
[recon1Image,recon1InitLoss,recon1SolnLoss] = estimator.runEstimate(forwardExcitationsToStimulusUse * scaleFactor, ...
    'maxIter', 500, 'display', 'iter', 'gpu', false);
[recon1Scene, ~, recon1ImageLinear] = sceneFromFile(gammaCorrection(recon1Image, theForwardDisplay), 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
recon1Scene = sceneSet(recon1Scene, 'fov', fieldSizeDegs);
[recon1NegLogPrior,~,recon1NegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, recon1ImageLinear(:));
visualizeScene(recon1Scene, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(outputDir,'Recon1.jpg'),'jpg');

% Start at stimulus image.  We do this to check if the random starting
% point gets stuck in a local minimum.  If we find much of this, we will
% need to think of a better way to initialize but that doesn't depend on
% knowing the stimulus.  Maybe start at the prior mean?
[recon2Image,recon2InitLoss,recon2SolnLoss] = estimator.runEstimate(forwardExcitationsToStimulusUse * scaleFactor, ...
    'maxIter', 500, 'display', 'iter', 'gpu', false, 'init', stimulusImageLinear(:));
[recon2Scene, ~, recon2ImageLinear] = sceneFromFile(gammaCorrection(recon2Image, theForwardDisplay), 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
recon2Scene = sceneSet(recon2Scene, 'fov', fieldSizeDegs);
[recon2NegLogPrior,~,recon2NegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, recon2ImageLinear(:));
visualizeScene(recon2Scene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(outputDir,'Recon2.jpg'),'jpg');

% Report back the better
if (-(recon1NegLogPrior+recon1NegLogLikely) > -(recon2NegLogPrior+recon2NegLogLikely))
    reconWhichStr = 'Recon1 (random start) better\n';
    reconNegLogPrior = recon1NegLogPrior;
    reconNegLogLikely = recon1NegLogLikely;
    reconImage = recon1Image;
    reconScene = recon1Scene;
    reconImageLinear = recon1ImageLinear;
else
    reconWhichStr = 'Recon2 (stimulus start) better\n';
    reconNegLogPrior = recon2NegLogPrior;
    reconNegLogLikely = recon2NegLogLikely;
    reconImage = recon2Image;
    reconScene = recon2Scene;
    reconImageLinear = recon2ImageLinear;
end

% Show reconstruction
[reconScene, ~, reconImageLinear] = sceneFromFile(gammaCorrection(recon2Image, theForwardDisplay), 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
reconScene = sceneSet(reconScene, 'fov', fieldSizeDegs);
visualizeScene(reconScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(outputDir,'Recon.jpg'),'jpg');

% Compute forward excitations from reconstruction
% And compare with stimulus exciations
forwardOI = oiCompute(reconScene,forwardOI);
if (reconstructfromRenderMatrix)
    title('Reconstruction from forward render matrix');
    forwardExcitationsToRecon = squeeze(forwardRenderMatrix*reconImageLinear(:));
else
    title('Reconstruction from forward ISETBio');
    forwardExcitationsToRecon = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));
end
figure; clf; hold on;
plot(forwardExcitationsToStimulusUse,forwardExcitationsToRecon,'ro','MarkerFaceColor','r','MarkerSize',10);
axis('square');
maxVal = max([forwardExcitationsToStimulusUse; forwardExcitationsToRecon]);
plot([0 maxVal],[0 maxVal],'k');
xlim([0 maxVal]); ylim([0 maxVal]);
xlabel('Excitations to stimulus');
ylabel('Excitations to reconstruction');
saveas(gcf,fullfile(outputDir,'StimulusVsReconExcitations.jpg'),'jpg');

%% Evaluate prior and likelihood of stimulus and reconstruction
[stimNegLogPrior,~,stimNegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, stimulusImageLinear(:));
txtFileName = fullfile(outputDir,'ReconProbInfo.txt');
if (exist(txtFileName,'file'))
    delete(txtFileName);
end
fid = fopen(txtFileName,'w');
fprintf(fid,'Stimulus: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
    -stimNegLogPrior,-stimNegLogLikely,-(stimNegLogPrior+stimNegLogLikely));
fprintf(fid,'Recon1: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
    -recon1NegLogPrior,-recon1NegLogLikely,-(recon1NegLogPrior+recon1NegLogLikely));
fprintf(fid,'Recon1 initial loss %0.6g; recon1 solution loss %0.6g; fractional difference (init less soln; should be pos): %0.6g\n', ...
    recon1InitLoss,recon1SolnLoss,(recon1InitLoss-recon1SolnLoss)/abs(recon1InitLoss));
fprintf(fid,'Recon2: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
    -recon2NegLogPrior,-recon2NegLogLikely,-(recon2NegLogPrior+recon2NegLogLikely));
fprintf(fid,'Recon2 initial loss %0.6g; recon2 solution loss %0.6g; fractional difference (init less soln; should be pos): %0.6g\n', ...
    recon2InitLoss,recon2SolnLoss,(recon2InitLoss-recon2SolnLoss)/abs(recon2InitLoss));
fprintf(fid,reconWhichStr);
fprintf(fid,'Each of the following should be *higher* for a valid reconstruction\n');
if (-stimNegLogPrior > -reconNegLogPrior)
    fprintf(fid,'\tReconstruction prior *lower* than stimulus\n');
else
    fprintf(fid,'\tReconstruction prior *higher* than stimulus\n');
end
if (-stimNegLogLikely > -reconNegLogLikely)
    fprintf(fid,'\tStimulus likelihood *higher* than reconstruction\n');
else
    fprintf(fid,'\tStimulus likelihood *lower* than reconstruction\n');
end
if (-(stimNegLogPrior+stimNegLogLikely) > -(reconNegLogPrior+reconNegLogLikely))
    fprintf(fid,'\tReconstruction neg objective *lower* than stimulus\n');
else
    fprintf(fid,'\tReconstruction neg objective *higher* than stimulus\n');
end
fclose(fid);

%% Save workspace
close all;
save(fullfile(outputDir,'xRunOutput.mat'));
