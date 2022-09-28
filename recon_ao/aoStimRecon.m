function aoStimRecon(displayName,sparsePriorStr,...
    forwardAORender, reconAORender, ...
    forwardDefocusDiopters, reconDefocusDiopters, ...
    stimSizeDegs,stimBgVal,stimRVal,stimGVal,stimBVal,...
    regPara, stride, forwardChrom, reconChrom, newCenter, nPixels, centerPixel)
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
%   08/19/22  dhb, chr  Edit to clarify and remove stochasticity 
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options
%   09/22/22  chr  Convert to its own dichrom file

%% Close existing figures
close all;

%% Point at directory with data files for this subproject
%
% This will allow us to load in project specific precomputed information.
% Also records initials of version editors, otherwise set to 'main'
aoReconDir = getpref('ISETImagePipeline','aoReconDir');

% Replaced versEditor to apply it to subDirectories as well
helpDir = '/helperFiles';
versEditor = 'dichromTests';

%% Setup / Simulation parameters
%
% Spatial parameters.  Common to forward and recon models
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

% Determine which method will be used for the reconstruction: ISETBIO or
% Render Matrix
reconstructfromRenderMatrix = true;
if (reconstructfromRenderMatrix)
    exciteSource = 'renderMatrix';
else
    exciteSource = 'isetbio';
end

% Establish if a random seed will be used when building the cone mosaic for
% the render matrices
forwardRandSeed = false;
reconRandSeed = false;
if (forwardRandSeed)
    forwardSeedStr = 'rand';
else
    forwardSeedStr = 'noRand';
end
if (reconRandSeed)
    reconSeedStr = 'rand';
else
    reconSeedStr = 'noRand';
end

% Optional chromaticities to be used in forward and recon cone mosaics,
% with normal being trichromatic. Can consider bumping this off to another
% function to declutter, return the replaceCones, startCones and newCones
[replaceForwardCones, forwardStartCones, ...
    forwardNewCones] = assignCones(forwardChrom);
[replaceReconCones, reconStartCones, ...
    reconNewCones] = assignCones(reconChrom);

%% Set render filennames
if (forwardAORender)
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters, forwardSeedStr, forwardChrom);
else
    forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters, forwardSeedStr, forwardChrom);
end
if (reconAORender)
    reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,reconPupilDiamMM,reconDefocusDiopters, reconSeedStr, reconChrom);
else
    reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,reconPupilDiamMM,reconDefocusDiopters, reconSeedStr, reconChrom);
end

%% Build render matrices/files or load from existing cache

% Build or grab foward cone mosaic and render 
if (buildNewForward || ~exist(fullfile(aoReconDir, helpDir, versEditor, forwardRenderStructureName),'file'))
    renderStructure = buildRenderStruct_dichrom(aoReconDir, helpDir, eccXDegs, eccYDegs, ...
        fieldSizeDegs, nPixels, forwardPupilDiamMM, forwardAORender, forwardDefocusDiopters, ...
        overwriteDisplayGamma, displayName, displayFieldName, displayGammaBits, ...
        displayGammaGamma, forwardRandSeed, replaceForwardCones, forwardStartCones, forwardNewCones);
    save(fullfile(aoReconDir, helpDir, versEditor, forwardRenderStructureName),'renderStructure');
    forwardRenderStructure = renderStructure; clear renderStructure;
else
    clear forwardRenderStructure;
    load(fullfile(aoReconDir, helpDir, versEditor, forwardRenderStructureName),'renderStructure');
    forwardRenderStructure = renderStructure; clear renderStructure; 
    grabRenderStruct(forwardRenderStructure, eccXDegs, eccYDegs, fieldSizeDegs, ...
        nPixels, forwardPupilDiamMM, forwardAORender, forwardDefocusDiopters)
end

% Build or grab recon cone mosaic and render 
if (buildNewRecon || ~exist(fullfile(aoReconDir, helpDir, versEditor, reconRenderStructureName),'file'))
    renderStructure = buildRenderStruct_dichrom(aoReconDir, helpDir, eccXDegs, eccYDegs, ...
        fieldSizeDegs, nPixels, reconPupilDiamMM, reconAORender, reconDefocusDiopters, ...
        overwriteDisplayGamma, displayName, displayFieldName, displayGammaBits, ...
        displayGammaGamma, reconRandSeed, replaceReconCones, reconStartCones, reconNewCones);
    save(fullfile(aoReconDir, helpDir, versEditor, reconRenderStructureName),'renderStructure');
    reconRenderStructure = renderStructure; clear renderStructure;
else
    clear reconRenderStructure;
    load(fullfile(aoReconDir, helpDir, versEditor, reconRenderStructureName),'renderStructure');
    reconRenderStructure = renderStructure; clear renderStructure; 
    grabRenderStruct(reconRenderStructure, eccXDegs, eccYDegs, fieldSizeDegs, ...
        nPixels, reconPupilDiamMM, reconAORender, reconDefocusDiopters)
end

% Set forward variables from loaded/built structure
forwardRenderMatrix = forwardRenderStructure.renderMatrix;
forwardConeMosaic = forwardRenderStructure.theConeMosaic;
forwardOI = forwardConeMosaic.PSF;

% Set recon variables from loaded/built structure
reconRenderMatrix = reconRenderStructure.renderMatrix;
reconConeMosaic = reconRenderStructure.theConeMosaic;
% reconOI = reconConeMosaic.PSF;

% Clear forward render structure
clear forwardRenderStructure;
clear reconRenderStructure;

%% Setup output directories
outputMainName = sprintf('%s_%s_%0.2f_%0.2f_%d_%d_%s_%s_%s', ...
    forwardAOStr,reconAOStr,forwardDefocusDiopters,reconDefocusDiopters,nPixels,fieldSizeMinutes,displayName,sparsePriorStr, versEditor);
outputSubName = sprintf('%0.1f_%0.4f_%d_%0.2f_%0.2f_%0.2f_%0.2f_%s_%s_%s',60*stimSizeDegs, regPara,stride,stimBgVal,stimRVal,stimGVal,stimBVal, exciteSource, forwardChrom, reconChrom);
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

% Shift the stimulus to be centered on desired values
idxXRange = (idxLB:idxUB) + newCenter(1);
idxYRange = (idxLB:idxUB) + newCenter(2);

% Check stimulus position. Ends function and deletes the created
% directory if the stimulus position exceeds bounds. 
if min(idxYRange) <= 0 || max(idxYRange) > nPixels ...
        || min(idxXRange) <= 0 || max(idxXRange) > nPixels
    warning(['Stimulus centered on ' int2str(newCenter' + centerPixel) ...
        ' exceeds bounds. Beginning next simulation']);
    rmdir(outputDir, 's');
    close all
    return
end

% Set image pixels
stimulusImageRGB = ones(nPixels, nPixels, 3) * stimBgVal;
stimulusImageRGB(idxYRange, idxXRange, 1) = stimRVal;
stimulusImageRGB(idxYRange, idxXRange, 2) = stimGVal;
stimulusImageRGB(idxYRange, idxXRange, 3) = stimBVal;

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

% Compare with ISETBio
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

%% Run reconstruction
%
% Load prior
prior = load(fullfile(aoReconDir, helpDir, sparsePriorName));

% Construct image estimator
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

maxReconIterations = 500;
specifiedStarts = {};
specifiedStarts{1} = 0.5*ones(length(stimulusImageLinear(:)), 1);
[multistartStruct,~,reconIndex] = estimator.runMultistartEstimate(forwardExcitationsToStimulusUse * scaleFactor, ...
    'maxIter', maxReconIterations, 'display', 'iter', 'gpu', false, ...
    'nWhiteStart', 1, 'nPinkStart', 1, ...
    'nSparsePriorPatchStart', 1, 'sparsePrior', prior, ...
    'specifiedStarts', specifiedStarts);

% Diagnose reconstructions
% txtFileName = fullfile(outputDir,'ReconProbInfo.txt');
% if (exist(txtFileName,'file'))
%     delete(txtFileName);
% end
% fid = fopen(txtFileName,'w');

% Evaluate stimulus
[stimNegLogPrior,~,stimNegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, stimulusImageLinear(:));
stimLoss = stimNegLogPrior + stimNegLogLikely;

for ii = 1:length(multistartStruct.initTypes)
    [initSceneTemp, ~, initImageLinearTemp] = sceneFromFile(gammaCorrection(multistartStruct.initImages{ii}, forwardConeMosaic.Display), 'rgb', ...
        meanLuminanceCdPerM2, forwardConeMosaic.Display);
    [reconSceneTemp, ~, reconImageLinearTemp] = sceneFromFile(gammaCorrection(multistartStruct.reconImages{ii}, forwardConeMosaic.Display), 'rgb', ...
        meanLuminanceCdPerM2, forwardConeMosaic.Display);
    theFig = figure; clf;
    set(theFig,'Position',[300 400 1150 780]);
    theAxes = subplot(3,2,1);
    visualizeScene(initSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    theAxes = subplot(3,2,2);
    visualizeScene(reconSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    subplot(3,2,3); hold on;
    minVal = 0.9*min([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)]);
    maxVal = 1.1*max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)]);
    plot(multistartStruct.coneVec,multistartStruct.reconPreds(:,ii),'ro','MarkerFaceColor','r','MarkerSize',6);
    xlabel('Measured Excitations');
    ylabel('Predicted Exciations');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    axis('square');
    title(sprintf('Recon %d, init %s',ii,multistartStruct.initTypes{ii}));

    subplot(3,2,4);
    bar([1]', ...
        [multistartStruct.initLogPriors(ii)  ; ...
         multistartStruct.reconLogPriors(ii) ; ...
         -stimNegLogPrior]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Log Prior');
    axis('square');
    title('Init/Recon/Stim Log Priors');
    subplot(3,2,5);
    bar([1]', ...
        [multistartStruct.initLogLikelihoods(ii)  ; ...
         multistartStruct.reconLogLikelihoods(ii) ; ...
         -stimNegLogLikely]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Log Likelihood');
    axis('square');
    title('Init/Recon/Stim Log Likelihoods');
    subplot(3,2,6);
    bar([1]', ...
        [-multistartStruct.initLosses(ii)  ; ...
         -multistartStruct.reconLosses(ii) ; ...
         -stimLoss]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Neg Loss');
    axis('square');
    if (multistartStruct.reconLosses(ii) < stimLoss)
        if (multistartStruct.reconLosses(ii) < multistartStruct.initLosses(ii))
            title({'Init/Recon/Stim Neg Losses' ; 'Recon BETTER than Stim' ; 'Recon BETTER than Init'});
        else
            title({'Init/Recon/Stim Neg Losses' ; 'Recon BETTER than Stim' ; 'Recon WORSE than Init'});
        end
    else
        if (multistartStruct.reconLosses(ii) < multistartStruct.initLosses(ii))
            title({'Init/Recon/Stim Neg Losses' ; 'Recon WORSE than Stim' ; 'Recon BETTER than Init'});
        else
            title({'Init/Recon/Stim Neg Losses' ; 'Recon WORSE than Stim' ; 'Recon WORSE than Init'});
        end
    end

    % Save
    saveas(gcf,fullfile(outputDir,sprintf('Recon%dSummary.jpg',ii)),'jpg');
    if (ii == reconIndex)
        saveas(gcf,fullfile(outputDir,sprintf('ReconSummary.jpg',ii)),'jpg');
    end
end
% fclose(fid);

% Save best
[reconScene, ~, reconImageLinear] = sceneFromFile(gammaCorrection(multistartStruct.reconImages{reconIndex}, forwardConeMosaic.Display), 'rgb', ...
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

% fprintf(fid,'Stimulus: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
%     -stimNegLogPrior,-stimNegLogLikely,-(stimNegLogPrior+stimNegLogLikely));
% fprintf(fid,'Recon1: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
%     -recon1NegLogPrior,-recon1NegLogLikely,-(recon1NegLogPrior+recon1NegLogLikely));
% fprintf(fid,'Recon1 initial loss %0.6g; recon1 solution loss %0.6g; fractional difference (init less soln; should be pos): %0.6g\n', ...
%     recon1InitLoss,recon1SolnLoss,(recon1InitLoss-recon1SolnLoss)/abs(recon1InitLoss));
% fprintf(fid,'Recon2: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
%     -recon2NegLogPrior,-recon2NegLogLikely,-(recon2NegLogPrior+recon2NegLogLikely));
% fprintf(fid,'Recon2 initial loss %0.6g; recon2 solution loss %0.6g; fractional difference (init less soln; should be pos): %0.6g\n', ...
%     recon2InitLoss,recon2SolnLoss,(recon2InitLoss-recon2SolnLoss)/abs(recon2InitLoss));
% fprintf(fid,reconWhichStr);
% fprintf(fid,'Each of the following should be *higher* for a valid reconstruction\n');
% if (-stimNegLogPrior > -reconNegLogPrior)
%     fprintf(fid,'\tReconstruction prior *lower* than stimulus\n');
% else
%     fprintf(fid,'\tReconstruction prior *higher* than stimulus\n');
% end
% if (-stimNegLogLikely > -reconNegLogLikely)
%     fprintf(fid,'\tStimulus likelihood *higher* than reconstruction\n');
% else
%     fprintf(fid,'\tStimulus likelihood *lower* than reconstruction\n');
% end
% if (-(stimNegLogPrior+stimNegLogLikely) > -(reconNegLogPrior+reconNegLogLikely))
%     fprintf(fid,'\tReconstruction neg objective *lower* than stimulus\n');
% else
%     fprintf(fid,'\tReconstruction neg objective *higher* than stimulus\n');
% end

%% Save workspace
close all;
save(fullfile(outputDir,'xRunOutput.mat'));
