function aoStimRecon(pr,cnv)
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
if (~exist(cnv.renderDir ,'dir'))
    mkdir(cnv.renderDir );
end

% Sparse prior name
sparsePriorName = [pr.sparsePriorStr 'SparsePrior.mat'];

%% Grab render matrices/files 
%
% These need to have been precomputed, and will have been if our code logic
% is correct at the calling level.
%
% Grab foward cone mosaic and render matrix
if (~exist(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'file'))
    error('Forward render strucure not cached')
else
    clear forwardRenderStructure;
    load(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'renderStructure');
    forwardRenderStructure = renderStructure; clear renderStructure; 
    grabRenderStruct(forwardRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
        pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters);
end

% Grab recon cone mosaic and render matrix
if (~exist(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'file'))
    error('Recon render strucure not cached');
else
    clear reconRenderStructure;
    load(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'renderStructure');
    reconRenderStructure = renderStructure; clear renderStructure; 
    grabRenderStruct(reconRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
        pr.nPixels, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconDefocusDiopters);
end

% Set forward variables from loaded/built structure
forwardRenderMatrix = forwardRenderStructure.renderMatrix;
forwardConeMosaic = forwardRenderStructure.theConeMosaic;
forwardOI = forwardConeMosaic.PSF;

% Set recon variables from loaded/built structure
reconRenderMatrix = reconRenderStructure.renderMatrix;
reconConeMosaic = reconRenderStructure.theConeMosaic;
reconOI = reconConeMosaic.PSF;

% Clear raw forward render structure
clear forwardRenderStructure;
clear reconRenderStructure;

%% Setup output directories
if (~exist(cnv.outputDir,'dir'))
    mkdir(cnv.outputDir);
end

%% Show forward and recon cone mosaics
forwardConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(cnv.outputDir,'forwardMosaic.jpg'),'jpg');

reconConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(cnv.outputDir,'reconMosaic.jpg'),'jpg');

%% Generate an image stimulus
%
% Stimulus size in retinal degrees should not exceed 'cnv.fieldSizeDegs'
stimSizeFraction = pr.stimSizeDegs / cnv.fieldSizeDegs;
if (stimSizeFraction > 1)
    error('Stimulus size too big given field size');
end
idxLB = round(pr.nPixels * (0.5 - stimSizeFraction / 2));
if (idxLB < 1)
    idxLB = 1;
end
idxUB = round(pr.nPixels * (0.5 + stimSizeFraction / 2));
if (idxUB > pr.nPixels)
    idxUB = pr.nPixels;
end

% Shift the stimulus to be centered on desired values
idxXRange = (idxLB:idxUB) + pr.stimCenter(1);
idxYRange = (idxLB:idxUB) + pr.stimCenter(2);

% Check stimulus position. Ends function and deletes the created
% directory if the stimulus position exceeds bounds. 
if min(idxYRange) <= 0 || max(idxYRange) > pr.nPixels ...
        || min(idxXRange) <= 0 || max(idxXRange) > pr.nPixels
    warning(['Stimulus centered on ' int2str(pr.stimCenter' + pr.trueCenter) ...
        ' exceeds bounds. Beginning next simulation']);
    rmdir(cnv.outputDir, 's');
    close all
    return
end

% Set image pixels
stimulusImageRGB = ones(pr.nPixels, pr.nPixels, 3) * pr.stimBgVal;
stimulusImageRGB(idxYRange, idxXRange, 1) = pr.stimRVal;
stimulusImageRGB(idxYRange, idxXRange, 2) = pr.stimGVal;
stimulusImageRGB(idxYRange, idxXRange, 3) = pr.stimBVal;

% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);
stimulusScene = sceneSet(stimulusScene, 'fov', cnv.fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(cnv.outputDir,'Stimulus.jpg'),'jpg');

%% Compute forward retinal image and excitations using ISETBio
%
% We may or may not reconstruct from these
forwardOI = oiCompute(stimulusScene,forwardOI);
visualizeOpticalImage(forwardOI, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(cnv.outputDir,'forwardStimulusRetinalImage.jpg'),'jpg');
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
saveas(gcf,fullfile(cnv.outputDir,'ISETBioVsRenderMatrixExciations.jpg'),'jpg');

%% Choose which excitations to reconstruct form
if (pr.reconstructfromRenderMatrix)
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusRenderMatrix;
else
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusISETBio;
end

% Visualization of the cone response note that we are using
% 'activationRange', [0 max(coneExcitations)] to avoid confusions due to
% small stimulus
if (pr.reconstructfromRenderMatrix)
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
saveas(gcf,fullfile(cnv.outputDir,'forwardMosaicExcitations.jpg'),'jpg');

%% Run reconstruction
%
% Load prior
prior = load(fullfile(pr.aoReconDir, 'priors', sparsePriorName));

% Construct image estimator
estimator = PoissonSparseEstimator(reconRenderMatrix, inv(prior.regBasis), ...
    prior.mu', pr.regPara, pr.stride, [pr.nPixels pr.nPixels 3]);

% Estimate
%
% Scale excitations to take into account difference in forward and recon
% pupil sizes.  This helps keep things in range.
%
% DHB: I think I might have scale factor backwards, now that I look at it.
% If we use a large forward pupil and reconstruct with a small one, we want
% to scale the cone exciations passed to the recon to be smaller.  So it
% should be (cnv.reconPupilDiamMM)^2/(fowardPupilDiamMM)^2.  Leaving it for
% now, but as soon as the two pupil sizes differ, we'll want to check
% carefully.
meanLuminanceCdPerM2 = [];
scaleFactor = (cnv.forwardPupilDiamMM/cnv.reconPupilDiamMM)^2;

% Set up uniform field starts
specifiedStarts = {};
for uu = 1:length(pr.uniformStartVals)
    specifiedStarts{uu} = pr.uniformStartVals(uu)*ones(length(stimulusImageLinear(:)), 1);
end

% Run the estimator
[multistartStruct,~,reconIndex] = estimator.runMultistartEstimate(forwardExcitationsToStimulusUse * scaleFactor, ...
    'maxIter', pr.maxReconIterations, 'display', 'iter', 'gpu', false, ...
    'nWhiteStart', pr.whiteNoiseStarts, 'nPinkStart', pr.pinkNoiseStarts, ...
    'nSparsePriorPatchStart', pr.sparsePriorPatchStarts, 'sparsePrior', prior, ...
    'specifiedStarts', specifiedStarts);

% Diagnose reconstructions
% txtFileName = fullfile(cnv.outputDir,'ReconProbInfo.txt');
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

    % Show stimulus, init, and recon
    theFig = figure; clf;
    set(theFig,'Position',[300 400 1150 1150]);
    theAxes = subplot(3,3,1);
    visualizeScene(stimulusScene, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    theAxes = subplot(3,3,2);
    visualizeScene(initSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    theAxes = subplot(3,3,3);
    visualizeScene(reconSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);

    % Plot predicted versus stim excitations
    subplot(3,3,5); hold on;
    minVal = 0.9*min([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)]);
    maxVal = 1.1*max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)]);
    plot(multistartStruct.coneVec,multistartStruct.reconPreds(:,ii),'ro','MarkerFaceColor','r','MarkerSize',6);
    xlabel('Measured Excitations');
    ylabel('Predicted Exciations');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    axis('square');
    title({sprintf('Recon %d, init %s',ii,multistartStruct.initTypes{ii}) ; sprintf('Iters = %d',pr.maxReconIterations) });

    % Priors, likelihoods, and losses
    subplot(3,3,7);
    bar([1]', ...
        [multistartStruct.initLogPriors(ii)  ; ...
         multistartStruct.reconLogPriors(ii) ; ...
         -stimNegLogPrior]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Log Prior');
    axis('square');
    title('Init/Recon/Stim Log Priors');
    subplot(3,3,8);
    bar([1]', ...
        [multistartStruct.initLogLikelihoods(ii)  ; ...
         multistartStruct.reconLogLikelihoods(ii) ; ...
         -stimNegLogLikely]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Log Likelihood');
    axis('square');
    title('Init/Recon/Stim Log Likelihoods');
    subplot(3,3,9);
    bar([1]', ...
        [-multistartStruct.initLosses(ii)  ; ...
         -multistartStruct.reconLosses(ii) ; ...
         -stimLoss]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Neg Loss');
    axis('square');

    % The loss figure title gives a useful summary
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
    saveas(gcf,fullfile(cnv.outputDir,sprintf('Recon%dSummary.jpg',ii)),'jpg');

    % Save summary of best recon in its own file
    if (ii == reconIndex)
        saveas(gcf,fullfile(cnv.outputDir,sprintf('ReconSummary.jpg',ii)),'jpg');
    end
end
% fclose(fid);

% Save best
[reconScene, ~, reconImageLinear] = sceneFromFile(gammaCorrection(multistartStruct.reconImages{reconIndex}, forwardConeMosaic.Display), 'rgb', ...
        meanLuminanceCdPerM2, forwardConeMosaic.Display);
reconScene = sceneSet(reconScene, 'fov', cnv.fieldSizeDegs);
visualizeScene(reconScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
saveas(gcf,fullfile(cnv.outputDir,'Recon.jpg'),'jpg');

% Compute forward excitations from reconstruction
% And compare with stimulus exciations
forwardOI = oiCompute(reconScene,forwardOI);
if (pr.reconstructfromRenderMatrix)
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
saveas(gcf,fullfile(cnv.outputDir,'StimulusVsReconExcitations.jpg'),'jpg');

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
save(fullfile(cnv.outputDir,'xRunOutput.mat'));
