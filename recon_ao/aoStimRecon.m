function aoStimRecon(pr, cnv)
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

%% Close existing figures and set prior name
close all;

% Sparse prior name
sparsePriorName = [pr.sparsePriorStr 'SparsePrior.mat'];

%% Grab render matrices/files
%
% These need to have been precomputed, and will have been if our code logic
% is correct at the calling level.
%
% Grab foward cone mosaic and render matrix
if (~exist(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'file'))
    error('Forward render strucure not cached')
else
    clear forwardRenderStructure;
    load(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'renderStructure');
    forwardRenderStructure = renderStructure; clear renderStructure;
    grabRenderStruct(forwardRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
        pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters);
end

% Set forward variables from loaded/built structure. Scale matrix with display factor.
forwardRenderMatrix = forwardRenderStructure.renderMatrix*pr.displayScaleFactor;
forwardConeMosaic = forwardRenderStructure.theConeMosaic;
forwardOI = forwardConeMosaic.PSF;
clear forwardRenderStructure;

% Grab recon cone mosaic and render matrix
if (~exist(fullfile(cnv.reconRenderDirFull, cnv.renderName),'file'))
    error('Recon render strucure not cached');
else
    clear reconRenderStructure;
    load(fullfile(cnv.reconRenderDirFull, cnv.renderName),'renderStructure');
    reconRenderStructure = renderStructure; clear renderStructure;
    grabRenderStruct(reconRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
        pr.nPixels, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconDefocusDiopters);
end

% Set recon variables from loaded/built structure. Scale matrix with display factor.
reconRenderMatrix = reconRenderStructure.renderMatrix*pr.displayScaleFactor;
reconConeMosaic = reconRenderStructure.theConeMosaic;
reconOI = reconConeMosaic.PSF;
clear reconRenderStructure;

% Scale display to match scaling of render matrices above.
forwardConeMosaic.Display = displaySet(forwardConeMosaic.Display,'spd primaries',displayGet(forwardConeMosaic.Display,'spd primaries')*pr.displayScaleFactor);
forwardConeMosaic.Display.ambient = displayGet(forwardConeMosaic.Display,'black spd')*pr.displayScaleFactor;
reconConeMosaic.Display = displaySet(reconConeMosaic.Display,'spd primaries',displayGet(reconConeMosaic.Display,'spd primaries')*pr.displayScaleFactor);
reconConeMosaic.Display.ambient = displayGet(reconConeMosaic.Display,'black spd')*pr.displayScaleFactor;

% Current sanity check to ensure that cone type at a given position does
% not change between forward/recon. At some point may want to make this an
% optional flag as we explore the idea of error in cone type assignment
if ~isequal(forwardConeMosaic.Mosaic.lConeIndices, reconConeMosaic.Mosaic.lConeIndices) || ...
        ~isequal(forwardConeMosaic.Mosaic.mConeIndices, reconConeMosaic.Mosaic.mConeIndices) || ...
        ~isequal(forwardConeMosaic.Mosaic.sConeIndices, reconConeMosaic.Mosaic.sConeIndices)
    error(['Cone types in forward and recon mosaics do not match']);
end

%% Setup output directories if they don't exist
if (~exist(cnv.outputDirFull,'dir'))
    mkdir(cnv.outputDirFull);
end

%% Show and save forward and recon cone mosaics
forwardConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(cnv.outputDirFull,'forwardMosaic.tiff'),'tiff');

% Save a version of the image with the stimulus region superimposed
forwardConeMosaic.Mosaic.visualize();
scaleForOI = 1; %1.23;

% Visualize the center used for stimulus presentation
hold on; 
plot(pr.eccXDegs, pr.eccYDegs, '.', 'MarkerSize', 40);

% Create the box boundaries
xBounds = []; yBounds = [];
xBounds(1,:) = [pr.eccXDegs pr.eccXDegs];
yBounds(1,:) = [pr.eccYDegs pr.eccYDegs];
centerWidth = scaleForOI * (pr.stimSizeDegs/2);
xBounds(2,:) = [xBounds(1)-centerWidth, xBounds(2)+centerWidth];
yBounds(2,:) = [yBounds(1)-centerWidth, yBounds(2)+centerWidth];

% Outline the stimulus region
hold on
rectangle('Position', [xBounds(2,1) yBounds(2,1) ...
    (xBounds(2,2) - xBounds(2,1)) ...
    (yBounds(2,2) - yBounds(2,1))], 'LineWidth', 3)
saveas(gcf,fullfile(cnv.outputDirFull, 'forwardMosaicwithBox.tiff'), 'tiff');

% Show and save recon mosaic
reconConeMosaic.visualizeMosaic();
saveas(gcf,fullfile(cnv.outputDirFull,'reconMosaic.tiff'),'tiff');

%% Generate an image stimulus
%
% When pr.stimBgVal is a scalar or three vector, we construct a uniform field of
% appropriate size and then put the stimulus on top of it.
if (length(pr.stimBgVal) == 1 || length(pr.stimBgVal) == 3)
    % Set up the uniform field
    if (length(pr.stimBgVal) == 1)
        stimImagergb = ones(pr.nPixels, pr.nPixels, 3) * pr.stimBgVal;
    elseif length(pr.stimBgVal) == 3
        stimImagergb(:,:,1) = ones(pr.nPixels, pr.nPixels) * pr.stimBgVal(1);
        stimImagergb(:,:,2) = ones(pr.nPixels, pr.nPixels) * pr.stimBgVal(2);
        stimImagergb(:,:,3) = ones(pr.nPixels, pr.nPixels) * pr.stimBgVal(3);
    end

    % Create stimulus to populate the invidivual quadrants
    % Stimulus size in retinal degrees should not exceed 'cnv.fieldSizeDegs'
    stimSizeFraction = pr.stimSizeDegs / cnv.fieldSizeDegs;
    if (stimSizeFraction > 1)
        error('Stimulus size too big given field size');
    end

    % Use stimulus size in pixels passed in.  This is computed
    % to deal with pixel quantization in the calling routine.
    idxLB = (pr.trueCenter)-(pr.stimSizePixels-1)/2;
    idxUB = (pr.trueCenter)+(pr.stimSizePixels-1)/2;
    idxXRange = (idxLB:idxUB);
    idxYRange = (idxLB:idxUB);
    % Could add back in the pixel shift here

    % OLD STIM PLACEMENT CODE.
    % idxLB = round(pr.nPixels * (0.5 - stimSizeFraction / 2));
    % if (idxLB < 1)
    %     idxLB = 1;
    % end
    % idxUB = round(pr.nPixels * (0.5 + stimSizeFraction / 2));
    % if (idxUB > pr.nPixels)
    %     idxUB = pr.nPixels;
    % end
    %
    % Shift the stimulus to be centered on desired values
    % idxXRange = (idxLB:idxUB) + pr.stimCenter(1);
    % idxYRange = (idxLB:idxUB) + pr.stimCenter(2);

    % Check stimulus position. Ends function and deletes the created
    % directory if the stimulus position exceeds bounds.
    if min(idxYRange) <= 0 || max(idxYRange) > pr.nPixels ...
            || min(idxXRange) <= 0 || max(idxXRange) > pr.nPixels
        % warning(['Stimulus centered on ' int2str(pr.stimCenter' + pr.trueCenter) ...
        %     ' exceeds bounds. Beginning next simulation']);
        rmdir(cnv.outputDirFull, 's');
        close all
        error('Stimulus out of image bounds');
    end

    % Set image pixels, Here's a place to fix naming convention in the
    % pr structure.
    stimImagergb(idxYRange, idxXRange, 1) = pr.stimrVal;
    stimImagergb(idxYRange, idxXRange, 2) = pr.stimgVal;
    stimImagergb(idxYRange, idxXRange, 3) = pr.stimbVal;


    % Otherwise, treat passed pr.stimBgVal as an actual image
else
    stimImagergb = pr.stimBgVal;
    nPixelsCheck = size(stimImagergb,1);
    if (nPixelsCheck ~= pr.nPixels)
        error('Passed image does not have correct pixel dimension');
    end
end

% Leave mean luminance alone
meanLuminanceCdPerM2 = [];

% Then we scale by some factor on the linear space and gamma correct
stimulusImageRGB = gammaCorrection(stimImagergb*pr.inputImageScaleFactor, forwardConeMosaic.Display);

% Then we UNDO that gamma correction we just did to get another linear
% verison, which in theory should be the same as the original linear image
% times the input scale factor.
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, forwardConeMosaic.Display);

% Set the field size
stimulusScene = sceneSet(stimulusScene, 'fov', cnv.fieldSizeDegs);

% Save the stimulus image as entered
imwrite(stimulusImageRGB,fullfile(cnv.outputDirFull,'Stimulus.tiff'),'tiff');

% Use the compareRenderingEW file to render the stimulus across the proper
% display and collect pertinent wavelength information
stimInfo = compareRenderingEW(stimulusImageLinear, pr.displayName, ...
    pr.viewingDisplayName, idxXRange, 'wls', pr.wls);

imwrite(stimInfo.imageRGBAcrossDisplays, ...
    fullfile(cnv.outputDirFull,'StimulusDispCorrected.tiff'),'tiff');

%% Compute forward retinal image and excitations using ISETBio
%
% We may or may not reconstruct from these
forwardOI = oiCompute(stimulusScene,forwardOI);
forwardExcitationsToStimulusISETBio = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));

% Check forward exciations calculation another way.  Also shows another way
% to visualize the retinal image.
forwardExcitationsToStimulusCheck = forwardConeMosaic.compute(stimulusImageRGB);
if (max(abs(forwardExcitationsToStimulusCheck-forwardExcitationsToStimulusISETBio)) ~= 0)
    forwardConeMosaic.visualizeOI()
    error('Two ways of doing the same thing do not agree');
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
title('Mean excitations ISETBio and render matrix');
saveas(gcf,fullfile(cnv.outputDirFull,'ISETBioVsRenderMatrixExciations.tiff'),'tiff');

%% Choose which excitations to reconstruct form
if (pr.reconstructfromRenderMatrix)
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusRenderMatrix;
else
    forwardExcitationsToStimulusUse = forwardExcitationsToStimulusISETBio;
end

%% Add noise?
if (pr.addPoissonNoise)
    forwardExcitationsToStimulusUse = iePoisson(forwardExcitationsToStimulusUse);
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
    'outlinedConesWithIndices', [], ...
    'activation', reshape(forwardExcitationsToStimulusUse,1,1,length(forwardExcitationsToStimulusUse)), ...
    'activationRange', [0 max(forwardExcitationsToStimulusUse)], ...
    'plotTitle',  titleStr,'labelConesInActivationMap', false);
saveas(gcf,fullfile(cnv.outputDirFull,'forwardMosaicExcitations.tiff'),'tiff');
forwardConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', axesHandle, ...
    'outlinedConesWithIndices', [], ...
    'activation', reshape(forwardExcitationsToStimulusUse,1,1,length(forwardExcitationsToStimulusUse)), ...
    'activationRange', [0 max(forwardExcitationsToStimulusUse)], ...
    'plotTitle',  titleStr,'labelConesInActivationMap', true);
saveas(gcf,fullfile(cnv.outputDirFull,'forwardMosaicExcitationsTypes.tiff'),'tiff');

%% Run reconstruction
%
% Load prior
prior = load(fullfile(pr.aoReconDir, 'priors', sparsePriorName));

% Scale factor to take into account difference in forward and recon
% sizes. We adjust the recon render matrix by this factor so
% that pupil size doesn't affect reconstruction purely by the
% scalar change in excitations.
meanLuminanceCdPerM2 = [];
pupilSizeScaleFactor = (cnv.reconPupilDiamMM/cnv.forwardPupilDiamMM)^2;
reconRenderMatrixPupilScaled = reconRenderMatrix/pupilSizeScaleFactor;
clear reconRenderMatrix;

% Construct image estimator
estimator = PoissonSparseEstimator(reconRenderMatrixPupilScaled, inv(prior.regBasis), ...
    prior.mu', pr.regPara, pr.stride, [pr.nPixels pr.nPixels 3]);

% Set up uniform field starts
specifiedStarts = {};
for uu = 1:length(pr.uniformStartVals)
    clear temp
    temp(:,:,1) = pr.uniformStartVals(1,uu)* ones(pr.nPixels,pr.nPixels);
    temp(:,:,2) = pr.uniformStartVals(2,uu)* ones(pr.nPixels,pr.nPixels);
    temp(:,:,3) = pr.uniformStartVals(3,uu)* ones(pr.nPixels,pr.nPixels);
    specifiedStarts{uu} = temp(:);
end

% Start from stimulus?
if (pr.stimulusStart)
    specifiedStarts{length(specifiedStarts)+1} = stimulusImageLinear(:);
end

% Run the estimator
if (pr.boundedSearch)
    ub = 1;
else
    ub = 100;
end
[multistartStruct,~,reconIndex] = estimator.runMultistartEstimate(forwardExcitationsToStimulusUse, ...
    'maxIter', pr.maxReconIterations, 'display', 'iter', 'gpu', false, ...
    'nWhiteStart', pr.whiteNoiseStarts, 'nPinkStart', pr.pinkNoiseStarts, ...
    'nSparsePriorPatchStart', pr.sparsePriorPatchStarts, 'sparsePrior', prior, ...
    'specifiedStarts', specifiedStarts, ...
    'ub', ub);

% Evaluate stimulus
[stimNegLogPrior,~,stimNegLogLikely] = ...
    estimator.evalEstimate(forwardExcitationsToStimulusUse, stimulusImageLinear(:));
stimLoss = stimNegLogPrior + stimNegLogLikely;

% Get information we need to render scenes from their spectra through
% the display.
theXYZStruct = load('T_xyz1931');
wls = oiGet(forwardOI,'wave');
T_xyz = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
M_rgbToxyz = T_xyz*displayGet(forwardConeMosaic.Display,'spd primaries')*(wls(2)-wls(1));
% M_rgbToxyz1 = displayGet(forwardConeMosaic.Display,'rgb2xyz');
M_xyzTorgb = inv(M_rgbToxyz);



% % Let's make sure we understand how to get the scene.  This all seems
% % to check out, and is now commented out.
%
% stimulusSceneEnergy = sceneGet(stimulusScene,'energy');
% [stimulusSceneEnergyCalFormat,m,n] = ImageToCalFormat(stimulusSceneEnergy);
% displayPrimaries = displayGet(forwardConeMosaic.Display,'spd primaries');
% stimulusSceneLinearrgbCalFormat = displayPrimaries\stimulusSceneEnergyCalFormat;
% figure; imshow(gammaCorrection(CalFormatToImage(stimulusSceneLinearrgbCalFormat,m,n),forwardConeMosaic.Display));
%
% % Now compute XYZ from energy spectra and use that to get RGB.
% % This also works.
%
% stimulusSceneXYZCalFormat = T_xyz*stimulusSceneEnergyCalFormat*(wls(2)-wls(1));
% stimulusSceneLinearrgb1CalFormat = M_xyzTorgb*stimulusSceneXYZCalFormat;
% figure; imshow(gammaCorrection(CalFormatToImage(stimulusSceneLinearrgb1CalFormat,m,n),forwardConeMosaic.Display));
%
% % And yet one more way, which also looks good.
% [stimulusScenexyz,m,n] = ImageToCalFormat(sceneGet(stimulusScene,'xyz'));
% stimulusScenergb = M_xyzTorgb*stimulusScenexyz;
% fprintf('Min/max of stimluus image: %0.2f, %0.2f\n',min(stimulusScenergb(:)),max(stimulusScenergb(:)));
% stimulusScenergb(stimulusScenergb < 0) = 1;
% stimulusScenergb(stimulusScenergb > 1) = 1;
% stimulusSceneRGB1 = gammaCorrection(CalFormatToImage(stimulusScenergb,m,n),forwardConeMosaic.Display);
% figure; clf; imshow(stimulusSceneRGB1);


%% Summary plot of what happened
for ii = 1:length(multistartStruct.initTypes)

    % Handle bounded search for display.  We scale down the stimulus image
    % by the scale factor needed to bring the reconstruction into range
    % 0-1.
    reconScaleFactor(ii) = max(multistartStruct.reconImages{ii}(:));
    if (reconScaleFactor(ii) < 1)
        reconScaleFactor(ii) = 1;
    end
    stimulusRGBScaled{ii} = gammaCorrection(stimulusImageLinear/reconScaleFactor(ii), reconConeMosaic.Display);
    maxStimulusScaledR(ii) = max(max(stimulusRGBScaled{ii}(:,:,1)));
    maxStimulusScaledG(ii) = max(max(stimulusRGBScaled{ii}(:,:,2)));
    maxStimulusScaledB(ii) = max(max(stimulusRGBScaled{ii}(:,:,3)));

    % Set up initial scene, that is scene from the initialization.  We just want a look at this,
    % so we don't fuss much with scaling.
    [initSceneTemp, ~, initImageLinearTemp] = sceneFromFile(gammaCorrection(multistartStruct.initImages{ii}, forwardConeMosaic.Display), 'rgb', ...
        meanLuminanceCdPerM2, forwardConeMosaic.Display);

    % Get the reconstruction as RGB image and find maxima.  This is scaled
    % down as needed so that the recon RGB is in range 0-1.
    reconScaledRGB{ii} = gammaCorrection(multistartStruct.reconImages{ii}/reconScaleFactor(ii), reconConeMosaic.Display);
    maxReconScaledR(ii) = max(max(reconScaledRGB{ii}(:,:,1)));
    maxReconScaledG(ii) = max(max(reconScaledRGB{ii}(:,:,2)));
    maxReconScaledB(ii) = max(max(reconScaledRGB{ii}(:,:,3)));

    % Adjust recon scene back up for recon scale factor and adjust for pupil size
    % scale factor. This puts into into the same scaling as the stimulus
    % scene, adjusted for pupil diameter.
    [reconSceneTemp, ~, reconImageLinearTemp] = sceneFromFile(reconScaledRGB{ii}, 'rgb', ...
        meanLuminanceCdPerM2, reconConeMosaic.Display);
    reconSceneTemp = sceneSet(reconSceneTemp,'photons',sceneGet(reconSceneTemp,'photons')*reconScaleFactor(ii)/pupilSizeScaleFactor);
    reconSceneTemp = sceneSet(reconSceneTemp, 'fov', cnv.fieldSizeDegs);
    reconImageLinearTemp = reconImageLinearTemp*reconScaleFactor(ii);


    % Visualize recon after being corrected for Display
    reconImageLinear = multistartStruct.reconImages{ii};

    reconInfo = compareRenderingEW(reconImageLinear, pr.displayName, ...
        pr.viewingDisplayName, idxXRange, 'wls', pr.wls);

    imwrite(reconInfo.imageRGBAcrossDisplays, fullfile(cnv.outputDirFull,'ReconDispCorrected.tiff'),'tiff');


    % Include portion for scaling of the corrected stim and recon images
    % based on the max across both, upscaled so the max is reset to 1


    %
    %     cfvStim.scaleFactor(ii) = max([max(cfvStim.imageRGBNoGamma(:)) max(cfvRecon.imageRGBNoGamma(:))]);
    %     cfvRecon.scaleFactor(ii) = max([max(cfvStim.imageRGBNoGamma(:)) max(cfvRecon.imageRGBNoGamma(:))]);
    %     cfvStim.stimulusRGBScaled{ii} = gammaCorrection(cfvStim.imageRGBNoGamma/cfvStim.scaleFactor(ii), cfvStim.viewingDisplay);
    %     cfvRecon.reconScaledRGB{ii} = gammaCorrection(cfvRecon.imageRGBNoGamma/cfvRecon.scaleFactor(ii), cfvRecon.viewingDisplay);


    % Get forward and reconstruction OI's computed on reconstruction.  Take
    % difference in pupil size into account with reconOI.
    forwardOI = oiCompute(stimulusScene,forwardOI);
    forwardOIToReconTemp = oiCompute(reconSceneTemp,forwardOI);
    reconOIToReconTemp = oiCompute(reconSceneTemp,reconOI);

    % Get recon excitations to stimulus
    reconExcitationsToStimulusTemp = reconRenderMatrixPupilScaled*stimulusImageLinear(:);

    % Render OIs on display as best we can.  Getting the scale factor right
    % from first principles is hard because it depends on the radiance ->
    % irradiance conversion and we don't want to unpack that here.  So we
    % scale both OI's we want to visualize to the same range.  We also
    % assume they are both the same size.
    [forwardOIxyz,m,n] = ImageToCalFormat(oiGet(forwardOI,'xyz'));
    forwardOIrgb = M_xyzTorgb*forwardOIxyz;
    forwardOITitleStr = {sprintf('Min/max (arb units) of forward OI image:  %0.4f, %0.4f\n',min(forwardOIrgb(:)),max(forwardOIrgb(:)))};
    [reconOIxyz,m,n] = ImageToCalFormat(oiGet(reconOIToReconTemp,'xyz'));
    reconOIrgb = M_xyzTorgb*reconOIxyz;
    reconOITitleStr = {sprintf('Min/max (arb units) of recon OI image: %0.4f, %0.3f\n',min(reconOIrgb(:)),max(reconOIrgb(:)))};
    oiScaleFactor = max([forwardOIrgb(:) ; reconOIrgb(:)]);
    forwardOIrgb = forwardOIrgb/oiScaleFactor;
    forwardOIrgb(forwardOIrgb < 0) = 0;
    forwardOIRGB = gammaCorrection(CalFormatToImage(forwardOIrgb,m,n),forwardConeMosaic.Display);
    reconOIrgb = reconOIrgb/oiScaleFactor;
    reconOIrgb(reconOIrgb < 0) = 0;
    reconOIRGB = gammaCorrection(CalFormatToImage(reconOIrgb,m,n),reconConeMosaic.Display);

    % Set up summary plot
    theFig = figure; clf;
    set(theFig,'Position',[100 400 2500 1500]);

    % Initial image
    theAxes = subplot(3,7,21);
    imshow(gammaCorrection(multistartStruct.initImages{ii}, forwardConeMosaic.Display));
    title({sprintf('Recon %d, init %s',ii,multistartStruct.initTypes{ii}) ; sprintf('Iters = %d',pr.maxReconIterations) });

    % Visualize stimulus
    theAxes = subplot(3,7,1);
    imshow(stimulusRGBScaled{ii});
    if (length(pr.stimBgVal) > 3)
        title({sprintf('Stimulus Image, input scale %0.4f',pr.inputImageScaleFactor) ; 'Scaled with recon' ; sprintf('Max scaled (image) RGB: %0.4f, %0.4f, %0.4f',maxStimulusScaledR(ii),maxStimulusScaledG(ii),maxStimulusScaledB(ii)) ; pr.imageName});
    else
        title({sprintf('Stimulus Image, input scale %0.4f',pr.inputImageScaleFactor)  ; 'Scaled with recon' ; sprintf('Max scaled (image) RGB: %0.4f, %0.4f, %0.4f',maxStimulusScaledR(ii),maxStimulusScaledG(ii),maxStimulusScaledB(ii)) ; sprintf('%0.4f, %0.4f, %0.4f, %0.4f',pr.stimBgVal(1),pr.stimrVal,pr.stimgVal,pr.stimbVal)});
    end
    %     if (ii == reconIndex)
    %         imwrite(stimulusRGBScaled{ii},fullfile(cnv.outputDirFull,'StimulusScaled.tiff'),'tiff');
    %     end

    % Visualize stimulus after being corrected for Display
    theAxes = subplot(3,7,7);
    imshow(stimInfo.imageRGBAcrossDisplays);
    %     title({sprintf('Stim on %s, viewed on %s', pr.displayName, pr.viewingDisplayName); ...
    %         sprintf('Min: %0.2f, %0.2f, %0.2f',cfvStim.bounds(1,1),cfvStim.bounds(1,2),cfvStim.bounds(1,3)); ...
    %         sprintf('Max: %0.2f, %0.2f, %0.2f',cfvStim.bounds(2,1),cfvStim.bounds(2,2),cfvStim.bounds(2,3))})


    % Contour plot of forward PSF
    theAxes = subplot(3,7,2);
    cmap = brewermap(1024,'blues');
    alpha = 0.75;
    contourLineColor = [0.2 0.2 0.2];
    zLevels = 0.05:0.15:0.95;
    psfSupportTemp = opticsGet(oiGet(forwardOI,'optics'),'psf support');
    psfPolyTemp = opticsGet(oiGet(forwardOI,'optics'),'psf data');
    psfTemp = zeros(size(squeeze(psfPolyTemp(:,:,1))));
    for ww = 1:size(T_xyz,2)
        psfTemp = psfTemp+T_xyz(2,ww)*squeeze(psfPolyTemp(:,:,ww));
    end
    psfDataStruct = struct(...
        'supportXdegs', psfSupportTemp{1}(1,:), ...  % [1 x xPoints] spatial support vector (x)
        'supportYdegs', psfSupportTemp{2}(:,1), ...  % [mPoints x 1] spatial support vector (y)
        'data', psfTemp ...        % [mPoints x nPoints] 2D PSF
        );
    cMosaic.semiTransparentContourPlot(theAxes, ...
        psfDataStruct.supportXdegs,...
        psfDataStruct.supportYdegs, ...
        psfDataStruct.data/max(psfDataStruct.data(:)), ...
        zLevels, cmap, alpha, contourLineColor, ...
        'lineWidth', 1.5);
    % xlim([min(psfSupportTemp{1}(1,:))   max(psfSupportTemp{1}(1,:))]);
    % ylim([min(psfSupportTemp{2}(:,1)) ; max(psfSupportTemp{2}(:,1))]);
    xlim([-0.05 0.05]); ylim([-0.05 0.05]); axis('square');
    title({sprintf('%s',cnv.forwardAOStrPlt); sprintf('%0.2fDD, %s', pr.forwardDefocusDiopters, cnv.forwardLCAStr); 'Forward lum weighted PSF'});
    xlabel('X (degs)'); ylabel('Y (degs)');

    % Optical image of stimulus
    theAxes = subplot(3,7,3);
    % visualizeOpticalImage(forwardOI, 'axesHandle',theAxes,'avoidAutomaticRGBscaling', true, 'displayRadianceMaps', false);
    imshow(forwardOIRGB);
    title(forwardOITitleStr);

    % Show forward mosaic
    theAxes = subplot(3,7,4);
    figureHandle = theFig;
    forwardConeMosaic.visualizeMosaic(figureHandle,theAxes);
    title('Forward Mosaic');

    % Forward excitations used for recon in mosaic form
    drawnow
    theAxes = subplot(3,7,5);
    figureHandle = theFig;
    forwardConeMosaic.Mosaic.visualize(...
        'figureHandle', figureHandle, ...
        'axesHandle', theAxes, ...
        'withSuperimposedOpticalImage', forwardOI, ...
        'outlinedConesWithIndices', [], ...
        'plotTitle','Forward OI on Forward Mosaic','superimposedOIAlpha',0.7);
    if (ii == reconIndex)
        tempFig = figure; clf;
        forwardConeMosaic.Mosaic.visualize(...
            'figureHandle', tempFig, ...
            'axesHandle', [], ...
            'withSuperimposedOpticalImage', forwardOI, ...
            'outlinedConesWithIndices', [], ...
            'plotTitle','Forward OI on Forward Mosaic','superimposedOIAlpha',0.7);
        saveas(tempFig,fullfile(cnv.outputDirFull,sprintf('forwardOIOnForwardMosaic.tiff',ii)),'tiff');
        close(tempFig);
        figure(theFig);
    end

    theAxes = subplot(3,7,6);
    drawnow
    figureHandle = theFig;
    forwardConeMosaic.Mosaic.visualize(...
        'figureHandle', figureHandle, ...
        'axesHandle', theAxes, ...
        'outlinedConesWithIndices', [], ...
        'activation', reshape(multistartStruct.coneVec,1,1,length(forwardExcitationsToStimulusUse)), ...
        'activationRange', 1.1*[0 max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)])], ...
        'plotTitle',  'Forward excitations','labelConesInActivationMap', false);

    theAxes = subplot(3,7,8);
    drawnow
    %visualizeScene(reconSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    imshow(reconScaledRGB{ii});
    if (pr.boundedSearch)
        title({sprintf('Reconstructed Image, reg %0.5f',pr.regPara) ; sprintf('Max scaled (image) RGB: %0.4f, %0.4f, %0.4f',maxReconScaledR(ii),maxReconScaledG(ii),maxReconScaledB(ii)) ; 'Bounded search' ; sprintf('Recon scale factor %0.3g',reconScaleFactor(ii))});
    else
        title({sprintf('Reconstructed Image, reg %0.5f',pr.regPara) ; sprintf('Max scaled (image) RGB: %0.4f, %0.4f, %0.4f',maxReconScaledR(ii),maxReconScaledG(ii),maxReconScaledB(ii)) ; 'Unbounded search' ; sprintf('Recon scale factor %0.3g',reconScaleFactor(ii))});
    end


    % Show corrected recon image
    theAxes = subplot(3,7,14);
    imshow(reconInfo.imageRGBAcrossDisplays);
    %     title({sprintf('Recon on %s, viewed on %s', pr.displayName, pr.viewingDisplayName); ...
    %         sprintf('Min: %0.2f, %0.2f, %0.2f',cfvRecon.bounds(1,1),cfvRecon.bounds(1,2),cfvRecon.bounds(1,3)); ...
    %         sprintf('Max: %0.2f, %0.2f, %0.2f',cfvRecon.bounds(2,1),cfvRecon.bounds(2,2),cfvRecon.bounds(2,3))})


    % Contour plot of recon PSF
    theAxes = subplot(3,7,9);
    cmap = brewermap(1024,'blues');
    alpha = 0.75;
    contourLineColor = [0.2 0.2 0.2];
    zLevels = 0.05:0.15:0.95;
    psfSupportTemp = opticsGet(oiGet(reconOI,'optics'),'psf support');
    psfPolyTemp = opticsGet(oiGet(reconOI,'optics'),'psf data');
    psfTemp = zeros(size(squeeze(psfPolyTemp(:,:,1))));
    for ww = 1:size(T_xyz,2)
        psfTemp = psfTemp+T_xyz(2,ww)*squeeze(psfPolyTemp(:,:,ww));
    end
    psfDataStruct = struct(...
        'supportXdegs', psfSupportTemp{1}(1,:), ...  % [1 x xPoints] spatial support vector (x)
        'supportYdegs', psfSupportTemp{2}(:,1), ...  % [mPoints x 1] spatial support vector (y)
        'data', psfTemp ...        % [mPoints x nPoints] 2D PSF
        );
    cMosaic.semiTransparentContourPlot(theAxes, ...
        psfDataStruct.supportXdegs,...
        psfDataStruct.supportYdegs, ...
        psfDataStruct.data/max(psfDataStruct.data(:)), ...
        zLevels, cmap, alpha, contourLineColor, ...
        'lineWidth', 1.5);
    % xlim([min(psfSupportTemp{1}(1,:))   max(psfSupportTemp{1}(1,:))]);
    % ylim([min(psfSupportTemp{2}(:,1)) ; max(psfSupportTemp{2}(:,1))]);
    xlim([-0.05 0.05]); ylim([-0.05 0.05]); axis('square');
    title({sprintf('%s', cnv.reconAOStrPlt); sprintf('%0.2fDD, %s', pr.reconDefocusDiopters, cnv.reconLCAStr); 'Recon lum weighted PSF'});
    xlabel('X (degs)'); ylabel('Y (degs)');

    % Optical image of recon
    theAxes = subplot(3,7,10);
    % visualizeOpticalImage(reconOIToReconTemp, 'axesHandle',theAxes,'avoidAutomaticRGBscaling', true, 'displayRadianceMaps', false);
    imshow(reconOIRGB);
    title(reconOITitleStr);

    % Show recon mosaic
    theAxes = subplot(3,7,11);
    figureHandle = theFig;
    reconConeMosaic.visualizeMosaic(figureHandle,theAxes);
    title('Recon Mosaic');

    % Recon oi on recon mosaic
    drawnow
    theAxes = subplot(3,7,12);
    figureHandle = theFig;
    reconConeMosaic.Mosaic.visualize(...
        'figureHandle', figureHandle, ...
        'axesHandle', theAxes, ...
        'outlinedConesWithIndices', [], ...
        'withSuperimposedOpticalImage', reconOIToReconTemp, ...
        'plotTitle','Recon OI on Recon Mosaic','superimposedOIAlpha',0.7);
    if (ii == reconIndex)
        tempFig = figure; clf;
        reconConeMosaic.Mosaic.visualize(...
            'figureHandle', tempFig, ...
            'axesHandle', [], ...
            'outlinedConesWithIndices', [], ...
            'withSuperimposedOpticalImage', reconOIToReconTemp, ...
            'plotTitle','Recon OI on Recon Mosaic','superimposedOIAlpha',0.7);
        saveas(tempFig,fullfile(cnv.outputDirFull,sprintf('reconOIOnReconMosaic.tiff',ii)),'tiff');
        close(tempFig);
        figure(theFig);
    end

    % Recon excitations
    drawnow
    theAxes = subplot(3,7,13);
    figureHandle = theFig;
    reconConeMosaic.Mosaic.visualize(...
        'figureHandle', figureHandle, ...
        'axesHandle', theAxes, ...
        'outlinedConesWithIndices', [], ...
        'activation', reshape(multistartStruct.reconPreds(:,ii),1,1,length(forwardExcitationsToStimulusUse)), ...
        'activationRange', 1.1*[0 max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)])], ...
        'plotTitle',  'Recon excitations','labelConesInActivationMap', false);
    if (ii == reconIndex)
        tempFig = figure; clf;
        reconConeMosaic.Mosaic.visualize(...
            'figureHandle', tempFig, ...
            'axesHandle', [], ...
            'outlinedConesWithIndices', [], ...
            'activation', reshape(multistartStruct.reconPreds(:,ii),1,1,length(forwardExcitationsToStimulusUse)), ...
            'activationRange', 1.1*[0 max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)])], ...
            'plotTitle',  'Recon excitations','labelConesInActivationMap', false);
        saveas(tempFig,fullfile(cnv.outputDirFull,sprintf('reconMosaicExcitations.tiff',ii)),'tiff');
        clf;
        reconConeMosaic.Mosaic.visualize(...
            'figureHandle', tempFig, ...
            'axesHandle', [], ...
            'outlinedConesWithIndices', [], ...
            'activation', reshape(multistartStruct.reconPreds(:,ii),1,1,length(forwardExcitationsToStimulusUse)), ...
            'activationRange', 1.1*[0 max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)])], ...
            'plotTitle',  'Recon excitations','labelConesInActivationMap', true);
        saveas(tempFig,fullfile(cnv.outputDirFull,sprintf('reconMosaicExcitationsTypes.tiff',ii)),'tiff');
        close(tempFig);
        figure(theFig);
    end

    % Make sure excitations used match what comes back from multistart
    if (any(forwardExcitationsToStimulusUse ~= multistartStruct.coneVec))
        error('Inconsistency in excitations driving reconstruction');
    end

    % Compute recon excitations to stimulus and compare with
    % scaled forward excitations to stimulus.
    drawnow
    subplot(3,7,15); hold on;
    if (pr.reconstructfromRenderMatrix)
        title({'Recon excittions to stim' ; 'Excitations from render matrix'});
        forwardExcitationsToReconTemp = forwardRenderMatrix*reconImageLinearTemp(:);
    else
        title({'Recon excittions to stim' ; 'Excitations from ISETBio'});
        forwardExcitationsToReconTemp = forwardConeMosaic.Mosaic.compute(forwardOIToReconTemp, 'opticalImagePositionDegs', 'mosaic-centered');
        forwardExcitationsToReconTemp = squeeze(forwardExcitationsToReconTemp);
    end

    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.lConeIndices),reconExcitationsToStimulusTemp(reconConeMosaic.Mosaic.lConeIndices),'ro','MarkerFaceColor','r','MarkerSize',6); hold on;
    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.mConeIndices),reconExcitationsToStimulusTemp(reconConeMosaic.Mosaic.mConeIndices),'go','MarkerFaceColor','g','MarkerSize',6); hold on;
    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.sConeIndices),reconExcitationsToStimulusTemp(reconConeMosaic.Mosaic.sConeIndices),'bo','MarkerFaceColor','b','MarkerSize',6); hold on;

    axis('square');
    minVal = 0.9*min([forwardExcitationsToStimulusUse;reconExcitationsToStimulusTemp]);
    maxVal = 1.1*max([forwardExcitationsToStimulusUse;reconExcitationsToStimulusTemp]);
    plot([minVal maxVal],[minVal maxVal],'k');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    xlabel('Forward excitations to stimulus');
    ylabel('Recon excitations to stimulus');

    % Compute forward excitations from reconstruction
    % and compare with scaled stimulus excitations
    subplot(3,7,16); hold on;
    if (pr.reconstructfromRenderMatrix)
        title({'Forward excitations to recon' ; 'Excitations from render matrix'});
        forwardExcitationsToReconTemp = forwardRenderMatrix*reconImageLinearTemp(:);
    else
        title({'Forward excitations to recon' ; 'Excitations from ISETBio'});
        forwardExcitationsToReconTemp = (forwardConeMosaic.Mosaic.compute(forwardOIToReconTemp, 'opticalImagePositionDegs', 'mosaic-centered'));
        forwardExcitationsToReconTemp = squeeze(forwardExcitationsToReconTemp);
    end
    %plot(forwardExcitationsToStimulusUse*scaleFactor,forwardExcitationsToReconTemp,'ro','MarkerFaceColor','r','MarkerSize',6);

    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.lConeIndices),forwardExcitationsToReconTemp(forwardConeMosaic.Mosaic.lConeIndices),'ro','MarkerFaceColor','r','MarkerSize',6); hold on;
    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.mConeIndices),forwardExcitationsToReconTemp(forwardConeMosaic.Mosaic.mConeIndices),'go','MarkerFaceColor','g','MarkerSize',6); hold on;
    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.sConeIndices),forwardExcitationsToReconTemp(forwardConeMosaic.Mosaic.sConeIndices),'bo','MarkerFaceColor','b','MarkerSize',6); hold on;

    axis('square');
    minVal = 0.9*min([forwardExcitationsToStimulusUse; forwardExcitationsToReconTemp]);
    maxVal = 1.1*max([forwardExcitationsToStimulusUse; forwardExcitationsToReconTemp]);
    plot([minVal maxVal],[minVal maxVal],'k');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    xlabel('Forward excitations to stimulus');
    ylabel('Forward excitations to recon');

    % Compute recon excitations from reconstruction
    % and compare with scaled stimulus excitations
    subplot(3,7,17); hold on;
    reconExcitationsToReconCheck = reconRenderMatrixPupilScaled*reconImageLinearTemp(:);
    if (pr.reconstructfromRenderMatrix)
        title({'Recon excitations to recon' ; 'Excitations from render matrix'});
        reconExcitationsToReconTemp = reconExcitationsToReconCheck;
    else
        title({'Recon excitations to recon' ; 'Excitations from ISETBio'});
        reconExcitationsToReconTemp = (reconConeMosaic.Mosaic.compute(reconOIToReconTemp, 'opticalImagePositionDegs', 'mosaic-centered'));
        reconExcitationsToReconTemp = squeeze(reconExcitationsToReconTemp);
    end

    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.lConeIndices),reconExcitationsToReconTemp(reconConeMosaic.Mosaic.lConeIndices),'ro','MarkerFaceColor','r','MarkerSize',6); hold on;
    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.mConeIndices),reconExcitationsToReconTemp(reconConeMosaic.Mosaic.mConeIndices),'go','MarkerFaceColor','g','MarkerSize',6); hold on;
    plot(forwardExcitationsToStimulusUse(forwardConeMosaic.Mosaic.sConeIndices),reconExcitationsToReconTemp(reconConeMosaic.Mosaic.sConeIndices),'bo','MarkerFaceColor','b','MarkerSize',6); hold on;

    axis('square');
    minVal = 0.9*min([forwardExcitationsToStimulusUse; reconExcitationsToReconTemp]);
    maxVal = 1.1*max([forwardExcitationsToStimulusUse; reconExcitationsToReconTemp]);
    plot([minVal maxVal],[minVal maxVal],'k');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    xlabel('Forward excitations to stimulus');
    ylabel('Recon excitations to recon');


    % Check that we know what we are doing.  Small difference may be gamma
    % correction and inverse gamma correction between the two predictions
    if (max(abs(multistartStruct.reconPreds(:,ii)-reconExcitationsToReconCheck)./reconExcitationsToReconCheck) > 1e-3)
        figure; clf; hold on;
        plot(multistartStruct.reconPreds(:,ii),reconExcitationsToReconCheck,'ro','MarkerFaceColor','r','MarkerSize',10);
        axis('square')
        minVal = 0.9*min([multistartStruct.reconPreds(:,ii); reconExcitationsToReconCheck]);
        maxVal = 1.1*max([multistartStruct.reconPreds(:,ii); reconExcitationsToReconCheck]);
        plot([minVal maxVal],[minVal maxVal],'k');
        xlim([minVal maxVal]); ylim([minVal maxVal]);
        title('Excitations to recon in two ways');
        xlabel('Excitations from multistart struct');
        ylabel('Excitations from aoStimRecon')
        saveas(gcf,fullfile(cnv.outputDirFull,sprintf('ReconExcitationsCheckError%d.tiff',ii)),'tiff');
        figure(theFig);
    end

    % Priors, likelihoods, and losses
    subplot(3,7,18);
    bar([1]', ...
        [multistartStruct.initLogPriors(ii)  ; ...
        multistartStruct.reconLogPriors(ii) ; ...
        -stimNegLogPrior]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Log Prior');
    axis('square');

    % Prior summary
    if (multistartStruct.reconLogPriors(ii) > -stimNegLogPrior)
        if (multistartStruct.reconLogPriors(ii) > multistartStruct.initLogPriors(ii))
            title({'Init/Recon/Stim Log Priors' ; 'Recon BETTER than Stim' ; 'Recon BETTER than Init'});
        else
            title({'Init/Recon/Stim Log Priors' ; 'Recon BETTER than Stim' ; 'Recon WORSE than Init'});
        end
    else
        if (multistartStruct.reconLogPriors(ii) > multistartStruct.initLogPriors(ii))
            title({'Init/Recon/Stim Log Priors' ; 'Recon WORSE than Stim' ; 'Recon BETTER than Init'});
        else
            title({'Init/Recon/Stim Log Priors' ; 'Recon WORSE than Stim' ; 'Recon WORSE than Init'});
        end
    end

    % Likelihoods
    subplot(3,7,19);
    bar([1]', ...
        [multistartStruct.initLogLikelihoods(ii)  ; ...
        multistartStruct.reconLogLikelihoods(ii) ; ...
        -stimNegLogLikely]');
    set(gca,'XTickLabel',sprintf('Recon %d',ii))
    ylabel('Log Likelihood');
    axis('square');

    % Likelihood summary
    if (multistartStruct.reconLogLikelihoods(ii) > -stimNegLogLikely)
        if (multistartStruct.reconLogLikelihoods(ii) > multistartStruct.initLogLikelihoods(ii))
            title({'Init/Recon/Stim Log Likelihoods' ; 'Recon BETTER than Stim' ; 'Recon BETTER than Init'});
        else
            title({'Init/Recon/Stim Log Likelihoods' ; 'Recon BETTER than Stim' ; 'Recon WORSE than Init'});
        end
    else
        if (multistartStruct.reconLogLikelihoods(ii) > multistartStruct.initLogLikelihoods(ii))
            title({'Init/Recon/Stim Log Likelihoods' ; 'Recon WORSE than Stim' ; 'Recon BETTER than Init'});
        else
            title({'Init/Recon/Stim Log Likelihoods' ; 'Recon WORSE than Stim' ; 'Recon WORSE than Init'});
        end
    end

    % Loss value (negative, so plus is good)
    subplot(3,7,20);
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
    saveas(gcf,fullfile(cnv.outputDirFull,sprintf('Recon%dSummaryUpdate.tiff',ii)),'tiff');

    % Save summary of best recon in its own file
    if (ii == reconIndex)
        saveas(gcf,fullfile(cnv.outputDirFull,sprintf('ReconSummaryUpdate.tiff',ii)),'tiff');
    end
end

% Save best reconstruction image
imwrite(reconScaledRGB{reconIndex},fullfile(cnv.outputDirFull,'Recon.tiff'),'tiff');

%% Save workspace without really big variables
close all;
clear forwardRenderMatrix reconRenderMatrixPupilScaled reconSceneTemp forwardOI reconOIToReconTemp psfDataStruct forwardOIToReconTemp forwardOIRGB
clear estimator
clear reconScaledRGB stimulusRGBScaled reconOI psfTemp psfPolyTemp
clear reconImageLinearTemp psfSupportTemp initImageLinearTemp
clear tempFig theAxes theFig axesHandle temp initSceneTemp
% No longer clearing the coneMosaics for analysis but be aware might lead
% to a crash of storage. Also might not be necessary as of now. 
clear forwardConeMosaic reconConeMosaic
save(fullfile(cnv.outputDirFull,'xRunOutput.mat'), '-v7.3');
end