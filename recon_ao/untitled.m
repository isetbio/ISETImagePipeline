

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

for ii = 1:length(multistartStruct.initTypes)
    % Set up initial scene.
    [initSceneTemp, ~, initImageLinearTemp] = sceneFromFile(gammaCorrection(multistartStruct.initImages{ii}, forwardConeMosaic.Display), 'rgb', ...
        meanLuminanceCdPerM2, forwardConeMosaic.Display);

    % Get the reconstruction as RGB image and find maxima
    reconRGBTemp = gammaCorrection(multistartStruct.reconImages{ii}, reconConeMosaic.Display);
    maxReconR = max(max(reconRGBTemp(:,:,1)));
    maxReconG = max(max(reconRGBTemp(:,:,2)));
    maxReconB = max(max(reconRGBTemp(:,:,3)));
    [reconSceneTemp, ~, reconImageLinearTemp] = sceneFromFile(reconRGBTemp, 'rgb', ...
        meanLuminanceCdPerM2, reconConeMosaic.Display);

    % Get forward and reconstruction OI's computed on reconstruction
    forwardOI = oiCompute(stimulusScene,forwardOI);
    forwardOIToReconTemp = oiCompute(reconSceneTemp,forwardOI);
    reconOIToReconTemp = oiCompute(reconSceneTemp,reconOI);

    % Get recon excitations to stimulus
    reconExcitationsToStimulusTemp = reconRenderMatrix*stimulusImageLinear(:);

    % Render OIs on display as best we can.  Getting the scale factor right
    % from first principles is hard because it depends on the radiance ->
    % irradiance conversion and we don't want to unpack that here.  So we
    % scale both OI's we want to visualize to the same range.  We also
    % assume they are both the same size.
    [forwardOIxyz,m,n] = ImageToCalFormat(oiGet(forwardOI,'xyz'));
    forwardOIrgb = M_xyzTorgb*forwardOIxyz*scaleFactor;
    forwardOITitleStr = {'Min/max (arb units) of scaled (pupil) forward OI image:' ; sprintf(' %0.4f, %0.4f\n',min(forwardOIrgb(:)),max(forwardOIrgb(:)))};
    [reconOIxyz,m,n] = ImageToCalFormat(oiGet(reconOIToReconTemp,'xyz'));
    reconOIrgb = M_xyzTorgb*reconOIxyz;
    reconOITitleStr = {'Min/max (arb units) of recon OI image:' ; sprintf('Min/max of recon OI image: %0.2f, %0.2f\n',min(reconOIrgb(:)),max(reconOIrgb(:)))};
    oiScaleFactor = max([forwardOIrgb(:) ; reconOIrgb(:)]);
    forwardOIrgb = forwardOIrgb/oiScaleFactor;
    forwardOIrgb(forwardOIrgb < 0) = 0;
    forwardOIRGB = gammaCorrection(CalFormatToImage(forwardOIrgb,m,n),forwardConeMosaic.Display);
    reconOIrgb = reconOIrgb/oiScaleFactor;
    reconOIrgb(reconOIrgb < 0) = 0;
    reconOIRGB = gammaCorrection(CalFormatToImage(reconOIrgb,m,n),forwardConeMosaic.Display);

    % Set up summary plot
    theFig = figure; clf;
    set(theFig,'Position',[100 400 2500 1500]);

    % Initial image
    theAxes = subplot(3,6,6);
    %visualizeScene(initSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    imshow(gammaCorrection(multistartStruct.initImages{ii}, forwardConeMosaic.Display));
    title({sprintf('Recon %d, init %s',ii,multistartStruct.initTypes{ii}) ; sprintf('Iters = %d',pr.maxReconIterations) });

    % Visualize stimulus
    theAxes = subplot(3,6,1);
    % visualizeScene(stimulusScene, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    imshow(stimulusImageRGB);
    if (length(pr.stimBgVal) > 1)
        title({'Stimulus Image' ; pr.imageName});
    else
        title({'Stimulus Image' ; sprintf('%0.4f, %0.4f, %0.4f, %0.4f',pr.stimBgVal,pr.stimRVal,pr.stimGVal,pr.stimBVal)});
    end

    % Contour plot of forward PSF
    theAxes = subplot(3,6,2);
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
    xlim([-0.05 0.05]); ylim([-0.05 0.05]);
    title('Forward lum weighted PSF');
    xlabel('X (degs)'); ylabel('Y (degs)');

    % Optical image of stimulus
    theAxes = subplot(3,6,3);
    % visualizeOpticalImage(forwardOI, 'axesHandle',theAxes,'avoidAutomaticRGBscaling', true, 'displayRadianceMaps', false);
    imshow(forwardOIRGB);
    title(forwardOITitleStr);

    % Show forward mosaic
    theAxes = subplot(3,6,4);
    figureHandle = theFig; 
    forwardConeMosaic.visualizeMosaic(figureHandle,theAxes);
    title('Forward Mosaic');

    % Forward excitations used for recon in mosaic form
    theAxes = subplot(3,6,5);
    figureHandle = theFig; 
    forwardConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', theAxes, ...
    'activation', reshape(multistartStruct.coneVec,1,1,length(forwardExcitationsToStimulusUse)), ...
    'activationRange', 1.1*[0 max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)])], ...
    'plotTitle',  'Scaled (pupil) forward excitations');

    theAxes = subplot(3,6,7);
    %visualizeScene(reconSceneTemp, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true,'axesHandle',theAxes);
    imshow(gammaCorrection(multistartStruct.reconImages{ii}, forwardConeMosaic.Display));
    title({'Reconstructed Image' ; sprintf('MaxRGB: %0.4f, %0.4f, %0.4f',maxReconR,maxReconG,maxReconB)});

    % Contour plot of recon PSF
    theAxes = subplot(3,6,8);
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
    xlim([-0.05 0.05]); ylim([-0.05 0.05]);
    title('Recon lum weighted PSF');
    xlabel('X (degs)'); ylabel('Y (degs)');
    
    % Optical image of recon
    theAxes = subplot(3,6,9);
    % visualizeOpticalImage(reconOIToReconTemp, 'axesHandle',theAxes,'avoidAutomaticRGBscaling', true, 'displayRadianceMaps', false);
    imshow(reconOIRGB);
    title(reconOITitleStr);

    % Show recon mosaic
    theAxes = subplot(3,6,10);
    figureHandle = theFig; 
    reconConeMosaic.visualizeMosaic(figureHandle,theAxes);
    title('Recon Mosaic');

    % Recon excitations to recon in mosaic form
    theAxes = subplot(3,6,11);
    figureHandle = theFig; 
    reconConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', theAxes, ...
    'activation', reshape(multistartStruct.reconPreds(:,ii),1,1,length(forwardExcitationsToStimulusUse)), ...
    'activationRange', 1.1*[0 max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)])], ...
    'plotTitle',  'Recon excitations');

    % Make sure excitations used match what comes back from multistart
    if (any(forwardExcitationsToStimulusUse * scaleFactor ~= multistartStruct.coneVec))
        error('Inconsistency in excitations driving reconstruction');
    end

    % Compute recon excitations to stimulus and compare with
    % scaled forward excitations to stimulus.
    subplot(3,6, 13); hold on; 
    if (pr.reconstructfromRenderMatrix)
        title({'Recon excittions to stim' ; 'Excitations from render matrix'});
        forwardExcitationsToReconTemp = forwardRenderMatrix*reconImageLinearTemp(:);
    else
        title({'Recon excittions to stim' ; 'Excitations from ISETBio'});
        forwardExcitationsToReconTemp = squeeze(forwardConeMosaic.Mosaic.compute(forwardOIToReconTemp, 'opticalImagePositionDegs', 'mosaic-centered'));
    end
    plot(forwardExcitationsToStimulusUse*scaleFactor,reconExcitationsToStimulusTemp,'ro','MarkerFaceColor','r','MarkerSize',6);
    axis('square');
    %minVal = 0.9*min([forwardExcitationsToStimulusUse*scaleFactor; forwardExcitationsToReconTemp]);
    %maxVal = 1.1*max([forwardExcitationsToStimulusUse*scaleFactor; forwardExcitationsToReconTemp]);
    minVal = 0.9*min([forwardExcitationsToStimulusUse*scaleFactor;reconExcitationsToStimulusTemp]);
    maxVal = 1.1*max([forwardExcitationsToStimulusUse*scaleFactor;reconExcitationsToStimulusTemp]);
    plot([minVal maxVal],[minVal maxVal],'k');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    xlabel('Scaled (pupil) forward excitations to stimulus');
    ylabel('Recon excitations to stimulus');

    % Compute forward excitations from reconstruction
    % and compare with scaled stimulus excitations
    subplot(3,6,14); hold on; 
    if (pr.reconstructfromRenderMatrix)
        title({'Forward excitations to recon' ; 'Excitations from render matrix'});
        forwardExcitationsToReconTemp = forwardRenderMatrix*reconImageLinearTemp(:);
    else
        title({'Forward excitations to recon' ; 'Excitations from ISETBio'});
        forwardExcitationsToReconTemp = squeeze(forwardConeMosaic.Mosaic.compute(forwardOIToReconTemp, 'opticalImagePositionDegs', 'mosaic-centered'));
    end
    %plot(forwardExcitationsToStimulusUse*scaleFactor,forwardExcitationsToReconTemp,'ro','MarkerFaceColor','r','MarkerSize',6);
    plot(forwardExcitationsToStimulusUse,forwardExcitationsToReconTemp,'ro','MarkerFaceColor','r','MarkerSize',6);
    axis('square');
    %minVal = 0.9*min([forwardExcitationsToStimulusUse*scaleFactor; forwardExcitationsToReconTemp]);
    %maxVal = 1.1*max([forwardExcitationsToStimulusUse*scaleFactor; forwardExcitationsToReconTemp]);
    minVal = 0.9*min([forwardExcitationsToStimulusUse; forwardExcitationsToReconTemp]);
    maxVal = 1.1*max([forwardExcitationsToStimulusUse; forwardExcitationsToReconTemp]);
    plot([minVal maxVal],[minVal maxVal],'k');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    xlabel('Unscaled (pupil) forward excitations to stimulus');
    ylabel('Forward excitations to recon');

    % Compute recon excitations from reconstruction
    % and compare with scaled stimulus excitations
    subplot(3,6,15); hold on;
    reconExcitationsToReconCheck = reconRenderMatrix*reconImageLinearTemp(:);
    if (pr.reconstructfromRenderMatrix)
        title({'Recon excitations to recon' ; 'Excitations from render matrix'});
        reconExcitationsToReconTemp = reconExcitationsToReconCheck;
    else
        title({'Recon excitations to recon' ; 'Excitations from ISETBio'});
        reconExcitationsToReconTemp = squeeze(reconConeMosaic.Mosaic.compute(reconOIToReconTemp, 'opticalImagePositionDegs', 'mosaic-centered'));
    end
    plot(forwardExcitationsToStimulusUse*scaleFactor,reconExcitationsToReconTemp,'ro','MarkerFaceColor','r','MarkerSize',6);
    axis('square');
    minVal = 0.9*min([forwardExcitationsToStimulusUse*scaleFactor; reconExcitationsToReconTemp]);
    maxVal = 1.1*max([forwardExcitationsToStimulusUse*scaleFactor; reconExcitationsToReconTemp]);
    plot([minVal maxVal],[minVal maxVal],'k');
    xlim([minVal maxVal]); ylim([minVal maxVal]);
    xlabel('Scaled (pupil) excitations to stimulus');
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
        saveas(gcf,fullfile(cnv.outputDir,sprintf('ReconExcitationsCheckError%d.jpg',ii)),'jpg');
        figure(theFig);
    end

    % Priors, likelihoods, and losses
    subplot(3,6,16);
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
    subplot(3,6,17);
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
    subplot(3,6,18);
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