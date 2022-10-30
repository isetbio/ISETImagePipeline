% Get matrix that goes from XYZ to monitor rgb.  We need this to render the OI for the display
%
% Direct code.  Result is the same (yay!)
%   theXYZStruct = load('T_xyz1931');
%   wls = oiGet(forwardOI,'wave');
%   T_xyz = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
%   M_rgbToxyz = T_xyz*displayGet(forwardConeMosaic.Display,'spd')*(wls(2)-wls(1));
M_rgbToxyz = displayGet(forwardConeMosaic.Display,'rgb2xyz');
M_xyzTorgb = inv(M_rgbToxyz);

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

    % Render OIs on display as best we can.  Getting the scale factor right
    % from first principles is hard because it depends on the radiance ->
    % irradiance conversion and we don't want to unpack that here.  So we
    % scale both OI's we want to visualize to the same range.  We also
    % assume they are both the same size.
    [forwardOIxyz,m,n] = ImageToCalFormat(oiGet(forwardOI,'xyz'));
    forwardOIrgb = M_xyzTorgb*forwardOIxyz/255;
    fprintf('Min/max of forward OI image: %0.2f, %0.2f\n',min(forwardOIrgb(:)),max(forwardOIrgb(:)));
    [reconOIxyz,m,n] = ImageToCalFormat(oiGet(reconOIToReconTemp,'xyz'));
    reconOIrgb = M_xyzTorgb*reconOIxyz;
    fprintf('Min/max of recon OI image: %0.2f, %0.2f\n',min(reconOIrgb(:)),max(reconOIrgb(:)));
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

    [stimulusScenexyz,m,n] = ImageToCalFormat(sceneGet(stimulusScene,'xyz'));
    stimulusScenergb = M_xyzTorgb*stimulusScenexyz;
    fprintf('Min/max of stimluus image: %0.2f, %0.2f\n',min(stimulusScenergb(:)),max(stimulusScenergb(:)));
    stimulusScenergb(stimulusScenergb < 0) = 1;
    stimulusScenergb(stimulusScenergb > 1) = 1;
    stimulusSceneRGB1 = gammaCorrection(CalFormatToImage(stimulusScenergb,m,n),forwardConeMosaic.Display);

    % Optical image of stimulus
    theAxes = subplot(3,6,3);
    % visualizeOpticalImage(forwardOI, 'axesHandle',theAxes,'avoidAutomaticRGBscaling', true, 'displayRadianceMaps', false);
    imshow(forwardOIRGB);

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

    % Optical image of recon
    theAxes = subplot(3,6,9);
    % visualizeOpticalImage(reconOIToReconTemp, 'axesHandle',theAxes,'avoidAutomaticRGBscaling', true, 'displayRadianceMaps', false);
    imshow(reconOIRGB);

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
    'plotTitle',  'Scaled (pupil) recon excitations');

    % Make sure excitations used match what comes back from multistart
    if (any(forwardExcitationsToStimulusUse * scaleFactor ~= multistartStruct.coneVec))
        error('Inconsistency in excitations driving reconstruction');
    end

%     % Plot predicted from recon versus stim excitations
%     subplot(3,6,13); hold on;
%     minVal = 0.9*min([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)]);
%     maxVal = 1.1*max([multistartStruct.coneVec ; multistartStruct.reconPreds(:,ii)]);
%     plot(multistartStruct.coneVec,multistartStruct.reconPreds(:,ii),'ro','MarkerFaceColor','r','MarkerSize',6);
%     xlabel('Scaled (pupil) stimulus excitations');
%     ylabel('Recon excitations to recon');
%     xlim([minVal maxVal]); ylim([minVal maxVal]);
%     axis('square');

    % Compute forward excitations from reconstruction
    % and compare with scaled stimulus excitations
    subplot(3,6,14); hold on; hold on;
    if (pr.reconstructfromRenderMatrix)
        title('Recon from render matrix');
        forwardExcitationsToReconTemp = forwardRenderMatrix*reconImageLinearTemp(:);
    else
        title('Recon from ISETBio');
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
        title('Recon from render matrix');
        reconExcitationsToReconTemp = reconExcitationsToReconCheck;
    else
        title('Recon from ISETBio');
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