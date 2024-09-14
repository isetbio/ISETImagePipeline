% t_renderMonDisplayImage
%
% Description:
%   Shows how to visualize a colorimetric equivlanet of what
%   we put onto the monochromatic monitor, so that the displayed
%   image gives a sense of what it would look like to a person
%   who was looking at the monochromatic monitor.
%
%   Works by converting the image on the monochromatic monitor to
%   XYZ, and then finding settings on conventional monitor that
%   produce that XYZ.
%
%   Only approximate, because some images on the monochromatic monitor
%   are outside the gamut of the conventional monitor.  Still, probably
%   more accurate than just looking at the RGB values for the monochromatic
%   monitor as a directly viewed image.  The plots illustrate the
%   chromaticiites of what happend.
%
%   Also computes the equivalent wavelength of the stimuli.
%
%   Some of the calculations are done here, but this illustrates the use of
%   functions RGBRenderAcrossDisplays and RGBToEquivWavelength, which are
%   where the work would be done in applications of these ideas.
%
%   Note that some thought is needed about how to scale the images for
%   display.  Here a common scaling is applied to all of the images within
%   each monitor type.

% History:
%   03/26/2023  dhb  Wrote it from various pieces we already had.
%   10/11/2023  dhb  More diagnostics, and make linear vs gamma corrected
%                    more explict.
%               dhb  Add equivalent wavelength calc.
%   10/12/2023  dhb  Clean up.

% Clear
clear; close all;

% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
forwardDisplayName = 'mono';
renderDisplayName = 'conventional';
displayGammaBits = 12;
displayGammaGamma = 2;
overwriteDisplayGamma = true;

% The first scale factor makes the stimulus luminance about 2800 cd/m2.
% The second allows adjustment of relative intensity range of the two monitors.
overallDisplayScaleFactor = 2.5*(2800/4987)*194.2;
renderRelativeDisplayScaleFactor = 1;
baseDisplayScaleVector = [1 0.6 0.8];

% Render in sRGB (true) or with respect to ISETBio viewing display (false)
SRGB = false;

% Set image size
nPixels = 128;

% Set up display specific fields
switch (forwardDisplayName)
    case 'conventional'
        forwardDisplayFieldName = 'CRT12BitDisplay';
    case 'mono'
        forwardDisplayFieldName = 'monoDisplay';
    otherwise
        error('Unknown forward display specified');
end
switch (renderDisplayName)
    case 'conventional'
        renderDisplayFieldName = 'CRT12BitDisplay';
    case 'mono'
        renderDisplayFieldName = 'monoDisplay';
    otherwise
        error('Unknown render display specified');
end

% Get displays
forwardDisplayLoad = load(fullfile(aoReconDir, 'displays', [forwardDisplayName 'Display.mat']));
eval(['forwardDisplay = forwardDisplayLoad.' forwardDisplayFieldName ';']);
renderDisplayLoad = load(fullfile(aoReconDir, 'displays', [renderDisplayName 'Display.mat']));
eval(['renderDisplay = renderDisplayLoad.' renderDisplayFieldName ';']);
clear forwardDisplayLoad renderDisplayLoad
% Fix up gamma, wavelength, sampling and display scaling.
if (overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    forwardDisplay = displaySet(forwardDisplay,'gamma',gammaOutput(:,[1 1 1]));
    renderDisplay = displaySet(renderDisplay,'gamma',gammaOutput(:,[1 1 1]));
end
wls = (400:1:700)';
forwardDisplay = displaySet(forwardDisplay,'wave',wls);
renderDisplay = displaySet(renderDisplay,'wave',wls);
forwardDisplay = ScaleDisplayPrimaries(forwardDisplay,overallDisplayScaleFactor*baseDisplayScaleVector);
renderDisplay = ScaleDisplayPrimaries(renderDisplay,renderRelativeDisplayScaleFactor*overallDisplayScaleFactor*baseDisplayScaleVector);

% Make sure black is zero for our calculations
forwardAmbient = displayGet(forwardDisplay,'black spd');
if (any(forwardAmbient ~= 0))
    error('Forward display ambient not all zeros');
end
renderAmbient = displayGet(renderDisplay,'black spd');
if (any(renderAmbient ~= 0))
    error('Render ambient not all zero');
end

% Get information we need to render scenes from their spectra through
% the render display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mforward_rgbToXYZ = T_XYZ*displayGet(forwardDisplay,'spd primaries')*(wls(2)-wls(1));
Mforward_XYZTorgb = inv(Mforward_rgbToXYZ);
Mrender_rgbToXYZ = T_XYZ*displayGet(renderDisplay,'spd primaries')*(wls(2)-wls(1));
Mrender_XYZTorgb = inv(Mrender_rgbToXYZ);

% Report on primaries
rPrimaryLuminance = Mforward_rgbToXYZ(2,1);
gPrimaryLuminance = Mforward_rgbToXYZ(2,2);
fprintf('Forward primary luminances (r, g, b): %0.2f, %0.2f %0.2f\n',Mforward_rgbToXYZ(2,1),Mforward_rgbToXYZ(2,2),Mforward_rgbToXYZ(2,3));
fprintf('Render primary luminances (r, g, b): %0.2f, %0.2f %0.2f\n',Mrender_rgbToXYZ(2,1),Mrender_rgbToXYZ(2,2),Mrender_rgbToXYZ(2,3));

% Get input rgb values to analyze.  Can construct or use table.
useEquilumConstruct = false;
targetEqWls = [650 625 600 590 582 578 570 560 545];
targetLuminance = 2800;
if (useEquilumConstruct)
    % Compute as set of equally spaced r/(r+g) values that lead
    % to equal luminance stimuli of specified luminance.
    startWithLinearRGB = true;
    nEquiLumStimuli = 80;
    rOverRPlusG = linspace(1,0,nEquiLumStimuli);
    gOverRPlusG = 1-rOverRPlusG;
    gPrimaryAdjust = rPrimaryLuminance/gPrimaryLuminance;
    for nn = 1:nEquiLumStimuli
        rRaw = rOverRPlusG(nn);
        gRaw = gOverRPlusG(nn);
        gAdjust = gRaw*gPrimaryAdjust;
        equiLumrgbValuesRaw = [rRaw ; gAdjust; 0];
        equiLumLuminanceRaw(nn) = rRaw*rPrimaryLuminance + gAdjust*gPrimaryLuminance;
        equiLumrgbValues = (targetLuminance/equiLumLuminanceRaw(nn))*equiLumrgbValuesRaw;
        inputLinearrgbValues(:,nn) = equiLumrgbValues;
    end
else
    % Set a list of linear RGB to render as uniform swatches.  These
    % were calculated at some earlier point using the useEquilumContrast option
    % set to true, and then selected to have reasonable appearance spacing by
    % hand. We think these were selected from 40 stimuli generated in the
    % section above and then sampled.
    startWithLinearRGB = true;
    inputLinearrgbValues = [ ...
            0.5941    0.5159    0.3518    0.2658    0.2033    0.1720    0.1173    0.0704    0.0078
            0.0082    0.0355    0.0929    0.1230    0.1449    0.1558    0.1749    0.1913    0.2132
            0         0         0         0         0         0         0         0         0
         ];
        % (0.8/0.4615) * ...
        % [0.4615 0.3846 0.3077 0.2308 0.1538 0.1154 0.0769 0.0385 0.0000;
        % 0.0081 0.0242 0.0403 0.0565 0.0726 0.0807 0.0888 0.0968 0.1049;
        % 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005];
end
inputLinearrgbValues = [inputLinearrgbValues [0.4 0.4 0.0]' [0.4 0.4 0.4]'];
nExtraStim = 2;

% Loop over all the input values
for iii = 1:size(inputLinearrgbValues,2)
    if (~useEquilumConstruct)
        theFigure{iii} = figure; clf;
        set(gcf,'Position',[100 100 1200 450]);
    end

    % Get linear rgb
    theInputLinearrgb = inputLinearrgbValues(:,iii);

    % If we started with linear, gamma correct so we're all on the same page
    if (startWithLinearRGB)
        theInputGammaCorrectedRGB = gammaCorrection(theInputLinearrgb,forwardDisplay);
    end

    % Create an ISETBio scene.  Rescale input image
    % according to pr.inputImageScaleFactor.
    theImageRGB = ones(nPixels,nPixels,3);
    for cc = 1:3
        theImageRGB(:,:,cc) = theInputGammaCorrectedRGB(cc)*theImageRGB(:,:,cc);
    end
    meanLuminanceCdPerM2 = [];
    [theForwardScene, ~, theForwardImagergb{iii}] = sceneFromFile(theImageRGB, 'rgb', ...
        meanLuminanceCdPerM2, forwardDisplay);
    
    therForward(iii) = min(min(theForwardImagergb{iii}(:,:,1)));
    thegForward(iii) = min(min(theForwardImagergb{iii}(:,:,2)));
    thebForward(iii) = min(min(theForwardImagergb{iii}(:,:,3)));
    maxrForward(iii) = max(max(theForwardImagergb{iii}(:,:,1)));
    maxgForward(iii) = max(max(theForwardImagergb{iii}(:,:,2)));
    maxbForward(iii) = max(max(theForwardImagergb{iii}(:,:,3)));
    maxForward(iii) = max(theForwardImagergb{iii}(:));
    if (therForward(iii) ~= maxrForward(iii) | thegForward(iii) ~= maxgForward(iii) | thebForward(iii) ~= maxbForward(iii))
        error('Image not spatially uniform as expected');
    end
    fprintf('Input mono linear rgb: %0.4f, %0.4f, %0.4f\n',therForward(iii),thegForward(iii),thebForward(iii));
   
    % Render the scene from the forward display on the render display, to try to
    % match XYZ.
    samplePixel = round(nPixels^2/2);
    [theForwardImagergbCalFormat,m,n] = ImageToCalFormat(theForwardImagergb{iii});
    theForwardImageXYZCalFormat = Mforward_rgbToXYZ*theForwardImagergbCalFormat;
    theForwardImagexyYCalFormat = XYZToxyY(theForwardImageXYZCalFormat);
    theForwardLuminance(iii) = theForwardImageXYZCalFormat(2,samplePixel);
    theForwardxyY(:,iii) = theForwardImagexyYCalFormat(:,samplePixel);

    % Either use SRGB or the render display, as specified
    if (SRGB)
        theRenderImagergbCalFormat = XYZToSRGBPrimary(theForwardImageXYZCalFormat);
    else
        theRenderImagergbCalFormat = Mrender_XYZTorgb*theForwardImageXYZCalFormat;
    end
    theRenderImagergb{iii} = CalFormatToImage(theRenderImagergbCalFormat,m,n);
    theRenderImageXYZCalFormat = Mrender_rgbToXYZ*theRenderImagergbCalFormat;
    theRenderImagexyYCalFormat = XYZToxyY(theRenderImageXYZCalFormat);
    theRenderLuminance(iii) = theRenderImageXYZCalFormat(2,samplePixel);
    theRenderxyY(:,iii) = theForwardImagexyYCalFormat(:,samplePixel);

    % Plot chromaticities
    theForwardImagexyYCalFormat = XYZToxyY(theForwardImageXYZCalFormat);
    theForwardDisplayxyY = XYZToxyY(Mforward_rgbToXYZ);
    theRenderDisplayxyY = XYZToxyY(Mrender_rgbToXYZ);
    T_xyY = XYZToxyY(T_XYZ);
    if (~useEquilumConstruct)
        subplot(1,3,3); hold on
        plot(theForwardDisplayxyY(1,:),theForwardDisplayxyY(2,:),'rs','MarkerFaceColor','r','MarkerSize',12);
        plot(theRenderDisplayxyY(1,:),theRenderDisplayxyY(2,:),'gs','MarkerFaceColor','g','MarkerSize',12);
        plot(theForwardImagexyYCalFormat(1,:),theForwardImagexyYCalFormat(2,:),'go','MarkerFaceColor','g','MarkerSize',12);
        plot(T_xyY(1,:),T_xyY(2,:),'k','LineWidth',2);
        xlim([0 1]); ylim([0 1]);
    end

    % Truncate rendered image if needed
    therRender(iii) = min(min(theRenderImagergb{iii}(:,:,1)));
    thegRender(iii) = min(min(theRenderImagergb{iii}(:,:,2)));
    thebRender(iii) = min(min(theRenderImagergb{iii}(:,:,3)));
    maxrRender(iii) = max(max(theRenderImagergb{iii}(:,:,1)));
    maxgRender(iii) = max(max(theRenderImagergb{iii}(:,:,2)));
    maxbRender(iii) = max(max(theRenderImagergb{iii}(:,:,3)));
    maxRender(iii) = max(theRenderImagergb{iii}(:));
    if (therRender(iii) ~= maxrRender(iii) | thegRender(iii) ~= maxgRender(iii) | thebRender(iii) ~= maxbRender(iii))
        error('Image not spatially uniform as expected');
    end
    fprintf('Conventional display rgb linear: %0.4f, %0.4f, %0.4f\n',therRender(iii),thegRender(iii),thebRender(iii));
    theRenderImagergbTruncated{iii} = theRenderImagergb{iii};
    theRenderImagergbTruncated{iii}(theRenderImagergbTruncated{iii} < 0) = 0;

    % Convert truncated linear rgb image back to XYZ so we can look in
    % chromaticity coordinates, and add to plot
    theRenderImagergbTruncatedCalFormat{iii} = ImageToCalFormat(theRenderImagergbTruncated{iii});
    if (SRGB)
        theRenderImageXYZTruncatedCalFormat = SRGBPrimaryToXYZ(theRenderImagergbTruncatedCalFormat{iii});
    else
        theRenderImageXYZTruncatedCalFormat = Mrender_rgbToXYZ*theRenderImagergbTruncatedCalFormat{iii};
    end
    theRenderImagexyYTruncatedCalFormat = XYZToxyY(theRenderImageXYZTruncatedCalFormat);
    theRenderxyYTruncated(:,iii) = theRenderImagexyYTruncatedCalFormat(:,samplePixel);

    if (~useEquilumConstruct)
        subplot(1,3,3);
        plot(theRenderImagexyYTruncatedCalFormat(1,:),theRenderImagexyYTruncatedCalFormat(2,:),'bo','MarkerFaceColor','b','MarkerSize',12);
    end

    % Info
    fprintf('Forward stimulus xyY:          %0.4f, %0.4f, %0.2f cd/m2\n',theForwardxyY(1,iii),theForwardxyY(2,iii),theForwardxyY(3,iii));
    fprintf('Render stimulus xyY:           %0.4f, %0.4f, %0.2f cd/m2\n',theRenderxyY(1,iii),theRenderxyY(2,iii),theRenderxyY(3,iii));
    fprintf('Render stimulus truncated xyY: %0.4f, %0.4f, %0.2f cd/m2\n',theRenderxyYTruncated(1,iii),theRenderxyYTruncated(2,iii),theRenderxyYTruncated(3,iii));

    % Get the equivalent wavelength with our function.  Picking off the
    % second return argument for a one pixel input gives us just one number
    % back.
    theImageForEquivRGB = zeros(1,1,3);
    theImageForEquivRGB(1,1,:) = theInputGammaCorrectedRGB';
    [~,forwardEqWavelength(iii)]= RGBToEquivWavelength(theImageForEquivRGB,forwardDisplay,'wls',wls);
        
    % Do it with our function.
    %
    % Note that we've incorporated the renderDisplayScaleFactor into the
    % renderDisplay object already, so when we call the function we set the
    % factor to one.
    if (~SRGB)
        % Note that the sRGB RGB image comes back as a double in the range
        % [0-1].  Below we leave it as a uint8 in range [0-255].  Also note
        % that the scaling of the rgb image is not necessarily in a good
        % range here.  Below we illustrate how to apply a common scaling
        % prior to gamma correction.
        [outputImageRGBFromFunction,outputImagergbFromFunction] = RGBRenderAcrossDisplays(theImageRGB, forwardDisplay, renderDisplay, ...
            'viewingDisplayScaleFactor',1,'linearInput',false,'wls',wls,'verbose',true,'scaleToMax',false,'SRGB',SRGB);
    else
        % Same comment on scaling applies as for sRGB branch above.
        [outputImageRGBFromFunction,outputImagergbFromFunction] = RGBRenderAcrossDisplays(theImageRGB, forwardDisplay, [], ...
            'viewingDisplayScaleFactor',1,'linearInput',false,'wls',wls,'verbose',true,'scaleToMax',false,'SRGB',SRGB);
    end

    % Reality check.
    if (max(abs(theRenderImagergbTruncated{iii}(:)-outputImagergbFromFunction(:))) > 1e-6)
        error('Not recreating same linear rgb output image through our function');
    end

    % Get the equivalent wavelength with our function.  Picking off the
    % second return argument for a one pixel input gives us just one number
    % back.
    theImageForEquivRGB = zeros(1,1,3);
    theImageForEquivRGB(1,1,:) = theInputGammaCorrectedRGB';
    [~,forwardEqWavelength(iii)]= RGBToEquivWavelength(theImageForEquivRGB,forwardDisplay,'wls',wls);
    fprintf('Equivalent wavelength: %d\n',forwardEqWavelength(iii));

    % Add to chromaticity plot.
    %
    % This plot shows:
    %   open red squares:   gamut of forward monitor
    %   open green squares: gamut of render monitor (sRGB gamut not shown)
    %   black cross:        EEW chromaticity
    %   cyan cross:         Equiv wavelength chromaticity (hard to see)
    %   black line:         Connects EEW and equiv wl points
    %   green circle:       original image chromacitity
    %   blue circle:        Rendered image chromaticity (might look black)
    index = find(wls == forwardEqWavelength(iii));

    if (~useEquilumConstruct)
        subplot(1,3,3);
        eewxy = [0.33 0.33]';
        plot(T_xyY(1,index),T_xyY(2,index),'cx','MarkerSize',16);
        plot(eewxy(1),eewxy(2),'kx','MarkerSize',16);
        plot([eewxy(1) T_xyY(1,index)],[eewxy(2) T_xyY(2,index)],'k','LineWidth',2);
        title(sprintf('Equiv wavelength: %d',forwardEqWavelength(iii)));
        drawnow;
    end
    fprintf('\n')
end


%% Compute background rgb
backgroundxyY = [0.33 0.33 26]';
backgroundXYZ = xyYToXYZ(backgroundxyY);
backgroundrgbLinear = Mforward_XYZTorgb*backgroundXYZ;
fprintf('Forward stimulus background xyY: %0.4f, %0.4f, %0.2f cd/m2\n',backgroundxyY(1),backgroundxyY(2),backgroundxyY(3));
fprintf('Forward stimulus rgbLinear: %0.5f, %0.5f, %0.5f\n',backgroundrgbLinear(1),backgroundrgbLinear(2),backgroundrgbLinear(3));

% Scale the images jointly and show, and write out each image for happy use
% later.
outputDir = fullfile(aoReconDir, 'displays','tutorialOutput');
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
maxForwardAll = max(maxForward);
maxRenderAll = max(maxRender);
for iii = 1:size(inputLinearrgbValues,2)
    if (~useEquilumConstruct)
        figure(theFigure{iii});

        % Show forward image
        subplot(1,3,1);
        theForwardImageRGB = gammaCorrection(theForwardImagergb{iii}/maxForwardAll, forwardDisplay);
        imshow(theForwardImageRGB);
        title({ 'Direct on Forward' ; sprintf('rgb: %0.2f %0.2f %0.2f', ...
            therForward(iii)/maxForwardAll,thegForward(iii)/maxForwardAll,thebForward(iii)/maxForwardAll) });

        % Gamma correct and show/write rendered image
        subplot(1,3,2);
        if (SRGB)
            theRenderImageTruncatedRGB = SRGBGammaCorrect(theRenderImagergbTruncated{iii}/maxRenderAll,0);
            imshow(uint8(theRenderImageTruncatedRGB));
            imwrite(uint8(theRenderImageTruncatedRGB),fullfile(outputDir,sprintf('Stimulus%d_SRGB.tiff',iii)),'tiff');
        else
            theRenderImageTruncatedRGB = gammaCorrection(theRenderImagergbTruncated{iii}/maxRenderAll, renderDisplay);
            imshow(theRenderImageTruncatedRGB);
            imwrite(uint8(theRenderImageTruncatedRGB),fullfile(outputDir,sprintf('Stimulus%d_ConvMonitor.tiff',iii)),'tiff');
        end
        title({ 'Render on Render' ; sprintf('Stimulus %d',iii) ; sprintf('rgb: %0.4f %0.4f %0.4f', ...
            therRender(iii)/maxRenderAll,thegRender(iii)/maxRenderAll,thebRender(iii)/maxRenderAll) });
    end
end

% If constructing, pull out stimuli that come as close as possible to
% producing the target equivalent wavelengths
if (useEquilumConstruct)
    for www = 1:length(targetEqWls)
        [~,minIndex] = min(abs(targetEqWls(www)-forwardEqWavelength(1:end-nExtraStim)));
        constructedInputLinearrgb(:,www) = inputLinearrgbValues(:,minIndex(1));
        constructedEqWls(www) = forwardEqWavelength(minIndex(1));
        constructedLuminances(www) = theForwardxyY(3,minIndex(1));
    end

    % Quick dump as a check
    targetEqWls
    constructedEqWls
    constructedLuminances
    constructedInputLinearrgb
end




