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
%   Only approximate, because many images on the monochromatic monitor
%   are outside the gamut of the conventional monitor.  Still, probably
%   more accurate than just looking at the RGB values for the monochromatic
%   monitor as a directly viewed image.

% History:
%   03/26/2023  dhb  Wrote it fredit Computeom various pieces we already had.
%   10/11/2023  dhb  More diagnostics, and make linear vs gamma corrected
%                    more explict.

% Clear
clear; close all;

% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
forwardDisplayName = 'mono';
reconDisplayName = 'conventional';
displayGammaBits = 12;
displayGammaGamma = 2;
overwriteDisplayGamma = true;
reconDisplayScaleFactor = 1;
SRGB = true;

% Set image size
nPixels = 128;

% Set a list of linear RGB
startWithLinearRGB = true;
inputLinearrgbValues = ...
   [0.1268    0.0971    0.0739   (0.0739 + 0.0536)/2  0.0536   0.0335    0.0255    0.0184    0.0123    0.0077    0.0040    0.0016    0.0001    0.0001 ;
    0.0240    0.0270    0.0299   (0.0299 + 0.0328)/2  0.0328   0.0365    0.0385    0.0404    0.0426    0.0445    0.0467    0.0489    0.0511    0.0517 ;
    0.0001    0.0001    0.0001   0.0001                        0.001     0.0001    0.0001    0.0001    0.0001    0.0001    0.0001    0.0001    0.0001    0.0001];

% Set up display specific fields
switch (forwardDisplayName)
    case 'conventional'
        forwardDisplayFieldName = 'CRT12BitDisplay';
    case 'mono'
        forwardDisplayFieldName = 'monoDisplay';
    otherwise
        error('Unknown forward display specified');
end
switch (reconDisplayName)
    case 'conventional'
        reconDisplayFieldName = 'CRT12BitDisplay';
    case 'mono'
        reconDisplayFieldName = 'monoDisplay';
    otherwise
        error('Unknown recon display specified');
end

% Get displays
forwardDisplayLoad = load(fullfile(aoReconDir, 'displays', [forwardDisplayName 'Display.mat']));
eval(['forwardDisplay = forwardDisplayLoad.' forwardDisplayFieldName ';']);
reconDisplayLoad = load(fullfile(aoReconDir, 'displays', [reconDisplayName 'Display.mat']));
eval(['reconDisplay = reconDisplayLoad.' reconDisplayFieldName ';']);
clear forwardDisplayLoad reconDisplayLoad

% Fix up gamma and wavelengths sampling
if (overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    forwardDisplay.gamma = gammaOutput(:,[1 1 1]);
    reconDisplay.gamma = gammaOutput(:,[1 1 1]);
end
wls = (400:1:700)';
forwardDisplay = displaySet(forwardDisplay,'wave',wls);
reconDisplay = displaySet(reconDisplay,'wave',wls);

% Scale recon display primaries to try to keep things in range
reconDisplay = displaySet(reconDisplay,'spd primaries',displayGet(reconDisplay,'spd primaries')*reconDisplayScaleFactor);

% Get information we need to render scenes from their spectra through
% the recon display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mforward_rgbToXYZ = T_XYZ*displayGet(forwardDisplay,'spd primaries')*(wls(2)-wls(1));
Mforward_XYZTorgb = inv(Mforward_rgbToXYZ);
Mrecon_rgbToXYZ = T_XYZ*displayGet(reconDisplay,'spd primaries')*(wls(2)-wls(1));
Mrecon_XYZTorgb = inv(Mrecon_rgbToXYZ);

% Report on primaries
fprintf('Forward primary luminances (r, g, b): %0.2f, %0.2f %0.2f\n',Mforward_rgbToXYZ(2,1),Mforward_rgbToXYZ(2,2),Mforward_rgbToXYZ(2,3));
fprintf('Recon primary luminances (r, g, b): %0.2f, %0.2f %0.2f\n',Mrecon_rgbToXYZ(2,1),Mrecon_rgbToXYZ(2,2),Mrecon_rgbToXYZ(2,3));

% Loop over all the input values
for iii = 1:size(inputLinearrgbValues,2)
    theFigure{iii} = figure; clf;
    set(gcf,'Position',[100 100 1200 450]);

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
    
    minrForward(iii) = min(min(theForwardImagergb{iii}(:,:,1)));
    mingForward(iii) = min(min(theForwardImagergb{iii}(:,:,2)));
    minbForward(iii) = min(min(theForwardImagergb{iii}(:,:,3)));
    maxrForward(iii) = max(max(theForwardImagergb{iii}(:,:,1)));
    maxgForward(iii) = max(max(theForwardImagergb{iii}(:,:,2)));
    maxbForward(iii) = max(max(theForwardImagergb{iii}(:,:,3)));
    maxForward(iii) = max(theForwardImagergb{iii}(:));
    fprintf('Min input mono linear: %0.4f, %0.4f, %0.4f\n',minrForward(iii),mingForward(iii),minbForward(iii));
    fprintf('Max input mono linear: %0.4f, %0.4f, %0.4f\n',maxrForward(iii),maxgForward(iii),maxbForward(iii));
   
    % Render the scene from the forward display on the recon display, to try to
    % match XYZ.
    [theForwardImagergbCalFormat,m,n] = ImageToCalFormat(theForwardImagergb{iii});
    theForwardImageXYZCalFormat = Mforward_rgbToXYZ*theForwardImagergbCalFormat;

    % Either use SRGB or the recon display, as specified
    if (SRGB)
        theReconImagergbCalFormat = XYZToSRGBPrimary(theForwardImageXYZCalFormat);
    else
        theReconImagergbCalFormat = Mrecon_XYZTorgb*theForwardImageXYZCalFormat;
    end
    theReconImagergb{iii} = CalFormatToImage(theReconImagergbCalFormat,m,n);

    % Check/plot chromaticities
    theReconImageXYZChkCalFormat = Mrecon_rgbToXYZ*theReconImagergbCalFormat;
    theForwardImagexyYCalFormat = XYZToxyY(theForwardImageXYZCalFormat);
    theReconImagexyYChkCalFormat = XYZToxyY(theReconImageXYZChkCalFormat);
    theForwardxyY = XYZToxyY(Mforward_rgbToXYZ);
    theReconxyY = XYZToxyY(Mrecon_rgbToXYZ);
    T_xyY = XYZToxyY(T_XYZ);
    subplot(1,3,3); hold on
    plot(theForwardxyY(1,:),theForwardxyY(2,:),'rs','MarkerFaceColor','r','MarkerSize',12);
    plot(theReconxyY(1,:),theReconxyY(2,:),'gs','MarkerFaceColor','g','MarkerSize',12);
    plot(theForwardImagexyYCalFormat(1,:),theForwardImagexyYCalFormat(2,:),'ro','MarkerFaceColor','r','MarkerSize',8);
    plot(theReconImagexyYChkCalFormat(1,:),theReconImagexyYChkCalFormat(2,:),'go','MarkerFaceColor','g','MarkerSize',6);
    plot(T_xyY(1,:),T_xyY(2,:),'k','LineWidth',2);
    xlim([0 1]); ylim([0 1]);

    % Truncate rendered image if needed
    minrRecon(iii) = min(min(theReconImagergb{iii}(:,:,1)));
    mingRecon(iii) = min(min(theReconImagergb{iii}(:,:,2)));
    minbRecon(iii) = min(min(theReconImagergb{iii}(:,:,3)));
    maxrRecon(iii) = max(max(theReconImagergb{iii}(:,:,1)));
    maxgRecon(iii) = max(max(theReconImagergb{iii}(:,:,2)));
    maxbRecon(iii) = max(max(theReconImagergb{iii}(:,:,3)));
    maxRecon(iii) = max(theReconImagergb{iii}(:));
    fprintf('Min conventional display linear: %0.4f, %0.4f, %0.4f\n',minrRecon(iii),mingRecon(iii),minbRecon(iii));
    fprintf('Max conventional display linear: %0.4f, %0.4f, %0.4f\n',maxrRecon(iii),maxgRecon(iii),maxbRecon(iii));
    theReconImagergbTruncated{iii} = theReconImagergb{iii};
    theReconImagergbTruncated{iii}(theReconImagergbTruncated{iii} < 0) = 0;
    theReconImageTruncatedRGBTemp = gammaCorrection(theReconImagergbTruncated{iii}, reconDisplay);

    theReconImagergbTruncatedCalFormat{iii} = ImageToCalFormat(theReconImagergbTruncated{iii});
    if (SRGB)
        theReconImageXYZTruncatedCalFormat = SRGBPrimaryToXYZ(theReconImagergbTruncatedCalFormat{iii});
    else
        theReconImageXYZTruncatedCalFormat = Mrecon_rgbToXYZ*theReconImagergbTruncatedCalFormat{iii};
    end
    theReconImagexyYTruncatedCalFormat = XYZToxyY(theReconImageXYZTruncatedCalFormat);
    subplot(1,3,3);
    plot(theReconImagexyYTruncatedCalFormat(1,:),theReconImagexyYTruncatedCalFormat(2,:),'bo','MarkerFaceColor','b','MarkerSize',4);

    % Do it with our function.
    %
    % Note that we've incorporated the viewingDisplayScaleFactor into the
    % reconDisplay object already, so when we call the function we set the
    % factor to one.
    [outputImageRGBFromFunction,outputImagergbFromFunction] = rgb2aoDisplay(theImageRGB, forwardDisplay, reconDisplay, ...
        'viewingDisplayScaleFactor',1,'linearInput',false,'wls',wls,'verbose',true);
    if (max(abs(theReconImageTruncatedRGBTemp(:)-outputImageRGBFromFunction(:))) > 1e-6)
        %error('Not recreating same output image through our function');
    end

    % Get the equivalent wavelength with our function.  Picking off the
    % second return argument for a one pixel input gives us just one number
    % back.
    theImageForEquivRGB = zeros(1,1,3);
    theImageForEquivRGB(1,1,:) = theInputGammaCorrectedRGB';
    [~,eqWavelengthCal]= rgb2EquivWavelength(theImageForEquivRGB,forwardDisplay,'wls',wls);

    % Add to chromaticity plot.
    index = find(wls == eqWavelengthCal);
    subplot(1,3,3);
    eewxy = [0.33 0.33]';
    plot(T_xyY(1,index),T_xyY(2,index),'cx','MarkerSize',16);
    plot(eewxy(1),eewxy(2),'kx','MarkerSize',16);
    plot([eewxy(1) T_xyY(1,index)],[eewxy(2) T_xyY(2,index)],'k','LineWidth',2);
    title(sprintf('Equiv wavelength: %d',eqWavelengthCal));
    drawnow;

    fprintf('\n')
end

% Scale the images jointly and show, and write out each image for happy use
% later.
outputDir = fullfile(aoReconDir, 'displays','tutorialOutput');
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
maxForwardAll = max(maxForward);
maxReconAll = max(maxRecon);
for iii = 1:size(inputLinearrgbValues,2)
    figure(theFigure{iii});

    % Show forward image
    subplot(1,3,1);
    theForwardImageRGB = gammaCorrection(theForwardImagergb{iii}/maxForwardAll, forwardDisplay);
    imshow(theForwardImageRGB);
    title({ 'Direct on Forward' ; sprintf('rgb: %0.2f %0.2f %0.2f', ...
        minrForward(iii)/maxForwardAll,mingForward(iii)/maxForwardAll,minbForward(iii)/maxForwardAll) });

    % Gamma correct and show recon image
    subplot(1,3,2);
    if (SRGB)
        theReconImageTruncatedRGB = SRGBGammaCorrect(theReconImagergbTruncated{iii}/maxReconAll,0);
        imshow(uint8(theReconImageTruncatedRGB));
        imwrite(uint8(theReconImageTruncatedRGB),fullfile(outputDir,sprintf('Stimulus%d.tiff',iii)),'tiff');
    else
        theReconImageTruncatedRGB = gammaCorrection(theReconImagergbTruncated{iii}/maxReconAll, reconDisplay);
        imshow(theReconImageTruncatedRGB);
        imwrite(uint8(theReconImageTruncatedRGB),fullfile(outputDir,sprintf('Stimulus%d.tiff',iii)),'tiff');
    end
    title({ 'Render on Recon' ; sprintf('Stimulus %d',iii) ; sprintf('rgb: %0.4f %0.4f %0.4f', ...
        minrRecon(iii)/maxReconAll,mingRecon(iii)/maxReconAll,minbRecon(iii)/maxReconAll) });
  
end