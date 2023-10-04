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
%   03/26/2023  dhb  Wrote it from various pieces we already had.

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
reconDisplayScaleFactor = 3;

% Set image size and color
startWithLinearRGB = true;

nPixels = 128;
maxVal = 0.8;
bgVal = 0.2;
if (startWithLinearRGB)
    theInputLinearRGB = maxVal*[0 1 0] + bgVal*[1 1 1];
else
    theInputGammaCorrectedRGB = maxVal*[0 1 0] + bgVal*[1 1 1];
end

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
wls = (400:10:700)';
forwardDisplay = displaySet(forwardDisplay,'wave',wls);
reconDisplay = displaySet(reconDisplay,'wave',wls);

% Scale recon display primaries to try to keep things in range
reconDisplay = displaySet(reconDisplay,'spd primaries',displayGet(reconDisplay,'spd primaries')*reconDisplayScaleFactor);

% If we started with linear, gamma correct so we're all on the same
% page
if (startWithLinearRGB)
    theInputGammaCorrectedRGB = gammaCorrection(theInputLinear,forwardDisplay);
end

% Create an ISETBio scene.  Rescale input image
% according to pr.inputImageScaleFactor.
theImageRGB = ones(nPixels,nPixels,3);
for cc = 1:3
    theImageRGB(:,:,cc) = theInputGammaCorrectedRGB(cc)*theImageRGB(:,:,cc);
end
meanLuminanceCdPerM2 = [];
[theForwardScene, ~, theForwardImagergb] = sceneFromFile(theImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, forwardDisplay);
theForwardImageRGB = gammaCorrection(theForwardImagergb, forwardDisplay);
figure; clf; imshow(theForwardImageRGB);


% [theForwardScene, ~, theForwardImageLinear] = sceneFromFile(theForwardImageRGB, 'rgb', ...
%     meanLuminanceCdPerM2, forwardConeMosaic.Display);
% imwrite(theForwardImageRGB,fullfile(cnv.outputDir,'Stimulus.tiff'),'tiff');


% Get information we need to render scenes from their spectra through
% the recin display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mforward_rgbToXYZ = T_XYZ*displayGet(forwardDisplay,'spd primaries')*(wls(2)-wls(1));
Mforward_XYZTorgb = inv(Mforward_rgbToXYZ);
Mrecon_rgbToXYZ = T_XYZ*displayGet(reconDisplay,'spd primaries')*(wls(2)-wls(1));
Mrecon_XYZTorgb = inv(Mrecon_rgbToXYZ);

% Render the scene from the forward display on the recon display, to try to
% match XYZ.
[theForwardImagergbCalFormat,m,n] = ImageToCalFormat(theForwardImagergb);
theForwardImageXYZCalFormat = Mforward_rgbToXYZ*theForwardImagergbCalFormat;
theReconImagergbCalFormat = Mrecon_XYZTorgb*theForwardImageXYZCalFormat;
theReconImagergb = CalFormatToImage(theReconImagergbCalFormat,m,n);

% Truncate
minr = min(min(theReconImagergb(:,:,1)));
ming = min(min(theReconImagergb(:,:,2)));
minb = min(min(theReconImagergb(:,:,3)));
maxr = max(max(theReconImagergb(:,:,1)));
maxg = max(max(theReconImagergb(:,:,2)));
maxb = max(max(theReconImagergb(:,:,3)));
fprintf('Min: %0.2f, %0.2f, %0.2f\n',minr,ming,minb);
fprintf('Max: %0.2f, %0.2f, %0.2f\n',maxr,maxg,maxb);
theReconImagergbTruncated = theReconImagergb;
theReconImagergbTruncated(theReconImagergbTruncated < 0) = 0;
theReconImagergbTruncated(theReconImagergbTruncated > 1) = 1;

theReconImageTruncatedRGB = gammaCorrection(theReconImagergbTruncated, reconDisplay);
figure; clf; imshow(theReconImageTruncatedRGB);
