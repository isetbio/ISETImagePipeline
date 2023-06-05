function [outputImageRGB] = rgb2aoDisplay(inputImageRGB, startDisplay, viewingDisplay)
% Synopsis:
%    Correct input image to a form that better approximates how it would be
%    displayed on the start Display Monitor
%    inputImageRGB: RGB image matrix prior to gamma correction
%    startDisplayName: 'mono' or 'conventional'
%    viewingDisplayName: 'mono' or 'conventional'
%    viewingDisplayScaleFactor: Some convenience scaling, value = 3 works
%    aoReconDir, displayGammaBits, displayGammaGamma: Params from pr struct
%
% Description:
%    Input a linear image specifying which monitor the image is being viewed under
%    (ViewingDisp) and apply adjustments based on which monitor the
%    calculations were run under (startDisp). Ex: [1 1 0] RGB image appears
%    yellow under mono monitor for calculations but orangish under conventional.
%
% See also: t_renderMonoDisplayImage, aoStimRecon, aoStimReconRunMany

% History:
%   04/04/23  chr  Made callable function from t_renderMonoDisplayImage

%
% startDisplayName = 'conventional';
% viewingDisplayName = 'mono';
viewingDisplayScaleFactor = 1;

% Convenience assumption that all displays operate on the same wave structure
wls = (400:10:700)';
startDisplay = displaySet(startDisplay,'wave',wls);
viewingDisplay = displaySet(viewingDisplay,'wave',wls);

% Scale recon display primaries to try to keep things in range
viewingDisplay = displaySet(viewingDisplay,'spd primaries',displayGet(viewingDisplay,'spd primaries')*viewingDisplayScaleFactor);

meanLuminanceCdPerM2 = [];
[~, ~, startImageLinear] = sceneFromFile(inputImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, startDisplay);

% Get information we need to render scenes from their spectra through
% the viewing display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mstart_rgbToXYZ = T_XYZ*displayGet(startDisplay,'spd primaries')*(wls(2)-wls(1));
Mstart_XYZTorgb = inv(Mstart_rgbToXYZ);
Mviewing_rgbToXYZ = T_XYZ*displayGet(viewingDisplay,'spd primaries')*(wls(2)-wls(1));
Mviewing_XYZTorgb = inv(Mviewing_rgbToXYZ);

% Render the scene from the start display on the viewing display, to try to
% match XYZ.
[thestartImagergbCalFormat,m,n] = ImageToCalFormat(startImageLinear);
thestartImageXYZCalFormat = Mstart_rgbToXYZ*thestartImagergbCalFormat;
theViewingImagergbCalFormat = Mviewing_XYZTorgb*thestartImageXYZCalFormat;
theViewingImagergb = CalFormatToImage(theViewingImagergbCalFormat,m,n);

% Truncate
% Choice for truncation might be important?
truncate = false;

if truncate
    theViewingImagergbTruncated = theViewingImagergb;
    theViewingImagergbTruncated(theViewingImagergbTruncated < 0) = 0;
    theViewingImagergbTruncated(theViewingImagergbTruncated > 1) = 1;
    outputImageRGB = theViewingImagergbTruncated;

else
    outputImageRGB = theViewingImagergb;
end

end
