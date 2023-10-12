function [outputImageRGB,outputImagergb] = RGBRenderAcrossDisplays(inputImageRGB, startDisplay, viewingDisplay, varargin)
% Correct image from one display for viewing on another
%
% Synopsis:
%    [outputImageRGB,outputImagergb] = rgb2aoDisplay(inputImageRGB, startDisplay, viewingDisplay)
%
% Description:
%    Input a gamma corrected RGG image with respect to passed startDisplay.
%    Generate a displayed image as close to metameric as possible to 
%    for the passed viewingDisplay.
%
%    Typically we'd be doing this with our mono display as start display
%    and our conventional display as the viewing display.
%
%    The tutorial t_renderMonoDisplayImage unpacks the calculation done
%    here, and also verifies that this routine returns the same answer as
%    the unpacked calculation.
%
% Inputs:
%    inputImageRGB  - Gamma corrected input image for startDisplay
%    startDisplay   - ISETBio display structure corresponding to the input image 
%    viewingDispla  - ISETBio display structure corresponding to the output image
%
% Outputs:
%    outputImageRGB - Gamma corrected output image
%    outputImagergb - Linear output image
%
% Optional key/value pairs
%    'linearInput'               - Input linear rather than gamma corrected
%                                  image. Default false.
%    'viewingDisplayScaleFactor' - Multiply viewing display primaries by
%                                  this.  Default 1.
%    'wls'                       - Spectral wavelength sampling, column
%                                  vector.  Default (400:10:700)'
%    'verbose'                   - Print out diagnostic info. Default false.
%
% See also: t_renderMonoDisplayImage, aoStimRecon, aoStimReconRunMany

% History:
%   04/04/23  chr  Made callable function from t_renderMonoDisplayImage
%   10/12/23  dhb  Code review. Fixing comments to indicate input is gamma
%                  corrected RGB
%             dhb  Return gamma corrected RGB as main return.  Second
%                  return is linear rgb.  (Original version returned linear
%                  rgb but with a name tha suggested gamma corrected.)
      
% Parse key value pairs
p = inputParser;
p.addParameter('linearInput', false, @islogical);
p.addParameter('viewingDisplayScaleFactor',1,@isnumeric);
p.addParameter('wls',(400:10:700)',@isnumeric);
p.addParameter('verbose',false,@islogical);
parse(p, varargin{:});

% Convenience assumption that all displays operate on the same wave structure
wls = p.Results.wls;
startDisplay = displaySet(startDisplay,'wave',wls);
viewingDisplay = displaySet(viewingDisplay,'wave',wls);

% Convert input to gamma corrected RGB?
if (p.Results.linearInput)
    inputImageRGB = gammaCorrection(inputImageRGB,forwardDisplay);
end

% Scale recon display primaries to try to keep things in range
viewingDisplay = displaySet(viewingDisplay,'spd primaries',displayGet(viewingDisplay,'spd primaries')*p.Results.viewingDisplayScaleFactor);

meanLuminanceCdPerM2 = [];
[~, ~, startImageLinear] = sceneFromFile(inputImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, startDisplay);
if (p.Results.verbose)
    minr = min(min(startImageLinear(:,:,1)));
    ming = min(min(startImageLinear(:,:,2)));
    minb = min(min(startImageLinear(:,:,3)));
    maxr = max(max(startImageLinear(:,:,1)));
    maxg = max(max(startImageLinear(:,:,2)));
    maxb = max(max(startImageLinear(:,:,3)));
    fprintf('\trgb2aoDisplay: min input linear: %0.2f, %0.2f, %0.2f\n',minr,ming,minb);
    fprintf('\trgb2aoDisplay: max input linear: %0.2f, %0.2f, %0.2f\n',maxr,maxg,maxb);
end

% Get information we need to render scenes from their spectra through
% the viewing display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mstart_rgbToXYZ = T_XYZ*displayGet(startDisplay,'spd primaries')*(p.Results.wls(2)-p.Results.wls(1));
Mstart_XYZTorgb = inv(Mstart_rgbToXYZ);
Mviewing_rgbToXYZ = T_XYZ*displayGet(viewingDisplay,'spd primaries')*(p.Results.wls(2)-p.Results.wls(1));
Mviewing_XYZTorgb = inv(Mviewing_rgbToXYZ);

% Render the scene from the start display on the viewing display, try to match XYZ.
[thestartImagergbCalFormat,m,n] = ImageToCalFormat(startImageLinear);
thestartImageXYZCalFormat = Mstart_rgbToXYZ*thestartImagergbCalFormat;
theViewingImagergbCalFormat = Mviewing_XYZTorgb*thestartImageXYZCalFormat;
theViewingImagergb = CalFormatToImage(theViewingImagergbCalFormat,m,n);
if (p.Results.verbose)
    minr = min(min(theViewingImagergb(:,:,1)));
    ming = min(min(theViewingImagergb(:,:,2)));
    minb = min(min(theViewingImagergb(:,:,3)));
    maxr = max(max(theViewingImagergb(:,:,1)));
    maxg = max(max(theViewingImagergb(:,:,2)));
    maxb = max(max(theViewingImagergb(:,:,3)));
    fprintf('\trgb2aoDisplay: min output linear: %0.2f, %0.2f, %0.2f\n',minr,ming,minb);
    fprintf('\trgb2aoDisplay: max output linear: %0.2f, %0.2f, %0.2f\n',maxr,maxg,maxb);
end

% Truncate.  We want the output image in gamut.
% 
% Truncate linear rgb into physical range. Max
% may still be larger than 1, though.
truncate = true;
if truncate
    theViewingImagergbTruncated = theViewingImagergb;
    theViewingImagergbTruncated(theViewingImagergbTruncated < 0) = 0;
    outputImagergb = theViewingImagergbTruncated;
else
    outputImagergb = theViewingImagergb;
end

% Gamma correct output 
outputImageRGB = gammaCorrection(outputImagergb, viewingDisplay);

end
