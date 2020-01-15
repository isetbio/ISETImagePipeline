function [cielab, cielabSpatial] = labDistance(image1, image2)

% Specify image size in degrees
imageSizeDeg = 1.0;
SCIELAB_LAB_Ver = '2000';

% Convert each image to cal format.
%
% This strings out the pixels as a 3 by nPixels matrix.
% It's easier to apply matrix operations in this format.
%
% Assume they all have same image size.
[image1_Cal,mPixels1,nPixels1] = ImageToCalFormat(image1);
[image2_Cal,mPixels2,nPixels2] = ImageToCalFormat(image2);
if (mPixels1 ~= mPixels2 || nPixels1 ~= nPixels2)
    error('Images being compared do not have same size');
end
if (mPixels1 ~= nPixels1)
    error('Images are not square');
end

% Convert RGB to XYZ.
%
% Assume each image is sRGB, which is probably close enough.
% Could match this to the assumptions we used in ISETBio, but
% let's not worry about that now.
image1XYZ_Cal = SRGBPrimaryToXYZ(SRGBGammaUncorrect(255*image1_Cal));
image2XYZ_Cal = SRGBPrimaryToXYZ(SRGBGammaUncorrect(255*image2_Cal));

% Get "white point" for LAB conversion
whiteXYZ = SRGBPrimaryToXYZ([1 1 1]');

% Convert images to LAB
%
% This is 1976 CIELAB
image1Lab_Cal = XYZToLab(image1XYZ_Cal,whiteXYZ);
image2Lab_Cal = XYZToLab(image2XYZ_Cal,whiteXYZ);

% Get CIELAB DE difference at each pixel
CIELAB_DEImage_Cal = image2Lab_Cal-image1Lab_Cal;

% A better metric is to use S-CIELAB
%
% Convert XYZ back to images
image1XYZ = CalFormatToImage(image1XYZ_Cal,mPixels1,nPixels1);
image2XYZ = CalFormatToImage(image2XYZ_Cal,mPixels2,nPixels2);
whiteXYZCell{1} = whiteXYZ;
whiteXYZCell{2} = whiteXYZ;

% Compute pixels per degree
pixelsPerDeg = mPixels1/imageSizeDeg;

% Set parameters
params.deltaEversion = SCIELAB_LAB_Ver;
params.sampPerDeg  = pixelsPerDeg;
params.imageFormat = 'xyz2';
params.filterSize  = pixelsPerDeg;
params.filters = [];

SCIELAB_DEImage = scielab(image1XYZ,image2XYZ,whiteXYZCell,params);
SCIELAB_DEs = ImageToCalFormat(SCIELAB_DEImage);

% Convert to image mean value
nPixels = size(CIELAB_DEImage_Cal,2);
CIELAB_DEs = zeros(1,nPixels);
for ii = 1:nPixels
    thePixel = CIELAB_DEImage_Cal(:,ii);
    thePixelDE = norm(thePixel);
    CIELAB_DEs(ii) = thePixelDE;
end

cielab = mean(CIELAB_DEs);
cielabSpatial = mean(SCIELAB_DEs);

end

