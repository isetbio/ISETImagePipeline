function converted = lmsVis(rgbImage, type, methodStr)
% RGB -> XYZ
[image_Cal, mPixels, nPixels] = ImageToCalFormat(SRGBGammaCorrect(rgbImage));
imageXYZ_Cal = SRGBPrimaryToXYZ(SRGBGammaUncorrect(image_Cal));

% Cal to images
imageXYZ = CalFormatToImage(imageXYZ_Cal, mPixels, nPixels);

% Convert
whiteXYZ = SRGBPrimaryToXYZ([1 1 1]');

if strcmp(methodStr, 'brettel')
    imageLMS = xyz2lms(imageXYZ, type, methodStr, whiteXYZ);
else
    imageLMS = xyz2lms(imageXYZ, type, methodStr);
end

imageXYZ = lms2xyz(imageLMS);

% Image to Cal
imageXYZ_Cal = ImageToCalFormat(imageXYZ);

% RGB image
rgb = XYZToSRGBPrimary(imageXYZ_Cal);
RGB = SRGBGammaCorrect(rgb);

RGB = CalFormatToImage(RGB, mPixels, nPixels);
converted = RGB ./ 255;