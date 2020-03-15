%% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

%% Load display setup & constant
imageSize = [128, 128, 3];

thisImageSet = 'ILSVRC';
imageName    = 'ILSVRC2017_test_00000021.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.4);

patch = sampleImage(image, imageSize(1));
patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 400);

%% Combine blur into render matrix
eccX = 10; eccY = 10;
retinaBlur = BlurConeResponse(eccX, eccY, 5, 'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

%% Visualization
[~, ~, linear, coneVec] = retinaBlur.compute(patch);

retinaBlur.visualizeOI();
retinaBlur.visualizeMosaic();
retinaBlur.visualizeExcitation();

%% Render matrix
renderBlur = retinaBlur.forwardRender(imageSize);

%% Show results
for idx = 1:size(inputLinear, 1)
    image = reshape(inputLinear(idx, :, :, :), imageSize);
    
    fltSD = 5;
    image = imgaussfilt(image, fltSD);
    image = invGammaCorrection(image, display.CRT12BitDisplay);
    
    figure();
    subplot(1, 2, 1);
    imshow(image, 'InitialMagnification', 400);
    
    subplot(1, 2, 2);
    output = reshape(outputImage(idx, :, :, :), imageSize);
    imshow(invGammaCorrection(output, display.CRT12BitDisplay), 'InitialMagnification', 400);
end
