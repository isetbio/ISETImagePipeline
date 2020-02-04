%% Load retina, render matrix, display setup
load('./retinaRenderM.mat');

projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

retina.Mosaic.integrationTime = retina.Mosaic.integrationTime * 5;
render = render * 5;

%% Sample test image
thisImageSet = 'ILSVRC';
imageName    = 'ILSVRC2017_test_00000021.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.5);

patch = sampleImage(image, imageSize(1));
patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 400);

%% Generate Cone Vec
[~, ~, patchLinear, patchConeVec] = retina.compute(patch);
showPlot = true;
if showPlot        
    figure();
    imshow(patch, 'InitialMagnification', 400);
    retina.visualizeOI();
    retina.visualizeExcitation();
end

figure();
scatter(patchConeVec, render * patchLinear(:));
grid on; hold on;
refPoint = [0, 15000];
plot(refPoint, refPoint);
axis square;
xlim(refPoint);
ylim(refPoint);

%% Reconstruction with Poisson likelihood and sparse prior
estimator = PoissonSparseEstimator(render, inv(regBasis), MU', 0.05, 4, imageSize);

meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

reconImage = estimator.estimate(response, 1.5e3, rand([prod(imageSize), 1]), true);

%% Results of reconstruction
figure();

subplot(1, 2, 1);
imshow(patch, 'InitialMagnification', 400); 
title('input');

subplot(1, 2, 2);
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 400); 
title('reconstruction');

%% Auto generation of unpacked montage for validation
retina.reconValidation(patch, reconImage, imageSize(1), response, estimator);

%% Cross-validation

