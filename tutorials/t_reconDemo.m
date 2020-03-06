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
image   = imresize(im2double(imread(fileDir)), 0.4);

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
meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

regPara = [1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1];
reconResult = zeros([length(regPara), imageSize]);

parfor idx = 1:length(regPara)
    
    estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara(idx), 4, imageSize);
    reconResult(idx, :, :, :) = invGammaCorrection(estimator.estimate(response, 1.2e3, rand([prod(imageSize), 1]), true), display.CRT12BitDisplay);        
    
end

%% Results of cross-validation
reconLoss = zeros(1, length(regPara));
figure();
for idx = 1:length(regPara)
    subplot(3, 3, idx);
    reconRGB = reshape(reconResult(idx, :, :, :), imageSize);
    imshow(reconRGB, 'InitialMagnification', 400);
    reconLoss(idx) = sqrt(mean((patch(:) - reconRGB(:)) .^ 2));
end

figure();
plot(log(regPara), reconLoss, '-ok', 'LineWidth', 1);
set(gca,'box','off');
set(gca,'TickDir','out');

xticks(log(regPara));
xticklabels(regPara);

%% Different images on prior - likelihood axis
estimator = PoissonSparseEstimator(render, inv(regBasis), MU', 0.05, 4, imageSize);
figure(1); hold on;
% Optimal
xPos = -estimator.prior(reconImage);
yPos = -estimator.likelihood(response, reconImage(:));
scatter(xPos, yPos, 100, 'k', 'filled');

total = xPos + yPos;

% Original
scatter(-estimator.prior(patchLinear), -estimator.likelihood(response, patchLinear(:)), 100, 'r', 'filled');

%% Null space analysis 
basisMtx = null(render);

%% Show image with null 1
image = 0.25 * rand(imageSize);
[vectRow, vectNull] = linearDecomp(render, basisMtx, image(:));
visualizeImage(vectRow, vectNull, imageSize);

imageNull = reshape(vectNull, imageSize);

addNull = reconImage + imageNull;
figure();
imshow(invGammaCorrection(addNull, display.CRT12BitDisplay), 'InitialMagnification', 400); 

figure(1);
scatter(-estimator.prior(addNull), -estimator.likelihood(response, addNull(:)), 100, 'y', 'filled');

%% Show image with null 2
image = 0.25 * rand(imageSize);
[vectRow, vectNull] = linearDecomp(render, basisMtx, image(:));
visualizeImage(vectRow, vectNull, imageSize);

imageNull = reshape(vectNull, imageSize);

addNull = reconImage - imageNull;
figure();
imshow(invGammaCorrection(addNull, display.CRT12BitDisplay), 'InitialMagnification', 400); 

figure(1);
scatter(-estimator.prior(addNull), -estimator.likelihood(response, addNull(:)), 100, 'y', 'filled');

%% Rand Image 
figure(2);
rndPatch = reconImage + 0.2 * rand(imageSize);
imshow(invGammaCorrection(rndPatch, display.CRT12BitDisplay), 'InitialMagnification', 400);

figure(1);
scatter(-estimator.prior(rndPatch), -estimator.likelihood(response,rndPatch(:)), 100, 'b', 'filled');

%% BW image
figure();
grayPatch = rgb2gray(patch);
imshow(grayPatch);

grayImage = zeros(imageSize);
grayImage(:, :, 1) = grayPatch;
grayImage(:, :, 2) = grayPatch;
grayImage(:, :, 3) = grayPatch;
imshow(grayImage, 'InitialMagnification', 400);

figure(1);
[~, ~, grayLinear, ~] = retina.compute(grayImage);
scatter(-estimator.prior(grayLinear), -estimator.likelihood(response, grayLinear(:)), 100, 'b', 'filled');

%% No Prior Recon
estimator = PoissonSparseEstimator(render, inv(regBasis), MU', 0.0, 4, imageSize);
reconImageNoise = estimator.estimate(response, 1.5e3, rand([prod(imageSize), 1]), true);

%% Add image to the plot
figure();
imshow(invGammaCorrection(reconImageNoise, display.CRT12BitDisplay), 'InitialMagnification', 400);

figure(1); hold on;
scatter(-estimator.prior(reconImageNoise), -estimator.likelihood(response,reconImageNoise(:)), 100, 'g', 'filled');

xaxisLim = xlim();
yaxisLim = ylim();
plot([xaxisLim(1), xaxisLim(2)], [total - xaxisLim(1), total - xaxisLim(2)], '-k', 'LineWidth', 1.5);
plot([xPos, xPos], [yaxisLim(1), yPos], '--k', 'LineWidth', 1.5);
plot([xaxisLim(1), xPos], [yPos, yPos], '--k', 'LineWidth', 1.5);

set(gca, 'TickDir', 'out');
set(gca, 'box', 'off');

%% Helper function
function [vectRow, vectNull] = linearDecomp(render, nullBasis, imageVec)
% Projection onto Null Space
coeffNull = nullBasis' * imageVec;
vectNull  = nullBasis * coeffNull;

% Orthogonal Complement
vectRow = imageVec - vectNull;

% Check
fprintf('Vector norm for row image: %.4f \n', norm(render * vectRow));
fprintf('Vector norm for null image: %.4f \n', norm(render * vectNull));
fprintf('Vector dot product: %.4f \n', vectRow' * vectNull);

end

function visualizeImage(vectRow, vectNull, imageSize)
figure();
subplot(1, 3, 1);
imshow(reshape(vectRow + vectNull, imageSize), 'InitialMagnification', 500);
title('Image');

subplot(1, 3, 2);
imshow(reshape(vectRow, imageSize), 'InitialMagnification', 500);
title('Row Space');

subplot(1, 3, 3);
imshow(reshape(vectNull, imageSize), 'InitialMagnification', 500);
title('Null Space');
end
