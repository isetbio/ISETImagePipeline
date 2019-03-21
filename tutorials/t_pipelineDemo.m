%% Generate dataset for learning render matrix
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
theDir = 'all_0.125_true_dataset';
dataBaseDir = getpref(projectName, 'dataDir');
dataFileIn = fullfile(dataBaseDir, thisImageSet, 'image_cifar_all.mat');
dataDirOut = fullfile(dataBaseDir, thisImageSet, theDir);
if (~exist(dataDirOut, 'dir'))
    mkdir(dataDirOut);
end

load(dataFileIn);
displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 0.125, 'display', display.CRT12BitDisplay);

imageSize = [32, 32, 3];
testImage = reshape(image_all(1, :), imageSize);
[~, ~, testLinearImage, testConeVec] = retina.compute(testImage);

nImage = 5e3;
allConeVec = zeros(nImage, length(testConeVec));
allLinearImage = zeros(nImage, length(testLinearImage(:)));

parfor idx = 1:nImage
    
    inputImage = reshape(image_all(idx, :), imageSize);   
    [~, ~, linearImage, coneExcitation, ~, ~, ~] = retina.compute(inputImage);
    
    allConeVec(idx, :) = coneExcitation;
    allLinearImage(idx, :) = linearImage(:);
    
end

%% Render matrix approximation
regEstimator = RegressionEstimator(allLinearImage, allConeVec);
renderMatrix = regEstimator.W';

%% Load linear image dataset
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');
basisInDir = fullfile(dataBaseDir, 'CIFAR_all');
load(fullfile(basisInDir, 'image_cifar_all_linear.mat'));

%% Render matrix test
nSample = 36;
for idx = 1:nSample
    subplot(6, 6, idx);
    testIdx =randi([1, size(image_all, 1)]);    
    inputImage = reshape(image_all(testIdx, :), imageSize);
    
    [~, ~, linearImage, gt, ~, ~, ~] = retina.compute(inputImage);
    recon = regEstimator.estimate(reshape(linearImage, [1, prod(imageSize)]));
    
    scatter(gt, recon); grid on; hold on;
    refPoint = [-500, 4000];
    plot(refPoint, refPoint);    
    axis square;
    xlim(refPoint);
    ylim(refPoint);       
end
suptitle('Render Matrix Approximation');

%% Gaussian Prior - Gaussian Likelihood Estimator (Ridge Regression)
% PCA Basis
[regBasis, mu] = computeBasisPCA(allLinearImage(1:9.8e4, :), 32, false);
estimatorRidge = RidgeGaussianEstimator(renderMatrix, regBasis, mu');

%% Laplace Prior - Gaussian Likelihood Estimator (LASSO Regression)
% Sparse Basis
load(fullfile(basisInDir, 'sparse_basis_linear', 'rica_color_3600.mat'));
[Z, U, SIG, MU] = whitening(allLinearImage(1:9.8e4, :), 'svd');

W = Mdl.TransformWeights;
regBasis = U * diag(sqrt(SIG)) * W;
estimatorLasso = LassoGaussianEstimator(renderMatrix, regBasis, MU');

%% Visualization
[~] = visualizeBasis(estimatorRidge.Basis, 32, size(estimatorRidge.Basis, 2), false);
[~] = visualizeBasis(estimatorLasso.Basis, 32, size(estimatorLasso.Basis, 2), false);

retina.visualizeMosaic();

%% Set some parameters for the estimation
estimatorRidge.setLambda(1);

estimatorLasso.dispOff();
estimatorLasso.setLambda(1);
estimatorLasso.setTolerance(1e-8);
estimatorLasso.setIterationLimit(5e4);

%% Simple evaluation
testIdx   = randi([9e4, size(image_all, 1)]);
testImage = reshape(image_all(testIdx, :), imageSize);
[~, ~, ~, coneRes] = retina.compute(testImage);
retina.visualizeOI();
retina.visualizeExcitation();

reconRidge = estimatorRidge.estimate(coneRes');
reconRidge = invGammaCorrection(reshape(reconRidge, imageSize), display.CRT12BitDisplay);

reconLasso = estimatorLasso.estimate(coneRes');
reconLasso = invGammaCorrection(reshape(reconLasso, imageSize), display.CRT12BitDisplay);

figure;
subplot(1, 3, 1);
imshow(testImage, 'InitialMagnification', 500);
title('Original');

subplot(1, 3, 2);
imshow(reconRidge, 'InitialMagnification', 500);
title('Gaussian Prior');

subplot(1, 3, 3);
imshow(reconLasso, 'InitialMagnification', 500);
title('Lasso Prior');
