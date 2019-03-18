%% Setup training and test dataset
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
dataBaseDir = getpref(projectName, 'dataDir');
dataInDir = fullfile(dataBaseDir, thisImageSet, 'all_1_true_dataset');

idxTrain = 1:1:1.8e4;
idxTest  = (2e4 - 2e3 + 1):1:2e4;

coneVecTr  = 'coneVector20k.mat';
coneVecTe  = 'coneVector20k.mat';

D = struct2cell(load(fullfile(dataInDir, coneVecTr)));
coneVecTr = D{1}; coneVecTr = coneVecTr(idxTrain, :); clear D;
D = struct2cell(load(fullfile(dataInDir, coneVecTe)));
coneVecTe = D{1}; coneVecTe = coneVecTe(idxTest, :); clear D;

imageTr = 'linearImage20k.mat';
imageTe = 'linearImage20k.mat';

D = struct2cell(load(fullfile(dataInDir, imageTr)));
imageTr = D{1}; imageTr = imageTr(idxTrain, :); clear D;
D = struct2cell(load(fullfile(dataInDir, imageTe)));
imageTe = D{1}; imageTe = imageTe(idxTest, :); clear D;

nTrain = size(imageTr, 1);
nTest  = size(imageTe, 1);
nDiag  = min(size(imageTr));

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));
imageSize = [32, 32, 3];

%% Estimate the render matrix
regEstimator = RegressionEstimator(imageTr, coneVecTr);

%% Load learned sparse basis
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');
basisInDir = fullfile(dataBaseDir, 'CIFAR_all');

load(fullfile(basisInDir, 'all_1_true_dataset', 'linearImage100k.mat'));
load(fullfile(basisInDir, 'sparse_basis_linear', 'rica_color_3600.mat'));

[Z, U, SIG, MU] = whitening(allLinearImage(1:9.8e4, :), 'svd'); 
clear allLinearImage;

%% Create the estimator
W = Mdl.TransformWeights;
regBasis = U * diag(sqrt(SIG)) * W;
renderMatrix = regEstimator.W';
estimatorLasso = LassoGaussianEstimator(renderMatrix, regBasis, MU');

[~] = visualizeBasis(estimatorLasso.Basis, 32, size(estimatorLasso.Basis, 2), false);

%% Sparsity check
coff = Z * W;

figure;
nPlot = 8;
for idx = 1:(nPlot * nPlot)
    subplot(nPlot, nPlot, idx);
    coffIdx = randi([1, size(estimatorLasso.Basis, 2)]);
    histogram(coff(:, coffIdx), 'Normalization', 'probability');    
end

figure;
kurtStat = kurtosis(coff);
histogram(kurtStat(kurtStat < 2e2) - 3);

%% Reconstruction test
evalObj = CrossValidation(coneVecTe, imageTe, nTest);
estimatorLasso.dispOn();

figure();
for idx = 1:8
    [recon, test] = evalObj.sampleTest(estimatorLasso, false);
    rgbTest  = invGammaCorrection(reshape(test, imageSize), display.CRT12BitDisplay);
    rgbRecon = invGammaCorrection(reshape(recon, imageSize), display.CRT12BitDisplay);

    subplot(4, 4, idx * 2 - 1);
    imshow(rgbTest, 'InitialMagnification', 500);
    subplot(4, 4, idx * 2);    
    imshow(rgbRecon, 'InitialMagnification', 500);    
end
