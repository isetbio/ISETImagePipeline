%% Setup training and test dataset
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
dataBaseDir = getpref(projectName, 'dataDir');
dataInDir = fullfile(dataBaseDir, thisImageSet, 'all_1_true_dataset');

coneVecTr  = 'coneVector20k.mat';
coneVecTe  = 'coneVector20k.mat';

idxTrain = 1:1:1.8e4;
idxTest  = (2e4 - 2e3 + 1):1:2e4;

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

%% Learning a sparse basis
% This step could take a really long time, consider using the saved result
% from shared dropbox.

% Whitening, SVD
[Z, U, SIG, MU] = whitening(imageTr, 'svd');

% RICA analysis
nBasis = 3072;
Mdl    = rica(Z, nBasis, 'IterationLimit', 1e3, 'VerbosityLevel', 1, 'GradientTolerance', 1e-4, 'StepTolerance', 1e-4);

% Visualization
W   = Mdl.TransformWeights;
[~] = visualizeBasis(U * diag(sqrt(SIG)) * W, 32, nBasis, false);

regBasis = U * diag(sqrt(SIG)) * W;
estimatorLasso = LassoGaussianEstimator(regEstimator.W', regBasis, MU');
