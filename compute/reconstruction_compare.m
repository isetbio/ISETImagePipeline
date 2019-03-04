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
nTest  = size(imageTe, 1) * 0.5;
nDiag  = min(size(imageTr));

%% Estimate the render matrix
% We want to use a linear approximation of the isetbio routine (which is
% linear indeed up to cone isomerization. The next few sections evaluate
% how good the linear estimation is (should be perfect).

regEstimator = RegressionEstimator(imageTr, coneVecTr);

%% Learn PCA basis
% We want to learn a PCA basis for image reconstruction, which should be
% equivalent to Gaussian markov image model, or 1/f statistics in the
% frequency domain.

[regBasis, mu] = computeBasisPCA(imageTr, 32);

%% Gaussian Prior - Poisson Likelihood Estimator
renderMatrix = regEstimator.W';
estimator = PoissonGaussianEstimator(renderMatrix, regBasis, mu');
estimator.setRegPara(2e3);
estimator.dispOff();

evalObj = CrossValidation(coneVecTe, imageTe, nTest);
[totalMSE, listMSE] = evalObj.evalTest(obj, estimator);

%% Gaussian Prior - Gaussian Likelihood Estimator (Ridge Regression)
renderMatrix = regEstimator.W';
estimator = RidgeGaussianEstimator(renderMatrix, regBasis, mu');

%% Cross-validation
evalObj = CrossValidation(coneVecTe, imageTe, nTest);

stepSize = floor(nDiag / 20);
estimator.setParaList(stepSize : stepSize : nDiag);
[~, ~] = evalObj.crossValidate(estimator);

%% Learning and Testing of Regression Estimator
% Add a column of one as feature vector, equivalent with adding bias in regression
coneVecTr = [coneVecTr, ones(nTrain, 1)];
coneVecTe = [coneVecTe, ones(nTest,  1)];

%% Learning the regression estimator
nDiag = min(size(coneVecTr));
regEstimator = RegressionEstimator(coneVecTr, imageTr);

%% Cross validation
evalObj = CrossValidation(coneVecTe, imageTe, nTest);

stepSize = floor(nDiag / 20);
regEstimator.setParaList(stepSize : stepSize : nDiag);
[paraList, mse] = evalObj.crossValidate(regEstimator);
