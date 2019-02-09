%% Setup training and test dataset
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
dataBaseDir = getpref(projectName, 'dataDir');

coneVecDir = 'all_1_true_coneVec';
coneVecTr  = 'coneVector25k.mat';
coneVecTe  = 'coneVector_hold_5k.mat';

ConeVecIn = fullfile(dataBaseDir, thisImageSet, coneVecDir);
D = struct2cell(load(fullfile(ConeVecIn, coneVecTr)));
coneVecTr = D{1}; clear D;
D = struct2cell(load(fullfile(ConeVecIn, coneVecTe)));
coneVecTe = D{1}; clear D;

imageFileIn = fullfile(dataBaseDir, thisImageSet);
imageTr = 'cifar_25k.mat';
imageTe = 'cifar_hold_5k.mat';

D = struct2cell(load(fullfile(imageFileIn, imageTr)));
imageTr = D{1}; clear D;
D = struct2cell(load(fullfile(imageFileIn, imageTe)));
imageTe = D{1}; clear D;

nTrain = size(imageTr, 1);
nTest  = size(imageTe, 1);

%% Learning the regression estimator
nDiag = min(size(coneVecTr));
regEstimator = RegressionEstimator(imageTr, coneVecTr);

%% Full evaluation
evalObj = CrossValidation(imageTe, coneVecTe, nTest);
[totalMSE, listMSE] = evalObj.evalTest(regEstimator);

