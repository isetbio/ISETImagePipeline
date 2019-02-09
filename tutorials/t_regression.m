%% Setup training and test dataset
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
dataBaseDir = getpref(projectName, 'dataDir');

coneVecDir = 'all_1_true_coneVec';
coneVecTr  = 'coneVector10k.mat';
coneVecTe  = 'coneVector_hold.mat';

ConeVecIn = fullfile(dataBaseDir, thisImageSet, coneVecDir);
D = struct2cell(load(fullfile(ConeVecIn, coneVecTr)));
coneVecTr = D{1}; clear D;
D = struct2cell(load(fullfile(ConeVecIn, coneVecTe)));
coneVecTe = D{1}; clear D;

imageFileIn = fullfile(dataBaseDir, thisImageSet);
imageTr = 'cifar_10k.mat';
imageTe = 'cifar_hold_5k.mat';

D = struct2cell(load(fullfile(imageFileIn, imageTr)));
imageTr = D{1}; clear D;
D = struct2cell(load(fullfile(imageFileIn, imageTe)));
imageTe = D{1}; clear D;

nTrain = size(imageTr, 1);
nTest  = size(imageTe, 1);

%% Learning the regression estimator
nDiag = min(size(coneVecTr));
regEstimator = RegressionEstimator(coneVecTr, imageTr, 'nDiag', round(nDiag * 0.5));

%% Simple testing
evalObj = CrossValidation(coneVecTe, imageTe, nTest);
for nSample = 1:20
    recon = evalObj.sampleTest(regEstimator);
    pause(0.5);
end

%% Full evaluation
[totalMSE, listMSE] = evalObj.evalTest(regEstimator);
histogram(listMSE); grid on;

%% Cross validation
stepSize = round(nDiag / 10);
regEstimator.setParaList(stepSize : stepSize : nDiag);
[paraList, mse] = evalObj.crossValidate(regEstimator);
