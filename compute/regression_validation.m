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
imageSize = [32, 32, 3];

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Add a column of one as feature vector, equivalent with adding bias in regression
coneVecTr = [coneVecTr, ones(nTrain, 1)];
coneVecTe = [coneVecTe, ones(nTest,  1)];

%% Learning the regression estimator
nDiag = min(size(coneVecTr));
regEstimator = RegressionEstimator(coneVecTr, imageTr, 'nDiag', round(nDiag * 0.5));

%% Simple testing
evalObj = CrossValidation(coneVecTe, imageTe, nTest);
for idx = 1:8
    [recon, test] = evalObj.sampleTest(regEstimator, false);
    rgbTest  = invGammaCorrection(reshape(test, imageSize), display.CRT12BitDisplay);
    rgbRecon = invGammaCorrection(reshape(recon, imageSize), display.CRT12BitDisplay);
                
    subplot(4, 4, idx * 2 - 1);
    imshow(rgbTest, 'InitialMagnification', 500);
    subplot(4, 4, idx * 2);    
    imshow(rgbRecon, 'InitialMagnification', 500);    
end

%% Cross validation
evalObj = CrossValidation(coneVecTe, imageTe, nTest);

regEstimator.setParaList(50 : 200 : 2050);
[paraList, mse] = evalObj.crossValidate(regEstimator);

%% Plot hyperparameter - loss
figure();
plot(paraList, mse, '-o', 'LineWidth', 2);

%% Set parameter (use the best cross-validated hyperparameter)
regEstimator.setRegPara(850);
