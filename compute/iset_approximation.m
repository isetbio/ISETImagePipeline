%% Setup training and test dataset
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
dataBaseDir = getpref(projectName, 'dataDir');
dataInDir = fullfile(dataBaseDir, thisImageSet, 'all_1_true_dataset');

coneVecTr  = 'coneVector10k.mat';
coneVecTe  = 'coneVector10k.mat';

D = struct2cell(load(fullfile(dataInDir, coneVecTr)));
coneVecTr = D{1}; coneVecTr = coneVecTr(1:10000, :); clear D;
D = struct2cell(load(fullfile(dataInDir, coneVecTe)));
coneVecTe = D{1}; coneVecTe = coneVecTe(9001:1e4, :); clear D;

imageTr = 'linearImage10k.mat';
imageTe = 'linearImage10k.mat';

D = struct2cell(load(fullfile(dataInDir, imageTr)));
imageTr = D{1}; imageTr = imageTr(1:10000, :); clear D;
D = struct2cell(load(fullfile(dataInDir, imageTe)));
imageTe = D{1}; imageTe = imageTe(9001:1e4, :); clear D;

nTrain = size(imageTr, 1);
nTest  = size(imageTe, 1);
nDiag  = min(size(imageTr));

%% Tricks
imageTr = [imageTr, ones(nTrain, 1)];
imageTe = [imageTe, ones(nTest, 1)];

%% Learning the regression estimator
regEstimator = RegressionEstimator(imageTr, coneVecTr);

%% Simple evaluation
evalObj = CrossValidation(imageTe, coneVecTe, nTest);
[recon, gt] = evalObj.sampleTest(regEstimator, false);

figure;
scatter(gt, recon); grid on;
xlabel('True Cone Mean Response');
ylabel('Estimated Cone Mean Response');
axis equal;
axis square;

%% Full evaluation
[totalMSE, listMSE] = evalObj.evalTest(regEstimator);

%% Cross validation
stepSize = round(nDiag / 20);
regEstimator.setParaList(stepSize : stepSize : nDiag);
[paraList, mse] = evalObj.crossValidate(regEstimator);

%% Reset regularization parameter
regEstimator.setRegPara(600);

%% Simple evaluation
nSample = 36;
for idx = 1:nSample
    subplot(6, 6, idx);
    [recon, gt] = evalObj.sampleTest(regEstimator, false);
    
    scatter(gt, recon); grid on; hold on;
    refPoint = [-500, 5000];
    plot(refPoint, refPoint);
    
    axis square;    
    xlim(refPoint);
    ylim(refPoint);
end

%% Learning the normalized regression estimator
nrmregEstimator = NrmRegressionEstimator(imageTr, coneVecTr);

%% Simple evaluation
evalObj = CrossValidation(imageTe, coneVecTe, nTest);

nSample = 36;
for idx = 1:nSample
    subplot(6, 6, idx);
    [recon, gt] = evalObj.sampleTest(nrmregEstimator, false);

    scatter(gt, recon); grid on; hold on;
    refPoint = [-500, 5000];
    plot(refPoint, refPoint);   
    
    axis square;    
    xlim(refPoint);
    ylim(refPoint);
end

%% Cross validation
stepSize = round(nDiag / 20);
nrmregEstimator.setParaList(stepSize : stepSize : nDiag);
evalObj.crossValidate(nrmregEstimator);
