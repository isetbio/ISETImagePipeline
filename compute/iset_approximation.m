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

%% Simple evaluation
nSample = 36;
for idx = 1:nSample
    subplot(6, 6, idx);
    [recon, gt] = evalObj.sampleTest(regEstimator, false);
    
    scatter(gt, recon); grid on; hold on;
    refPoint = [-500, 4000];
    plot(refPoint, refPoint);
    
    axis square;    
    xlim(refPoint);
    ylim(refPoint);
    xlabel('True Cone Mean Response');
    ylabel('Estimated Cone Mean Response');
end

%% Cross validation
stepSize = round(nDiag / 20);
regEstimator.setParaList(stepSize : stepSize : nDiag);
[paraList, mse] = evalObj.crossValidate(regEstimator);
