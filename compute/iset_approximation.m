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
nDiag  = min(size(imageTr));

%% Tricks
imageTr = [imageTr, ones(25000, 1)];
imageTe = [imageTe, ones(5000, 1)];

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
    xlabel('True Cone Mean Response');
    ylabel('Estimated Cone Mean Response');    
    axis square;    
    xlim(refPoint);
    ylim(refPoint);
end