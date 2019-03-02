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
% We want to use a linear approximation of the isetbio routine (which is
% linear indeed up to cone isomerization. The next few sections evaluate
% how good the linear estimation is (should be perfect).

regEstimator = RegressionEstimator(imageTr, coneVecTr);

%% Simple evaluation
evalObj = CrossValidation(imageTe, coneVecTe, nTest);
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
    
    if(idx == 34)
        xlabel('True Cone Mean Response');
    end
    
    if(idx == 13)
        ylabel('Estimated Cone Mean Response');
    end    
end

%% Full evaluation
[totalMSE, listMSE] = evalObj.evalTest(regEstimator);

listMSE = sort(listMSE);
histogram(listMSE(1:floor(nTest * 0.975)));

%% Learn PCA basis
% We want to learn a PCA basis for image reconstruction, which should be
% equivalent to Gaussian markov image model, or 1/f statistics in the
% frequency domain.

[regBasis, mu] = computeBasisPCA(imageTr, 32);

%% Gaussian Prior - Poisson Likelihood estimator
renderMatrix = regEstimator.W';
estimator = PoissonGaussianEstimator(renderMatrix, regBasis, mu');

%% Simple evaluation
figure();
evalObj = CrossValidation(coneVecTe, imageTe, nTest);
for idx = 1:8
    [recon, test] = evalObj.sampleTest(estimator, false);
    rgbTest  = invGammaCorrection(reshape(test, imageSize), display.CRT12BitDisplay);
    rgbRecon = invGammaCorrection(reshape(recon, imageSize), display.CRT12BitDisplay);
    
    subplot(4, 4, idx * 2 - 1);
    imshow(rgbTest, 'InitialMagnification', 500);
    subplot(4, 4, idx * 2);    
    imshow(rgbRecon, 'InitialMagnification', 500);    
end

%% Set hyperparameter
figure();
estimator.setRegPara(round(nDiag * 0.25));
for idx = 1:8
    [recon, test] = evalObj.sampleTest(estimator, false);
    rgbTest  = invGammaCorrection(reshape(test, imageSize), display.CRT12BitDisplay);
    rgbRecon = invGammaCorrection(reshape(recon, imageSize), display.CRT12BitDisplay);
    
    
    subplot(4, 4, idx * 2 - 1);
    imshow(rgbTest, 'InitialMagnification', 500);
    subplot(4, 4, idx * 2);    
    imshow(rgbRecon, 'InitialMagnification', 500);    
end