%% Generate dataset for learning render matrix
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
theDir = 'all_0.125_true_dataset';
dataBaseDir = getpref(projectName, 'dataDir');
dataFileIn = fullfile(dataBaseDir, thisImageSet, 'image_cifar_all.mat');
dataDirOut = fullfile(dataBaseDir, thisImageSet, theDir);
if (~exist(dataDirOut, 'dir'))
    mkdir(dataDirOut);
end

load(dataFileIn);
displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 0.125, 'display', display.CRT12BitDisplay);

testImage = reshape(image_all(1, :), [32, 32, 3]);
[~, ~, testLinearImage, testConeVec] = retina.compute(testImage);

nImage = 5e3;
allConeVec = zeros(nImage, length(testConeVec));
allLinearImage = zeros(nImage, length(testLinearImage(:)));

parfor idx = 1:nImage
    
    inputImage = reshape(image_all(idx, :), [32, 32, 3]);   
    [~, ~, linearImage, coneExcitation, ~, ~, ~] = retina.compute(inputImage);
    
    allConeVec(idx, :) = coneExcitation;
    allLinearImage(idx, :) = linearImage(:);
    
end

%% Render matrix approximation
regEstimator = RegressionEstimator(allLinearImage, allConeVec);

%% Load linear image dataset
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');
basisInDir = fullfile(dataBaseDir, 'CIFAR_all');
load(fullfile(basisInDir, 'image_cifar_all_linear.mat'));

%% Render matrix test
nSample = 36;
for idx = 1:nSample
    subplot(6, 6, idx);
    testIdx =randi([1, size(image_all, 1)]);    
    inputImage = reshape(image_all(testIdx, :), [32, 32, 3]);
    
    [~, ~, linearImage, gt, ~, ~, ~] = retina.compute(inputImage);
    recon = regEstimator.estimate(reshape(linearImage, [1, 32*32*3]));
    
    scatter(gt, recon); grid on; hold on;
    refPoint = [-500, 4000];
    plot(refPoint, refPoint);    
    axis square;
    xlim(refPoint);
    ylim(refPoint);       
end
suptitle('Render Matrix Approximation');