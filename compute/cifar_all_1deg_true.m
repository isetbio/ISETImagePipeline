projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
theDir = 'all_1_true_dataset';
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
                      'display', display.CRT12BitDisplay);

testImage = reshape(image_all(1, :), [32, 32, 3]);
[~, ~, testLinearImage, testConeVec] = retina.compute(testImage);

nImage = 1e5;
allConeVec = zeros(nImage, length(testConeVec));
allLinearImage = zeros(nImage, length(testLinearImage(:)));

parfor idx = 1:nImage
    
    inputImage = reshape(image_all(idx, :), [32, 32, 3]);    
    [~, ~, linearImage, coneExcitation, ~, ~, ~] = retina.compute(inputImage);
    
    allConeVec(idx, :) = coneExcitation;
    allLinearImage(idx, :) = linearImage(:);
    
end

linearImageFile = 'linearImage100k.mat';
save(fullfile(dataDirOut, linearImageFile), 'allLinearImage', '-v7.3');

coneVectorFile = 'coneVector100k.mat';
save(fullfile(dataDirOut, coneVectorFile), 'allConeVec', '-v7.3');
