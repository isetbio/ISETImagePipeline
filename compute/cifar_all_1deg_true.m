projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
theDir = 'all_1_true_coneVec';
dataBaseDir = getpref(projectName, 'dataDir');
dataFileIn = fullfile(dataBaseDir, thisImageSet, 'image_cifar_all.mat');
dataDirOut = fullfile(dataBaseDir, thisImageSet, theDir);
if (~exist(dataDirOut, 'dir'))
    mkdir(dataDirOut);
end

load(dataFileIn);
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true);

testImage = reshape(image_all(1, :), [32, 32, 3]);
[~, ~, testConeVec] = retina.compute(testImage);

nImage = 1e5;
allConeVec = zeros(nImage, length(testConeVec));

parfor idx = 1 : nImage
    
    inputImage = reshape(image_all(idx, :), [32, 32, 3]);    
    [~, ~, coneExcitation] = retina.compute(inputImage);
    
    allConeVec(idx, :) = coneExcitation;
    
end

coneVectorFile = 'coneVector100k.mat';
save(fullfile(dataDirOut, coneVectorFile), 'allConeVec', '-v7.3');