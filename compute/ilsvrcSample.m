%% Load dataset, sample large image patch
projectName = 'ISETImagePipeline';
thisImageSet = 'ILSVRC';
dataBaseDir = getpref(projectName, 'dataDir');

files = dir(fullfile(dataBaseDir, thisImageSet, '*.JPEG'));
dirOut = fullfile(dataBaseDir, 'ILSVRC_mat');

nSample = 6; nPatch = 128; count = 1;
imageData = zeros(nPatch, nPatch, 3, 5500 * nSample);

for file = files'
    image = im2double(imread(fullfile(file.folder, file.name)));
    if(size(image, 3) > 1)
        for idx = 1:nSample
            patch = sampleImage(image, nPatch);
            imageData(:, :, :, count) = patch;
            count = count + 1;
            % imshow(patch, 'InitialMagnification', 400);
        end
    end
end

%% Save file
imageData = imageData(:, :, :, 1:count-1);
save(fullfile(dirOut, 'imageDataLarge.mat'), 'imageData', '-v7.3');

%% Load dataset, sample small image patch and generate linear images
projectName = 'ISETImagePipeline';
thisImageSet = 'ILSVRC';
dataBaseDir = getpref(projectName, 'dataDir');

files = dir(fullfile(dataBaseDir, thisImageSet, '*.JPEG'));

nSample = 20; nPatch = 16; count = 1;
imageData = zeros(5500 * nSample, nPatch * nPatch * 3);

for file = files'
    image = im2double(imread(fullfile(file.folder, file.name)));
    if(size(image, 3) > 1)
        for idx = 1:nSample
            patch = sampleImage(image, nPatch);
            imageData(count, :) = patch(:);
            count = count + 1;
            % imshow(patch, 'InitialMagnification', 1000);
        end
    end
end

%% Save file
imageData = imageData(1:count-1, :);
save(fullfile(dirOut, 'imageDataSmall.mat'), 'imageData', '-v7.3');

%% Generate dataset for learning render matrix
displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 0.1, 'display', display.CRT12BitDisplay);

imageSize = [16, 16, 3];
testImage = reshape(imageData(1, :), imageSize);
[~, ~, testLinearImage, testConeVec] = retina.compute(testImage);

nImage = size(imageData, 1);
imageDataLinear = zeros(nImage, length(testLinearImage(:)));

parfor idx = 1:nImage    
    inputImage = reshape(imageData(idx, :), imageSize);
    [~, ~, linearImage, coneExcitation, ~, ~, ~] = retina.compute(inputImage);
        
    imageDataLinear(idx, :) = linearImage(:);
end
