% Construct Sparse-Coding Based Prior of Natural Color Images

%% Sample from the ILSVRC image dataset
% Path to the dataset
dropbox = false;
if dropbox
    dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
    files = dir(fullfile(dataBaseDir, 'ILSVRC_train', '*.JPEG'));
else
    % Special case for the server
    dataBase = '/home/lingqi/Data/denoiser_recon/utils/dataset/islvrc/train';
    files = dir(fullfile(dataBase, '*.JPEG'));
end

% There are 5400+ high resolulation natural images
% We will sample 20 small patches of size 16-by-16 from each of them
nImages = 5400; nSample = 20; nPatch = 16;
imageData = zeros(5400 * nSample, nPatch * nPatch * 3);

% Read and sample each images
% Show each samples or not
showImage = false;
count = 1;
for file = files'
    image = im2double(imread(fullfile(file.folder, file.name)));
    
    if(size(image, 3) == 3)
        for idx = 1:nSample
            patch = sampleImage(image, nPatch);
            imageData(count, :) = patch(:);
            
            count = count + 1;
            if showImage
                imshow(patch, 'InitialMagnification', 1000);
            end
        end
    end
end

imageData = imageData(1:count - 1, :);

%% Generate linear images (pixel values after display gamma function)
% RGB to Linear
% Load specs for a CRT display
display = displayCreate('CRT12BitDisplay');

imageSize = [16, 16, 3];
imageDataLinear = zeros(size(imageData));

% Generate linear images
parfor idx = 1:size(imageData, 1)
    inputImage = reshape(imageData(idx, :), imageSize);    
    [~, ~, linearImage] = sceneFromFile(inputImage, 'rgb', [], display);    
    imageDataLinear(idx, :) = linearImage(:);
end

%% Compute pixel statistics
figure();
subplot(1, 3, 1);
histogram(imageData(:), 40);
box off;
xlabel('Value Before');

subplot(1, 3, 2);
histogram(imageDataLinear(:), 40);
box off;
xlabel('Value After');

imageDataReshape = reshape(imageDataLinear, [size(imageDataLinear, 1), 16, 16, 3]);
averageImage = squeeze(mean(imageDataReshape, 1));

subplot(1, 3, 3);
imshow(gammaCorrection(averageImage, display));

%% Transformation to monochromatic display space
% RGB to AO
displayFile = load('monoDisplay.mat');
monoDisplay = displayFile.monoDisplay;
crtDisplay = displayCreate('CRT12BitDisplay');

% ?resample gamma table
displayGammaBits = 12;
displayGammaGamma = 2;
gammaInput = linspace(0,1,2^displayGammaBits)';
gammaOutput = gammaInput.^displayGammaGamma;
monoDisplay.gamma = gammaOutput(:, [1 1 1]);

% Generate linear images in monochromatic display space
imageDataMono = zeros(size(imageData));
parfor idx = 1:size(imageData, 1)
    inputImage = reshape(imageData(idx, :), imageSize);    
    imageMono = rgb2aoDisplay(inputImage, crtDisplay, monoDisplay);
    imageDataMono(idx, :) = imageMono(:);
end

%% Compute pixel statistics
figure();
subplot(1, 3, 1);
histogram(imageData(:), 40);
box off;
xlabel('Value Before');

subplot(1, 3, 2);
histogram(imageDataMono(:), 40);
box off;
xlabel('Value After');

imageDataReshape = reshape(imageDataMono, [size(imageDataMono, 1), 16, 16, 3]);
averageImage = squeeze(mean(imageDataReshape, 1));

subplot(1, 3, 3);
imshow(averageImage);

%% RICA Analysis to learn basis function
% Whitening, with PCA/SVD

% Choose dataset to build the prior for
% imageSet = imageDataLinear;
imageSet = imageDataMono;

[Z, U, SIG, mu] = whitening(imageSet, 'svd');

% RICA analysis
% RICA: Reconstruction Independent Component Analysis Algorithm
nBasis = 16 * 16 * 3;
result = rica(Z, nBasis, 'IterationLimit', 1e4, 'VerbosityLevel', 1, 'GradientTolerance', 1e-8, 'StepTolerance', 1e-8);
regBasis = U * diag(sqrt(SIG)) * result.TransformWeights;

%% Visualization of basis function
[~] = visualizeBasis(regBasis, 16, size(regBasis, 2), false);

%% Project Images onto the basis function
projBasis = inv(regBasis);
projs = projBasis * (imageData' - mu');

% plot the histogram for the coefficients
nRow = 5; nCol = 5;
figure();
nFig = 1;
for i = 1:nRow
    for j = 1:nCol
        subplot(nRow, nCol, nFig);
        histogram(projs(nFig, :));
        
        nFig = nFig + 1;
        xlim([-8, 8])
        set(gca, 'yscale', 'log'); box off;
    end
end

%% Variance 
projStd = std(projs, 0, 2);
normBasis = projBasis ./ projStd;
projs_norm = normBasis * (imageData' - mu');

% plot the histogram for the coefficients
nRow = 5; nCol = 5;
figure();
nFig = 1;
for i = 1:nRow
    for j = 1:nCol
        subplot(nRow, nCol, nFig);
        histogram(projs_norm(nFig, :));
        
        nFig = nFig + 1;
        xlim([-8, 8])
        set(gca, 'yscale', 'log'); box off;
    end
end

%% Save results
regBasis = inv(normBasis);
[~] = visualizeBasis(regBasis, 16, size(regBasis, 2), false);

save('./monoDisplayPrior.mat', 'mu', 'regBasis', 'imageSize');
