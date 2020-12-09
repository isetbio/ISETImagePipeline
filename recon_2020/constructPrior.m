%% Sample from the ILSVRC image dataset
% Path to the dataset
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
files = dir(fullfile(dataBaseDir, 'ILSVRC', '*.JPEG'));

% There are 5500 high resolulation natural images
% We will sample 20 small patches of size 16-by-16 from each of them
nImages = 5500; nSample = 20; nPatch = 16;
imageData = zeros(5500 * nSample, nPatch * nPatch * 3);

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
% Load specs for a CRT display
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 0.1, 'display', display.CRT12BitDisplay);

imageSize = [16, 16, 3];
imageDataLinear = zeros(size(imageData));

% Generate linear images
parfor idx = 1:size(imageData, 1)    
    inputImage = reshape(imageData(idx, :), imageSize);
    [~, ~, linearImage] = retina.compute(inputImage);
        
    imageDataLinear(idx, :) = linearImage(:);
end

%% RICA Analysis to learn basis function
% Whitening, with PCA/SVD
imageSet = imageDataLinear;
[Z, U, SIG, mu] = whitening(imageSet, 'svd');

% RICA analysis
% RICA: Reconstruction Independent Component Analysis
nBasis = 16 * 16 * 3;
result = rica(Z, nBasis, 'IterationLimit', 1e4, 'VerbosityLevel', 1, 'GradientTolerance', 1e-8, 'StepTolerance', 1e-8);
regBasis = U * diag(sqrt(SIG)) * result.TransformWeights;

%% Visualization of basis function
[~] = visualizeBasis(regBasis, 16, size(regBasis, 2), false);
