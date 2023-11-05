% Construct Sparse-Coding Based Prior of Natural Color Images

%% Sample from the ILSVRC image dataset
% Path to the dataset
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
files = dir(fullfile(dataBaseDir, 'ILSVRC_train', '*.JPEG'));

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
display = displayCreate('CRT12BitDisplay');
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 0.1, 'display', display);

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

%% Scale the basis to normalize variance
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

%% Fit exp distribution on the coefficients
allProj = projs(:);
index = randsample(length(allProj), 2.5e4);

xmax = 12;
figure(); subplot(1, 2, 1);
histogram(allProj(index), 100, 'Normalization', 'pdf');
xlim([-xmax, xmax]);
set(gca, 'yscale', 'log'); box off;
set(gca,'TickDir','out');

subplot(1, 2, 2);
histogram(abs(allProj(index)), 100, 'Normalization', 'pdf');
set(gca, 'yscale', 'log'); box off;

xlim([0, xmax]);
ylimits = ylim();

yyaxis right
dist = fitdist(abs(allProj), 'Exponential');
plot(0:0.1:xmax, dist.pdf(0:0.1:xmax), 'LineWidth', 2);
set(gca, 'yscale', 'log'); box off;
set(gca,'TickDir','out');
xlim([0, xmax]);
ylim(ylimits);
