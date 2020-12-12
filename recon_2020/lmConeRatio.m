% Effect of L/M Cone Proportion on Image Reconstruction

%% define constant
imageSize = [64, 64, 3];
display = load('display.mat');
prior   = load('sparsePrior.mat');
ratio   = [0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0];

%% generate a cone mosaic
% analysis with normal optics
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display.CRT12BitDisplay, 'pupilSize', 2.5);

%% load images
nImage = 10;
input = zeros([nImage, imageSize]);

fileType = '.jpeg';
for idx = 1:nImage
    fileName = strcat(num2str(idx), fileType);
    filePath = fullfile('.', 'images', fileName);
    image = imresize(im2double(imread(filePath)), 0.25);
    
    image = sampleImage(image, imageSize(1));
    image = image - min(image(:));
    image = image ./ max(image(:));
    
    [~, ~, linearImage] = retina.compute(image);
    input(idx, :, :, :) = linearImage;
end

%% manipulate the number of cones
% and generate the corresponding render matrix
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    if ratio(idx) < 0.95
        retina.reassignLCone(0.0, false);
        retina.reassignMCone(ratio(idx), true);
    else
        retina.reassignLCone(1 - ratio(idx), true);
    end
    
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    renderMtx = retina.forwardRender(imageSize, false, true, false);
    renderArray(idx) = {double(renderMtx)};
end

%% run reconstruction
regPara = 5e-4;
output = computeRecon(input, renderArray, prior, regPara, imageSize);

%% show results
plotResults(input, output, ratio, display, imageSize);

%% helper function
% compute image reconstruction for each mosaic
function output = computeRecon(input, renderArray, prior, regPara, imageSize)
nImage = size(input, 1);
output = zeros([length(renderArray), nImage, imageSize]);
parfor i = 1:length(renderArray)
    render = renderArray{i};
    estimator = ...
        PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regPara, 4, imageSize);
    
    for j = 1:nImage
        image = input(j, :, :, :);
        resp  = render * image(:);
        
        nIter = 500;
        recon = estimator.estimate(resp, nIter, rand([prod(imageSize), 1]), true, 1.0, 'final');
        output(i, j, :, :, :) = recon;
    end
end
end

% compute reconstruct error and show reconstructed images
function plotResults(input, output, ratio, display, imageSize)
display = display.CRT12BitDisplay;

% compute RMSE
nImage = size(input, 1);
rmse = zeros([length(ratio), nImage]);
for i = 1:length(ratio)
    for j = 1:nImage
        inputImage  = input(j, :, :, :);
        outputImage = output(i, j, :, :, :);
        rmse(i, j) = norm(inputImage(:) - outputImage(:));
    end
end

% plot RMSE
figure();
errorbar(ratio, mean(rmse, 2), std(rmse, 0, 2) / sqrt(nImage), '--ok', 'LineWidth', 2);

xticks(ratio);
grid off; box off; hold on;

% show original images
figure();
plotAxis = tight_subplot(length(ratio) + 1, nImage, [.01 .01], [.01 .01], [.01 .01]);
for idx = 1:nImage
    image = reshape(input(idx, :, :, :), imageSize);
    image = gammaCorrection(image, display);
    axes(plotAxis(idx));
    imshow(image, 'InitialMagnification', 200);
end

% show reconstructed images
for i = 1:length(ratio)
    for j = 1:nImage
        image = reshape(output(i, j, :, :, :), imageSize);
        image = gammaCorrection(image, display);
        axes(plotAxis(i * nImage + j));
        imshow(image, 'InitialMagnification', 200);
    end
end
end
