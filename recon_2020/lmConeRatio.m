% Effect of L/M Cone Proportion on Image Reconstruction

%% Define constant
imageSize = [64, 64, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');
ratio   = [0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0];

%% Generate a cone mosaic
% Analysis with normal optics
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display, 'pupilSize', 2.5);

%% Load images
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

%% Manipulate the number of cones
% and generate the corresponding render matrix
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    % Change the cone ratio of the mosaic using "reassign_" function
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

%% Run reconstruction
regPara = 1e-3;
output = computeRecon(input, renderArray, prior, regPara, imageSize);

%% Show results
plotResults(input, output, ratio, display, imageSize);

%% Helper function
% Compute image reconstruction for each mosaic
% See "imageRecon.mlx" for basics of reconstruction
function output = computeRecon(input, renderArray, prior, regPara, imageSize)
nImage = size(input, 1);
output = zeros([length(renderArray), nImage, imageSize]);
parfor i = 1:length(renderArray)
    render = renderArray{i};
    estimator = ...
        PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regPara, 4, imageSize);
    
    for j = 1:nImage
        % compute cone response to the images
        image = input(j, :, :, :);
        resp  = poissrnd(render * image(:));
        
        nIter = 500;
        recon = estimator.estimate(resp, nIter, rand([prod(imageSize), 1]), true, 1.0, 'final');
        output(i, j, :, :, :) = recon;
    end
end
end

% Compute reconstruct error and show reconstructed images
function plotResults(input, output, ratio, display, imageSize)

% Compute RMSE
nImage = size(input, 1);
rmse = zeros([length(ratio), nImage]);
for i = 1:length(ratio)
    for j = 1:nImage
        inputImage  = input(j, :, :, :);
        outputImage = output(i, j, :, :, :);
        rmse(i, j) = norm(inputImage(:) - outputImage(:));
    end
end

% Plot RMSE
figure();
errorbar(ratio, mean(rmse, 2), std(rmse, 0, 2) / sqrt(nImage), '--ok', 'LineWidth', 2);

xticks(ratio);
grid off; box off; hold on;

% Show original images
figure();
plotAxis = tight_subplot(length(ratio) + 1, nImage, [.01 .01], [.01 .01], [.01 .01]);
for idx = 1:nImage
    
    image = reshape(input(idx, :, :, :), imageSize);
    image = gammaCorrection(image, display);
    axes(plotAxis(idx));
    imshow(image, 'InitialMagnification', 200);
end

% Show reconstructed images
for i = 1:length(ratio)
    for j = 1:nImage
        image = reshape(output(i, j, :, :, :), imageSize);
        image = gammaCorrection(image, display);
        axes(plotAxis(i * nImage + j));
        imshow(image, 'InitialMagnification', 200);
    end
end
end
