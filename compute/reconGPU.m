%% Run Reconstruction with GPU Array

%% generate human cone mosaic with optics
imageSize = [128, 128, 3];
display = load('display.mat');
prior   = load('sparsePrior.mat');

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'pupilSize', 2.5);

render = retina.forwardRender(imageSize);

%% load images
nImage = 1;
input = zeros([nImage, imageSize]);

figure();
fileType = '.jpeg';
for idx = 1:nImage
    fileName = strcat(num2str(idx), fileType);
    filePath = fullfile('.', 'images', fileName);
    image = imresize(im2double(imread(filePath)), 0.5);
    
    image = sampleImage(image, imageSize(1));
    image = image - min(image(:));
    image = image ./ max(image(:));    
    input(idx, :, :, :) = image;
    
    subplot(2, 4, idx);
    imshow(image, 'InitialMagnification', 500);
end

%% reconstruction without GPU
profile on

regConst = 5e-4; stride = 4;
estimator = ...
    PoissonSparseEstimator(double(render), inv(prior.regBasis), ...
    prior.mu', regConst, stride, imageSize);

output = zeros([nImage, imageSize]);
for idx = 1:nImage
    [~, ~, ~, coneRespVec] = retina.compute(reshape(input(idx, :, :, :), imageSize));
    recon = estimator.runEstimate(coneRespVec, 'maxIter', 100, 'display', 'iter');
    output(idx, :, :, :) = gammaCorrection(recon, display.CRT12BitDisplay);
end

profile viewer

%% show reconstructed images
figure();
for idx = 1:nImage
    subplot(2, 4, idx);
    imshow(reshape(output(idx, :, :, :), imageSize), 'InitialMagnification', 500);
end

%% reconstruction with GPU array
profile on

regConst = 5e-4; stride = 4;

% transfer the render matrix to GPU (gpuArray)
% use single precision for better performance 
% call estimate function with gpu = true
estimator = ...
    PoissonSparseEstimator(gpuArray(single(render)), inv(prior.regBasis), ...
    prior.mu', regConst, stride, imageSize);

output = zeros([nImage, imageSize]);
for idx = 1:nImage
    [~, ~, ~, coneRespVec] = retina.compute(reshape(input(idx, :, :, :), imageSize));
    recon = estimator.runEstimate(coneRespVec, 'maxIter', 100, 'display', 'iter', 'gpu', true);
    output(idx, :, :, :) = gammaCorrection(recon, display.CRT12BitDisplay);
end

profile viewer

%% show reconstructed images
figure();
for idx = 1:nImage
    subplot(2, 4, idx);
    imshow(reshape(output(idx, :, :, :), imageSize), 'InitialMagnification', 500);
end
