%% Run Reconstruction with GPU Array

%% generate human cone mosaic with optics
imageSize = [100, 100, 3];
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
    image = imresize(im2double(imread(filePath)), 0.30);
    
    image = sampleImage(image, imageSize(1));
    image = image - min(image(:));
    image = image ./ max(image(:));    
    input(idx, :, :, :) = image;
    
    subplot(2, 4, idx);
    imshow(image, 'InitialMagnification', 400);
end
set(gcf,'position',[0, 0, 800, 400]);

%% reconstruction without GPU
profile on

regConst = 5e-4; stride = 4;
estimator = ...
    PoissonSparseEstimator(render, inv(prior.regBasis), ...
    prior.mu', regConst, stride, imageSize);

nIter = 500; bounded = true; optDisp = 'final';
output = zeros([nImage, imageSize]);
for idx = 1:nImage
    [~, ~, ~, coneRespVec] = retina.compute(reshape(input(idx, :, :, :), imageSize));
    recon = estimator.estimate(coneRespVec, nIter, rand([prod(imageSize), 1]), bounded, 1.0, optDisp);
    output(idx, :, :, :) = gammaCorrection(recon, display.CRT12BitDisplay);
end

profile viewer

%% reconstruction with GPU array
profile on

regConst = 5e-4; stride = 4;

% transfer the render matrix to GPU (gpuArray)
% call estimate function with gpu = true
estimator = ...
    PoissonSparseEstimator(gpuArray(render), inv(prior.regBasis), ...
    prior.mu', regConst, stride, imageSize);

nIter = 500; bounded = true; optDisp = 'final'; gpu = true;
output = zeros([nImage, imageSize]);
for idx = 1:nImage
    [~, ~, ~, coneRespVec] = retina.compute(reshape(input(idx, :, :, :), imageSize));
    recon = estimator.estimate(coneRespVec, nIter, rand([prod(imageSize), 1]), bounded, 1.0, optDisp, gpu);
    output(idx, :, :, :) = gammaCorrection(recon, display.CRT12BitDisplay);
end

profile viewer

%% show reconstructed images
figure();
for idx = 1:nImage
    subplot(2, 4, idx);
    imshow(reshape(output(idx, :, :, :), imageSize), 'InitialMagnification', 400);
end
set(gcf,'position',[0, 0, 800, 400]);
