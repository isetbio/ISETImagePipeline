%% setup
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

%% load matrix
load('retinaRender40.mat');
render = double(render40);

%% show image
load('./inputLinear.mat');
inputLinear = inputLinear([1, 3, 6, 7, 9, 10, 11, 12], :, :, :);

figure();
for idx = 1 : size(inputLinear, 1)
    subplot(2, 4, idx);
    image = reshape(inputLinear(idx, :, :, :), imageSize);
    imshow(gammaCorrection(image, display), 'InitialMagnification', 200);
end

%% run reconstruction
% Build an image reconstruction object with sparse prior
regConst = 0.0; stride = 4;
estimator = ...
    PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regConst, stride, imageSize);

% Run reconstruction on cone response to each images
% reconstructed images are in linear pixel space, need to
% gamma correct them before visulization
nIter = 800; optDisp = 'iter';
output = zeros(size(inputLinear));

parfor idx = 1 : size(inputLinear, 1)
    input = reshape(inputLinear(idx, :, :, :), imageSize);
    coneResp = render * input(:);
    recon = estimator.runEstimate(coneResp, 'maxIter', nIter, 'display', optDisp);
    output(idx, :, :, :) = gammaCorrection(recon, display);
end

%% show results
figure();
for idx = 1 : size(output, 1)
    subplot(2, 4, idx);
    reconImage = reshape(output(idx, :, :, :), imageSize);
    imshow(reconImage, 'InitialMagnification', 200);
end