%% Clear
clear; close all;

%% Sample images from the sparse prior
% load prior and display
prior = load('sparsePrior.mat');
display = displayCreate('CRT12BitDisplay');

%% Sample a small patch of image
stride = sqrt(size(prior.regBasis, 1) / 3);
imSize = [stride, stride, 3];
meanPatch = prior.mu';

% sample from an exponential distribution 
rndMu = exprnd(0.3980 * ones(size(meanPatch)));

% sample random binary sign
rndSign = rand(size(rndMu));
rndSign(rndSign > 0.5) = 1; rndSign(rndSign < 0.5) = -1;

% multiply the coefficients with the basis function
rndPatch = meanPatch + prior.regBasis * (rndMu .* rndSign);
imshow(gammaCorrection(reshape(rndPatch, imSize), display), ...
        'InitialMagnification', 1000);

%% Sample a larger image
imSize = [100 100 3];
xStep = ceil(imSize(1) / stride);
yStep = ceil(imSize(2) / stride);

fullImg = zeros(xStep * stride, yStep * stride, 3);

for x = 1 : xStep
    for y = 1 : yStep
        fullImg((x - 1) * stride + 1 : x * stride, ...
                (y - 1) * stride + 1: y * stride, :) ...
                = returnSample(prior);
    end
end

image = fullImg(1:imSize(1), 1:1:imSize(2), :);
image = imgaussfilt(image, 0.1);

imshow(gammaCorrection(image, display), ...
    'InitialMagnification', 500);

%% Try climbing uphill from the above
priorStride = 2;
tau = 2e-6;
gamma = tau;
sampleSteps = 200;

sample = lgvSampler(prior, 1, imSize,'burnIn',sampleSteps,'nStep',sampleSteps, ...
    'stride',priorStride,'gamma',gamma,'tau',tau,'seed',image(:));
figure;
imshow(gammaCorrection(reshape(sample, imSize), display));

%% 1/f pink noise
image = zeros(imSize);
for idx = 1:3
    image(:, :, idx) = spectrumSample(imSize(1:2));
end
figure;
imshow(gammaCorrection(reshape(image, imSize), display));

sample = lgvSampler(prior, 1, imSize,'burnIn',sampleSteps,'nStep',sampleSteps, ...
    'stride',priorStride,'gamma',gamma,'tau',tau,'seed',image(:));
figure;
imshow(gammaCorrection(reshape(sample, imSize), display));

%% Helper function that samples from sparse prior
function sample = returnSample(prior)

stride = sqrt(size(prior.regBasis, 1) / 3);
imSize = [stride, stride, 3];

% sample from an exponential distribution 
rndMu = exprnd(0.3980 * ones(size(prior.mu')));

% sample random binary sign
rndSign = rand(size(rndMu));
rndSign(rndSign > 0.5) = 1; rndSign(rndSign < 0.5) = -1;

% multiply the coefficients with the basis function
sample = prior.mu' + prior.regBasis * (rndMu .* rndSign);

% reshape image
sample = reshape(sample, imSize);

end

function image = spectrumSample(imSize)

imCenter = floor(imSize / 2);
amplitude = zeros(imSize);
phase = exp(1i * (rand(imSize) * 2 * pi));

for x = 1:imSize(1)
    for y = 1:imSize(2)
        amplitude(x, y) = ...
            1 / sqrt((x - imCenter(1)) ^ 2 + ...
                     (y - imCenter(2)) ^ 2 + 1);
    end
end

specturm = amplitude .* phase;
image = real(ifft2(ifftshift(specturm)));

image = image - min(image(:));
image = image / max(image(:));

end