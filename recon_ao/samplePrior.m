% Clear and close
clear; close all;

% Load prior and display
prior = load('sparsePrior.mat');
display = displayCreate('CRT12BitDisplay');

regPara = 1.0;
stride = 16;
imSize = [16, 16, 3];

tau = 2e-5;           % 1e-5
gamma = 2.5e-6;       % 2.5e-4
tau = 1e-5;
gamma = 2.5e-4;

% Construct image estimator
estm = PoissonSparseEstimator([], inv(prior.regBasis), ...
                        prior.mu', regPara, stride, imSize);

% We need the log-prior gradident
% term for the calculation
image = rand(imSize);
[~, grad] = estm.prior(image);

%% Correctness check: 
% Gradient ascent starting from a random init
nIter = 250;
stepSize = 50;
nPlot = nIter / stepSize;

figure();
% random init
image = rand(prod(imSize), 1);
tau = 1e-5;
for idx = 1:nIter
    % neg prior grad
    [~, grad] = estm.prior(reshape(image, imSize));

    % go uphill on the log-gradient
    image = image + tau * (-grad);

    % plot
    if mod(idx - 1, stepSize) == 0
        plotIdx = floor((idx - 1) / stepSize) + 1;
        subplot(1, nPlot, plotIdx);
        imshow(gammaCorrection(reshape(image, imSize), display));
        title(sprintf('n Iter = %d', idx));
        drawnow;
    end
end

%% Adding noise: Langevin algorithm
nIter = 4000;
stepSize = 1000;
nPlot = nIter / stepSize;

figure();

image = rand(prod(imSize), 1);
center = zeros(nIter, 1);
for idx = 1:nIter
    % neg prior grad
    [~, grad] = estm.prior(reshape(image, imSize));

    % go uphill on the log-gradient with *added noise*
    % https://en.wikipedia.org/wiki/Metropolis-adjusted_Langevin_algorithm
    image = image + tau * (-grad) + sqrt(2 * gamma) * normrnd(0, 1, prod(imSize), 1);

    % record the value of a pixel
    pxIdx = round(prod(imSize)/2);
    center(idx) = image(pxIdx);
 
    % plot
    if mod(idx - 1, stepSize) == 0
        plotIdx = floor((idx - 1) / stepSize) + 1;
        subplot(1, nPlot, plotIdx);
        imshow(gammaCorrection(reshape(image, imSize), display));
        title(sprintf('n Iter = %d', idx));
        drawnow;
    end
end

% plot ACF of the center pixel
[acf, lags] = autocorr(center, 'NumLags', nIter-1);
figure();
plot(lags, acf, 'LineWidth', 2);
box off;
drawnow;

%% Sampler function
nSample = 3;
sampleSteps = 4000;
samples = lgvSampler(prior, nSample, imSize,'burnIn',sampleSteps,'nStep',sampleSteps, ...
    'stride',stride,'gamma',gamma,'tau',tau);

figure();
for idx = 1:nSample
    subplot(2, nSample, idx);
    imshow(gammaCorrection(reshape(samples(idx, :, :, :), imSize), display));
end

%% Sample by drawing from exponential
patchSize = sqrt(size(prior.regBasis,1)/3);
meanPatch = prior.mu';
for idx = 1:nSample
    rndMu = exprnd(ones(size(prior.mu)))';
    rndSign = randi(2,size(rndMu)); rndSign(rndSign == 2) = -1;
    rndPatch = meanPatch + prior.regBasis*(rndMu .* rndSign);
    max(rndPatch)
    min(rndPatch)
    rndPatch(rndPatch < 0) = 0;
    subplot(2,nSample,idx+nSample); 
    imshow(gammaCorrection(reshape(rndPatch, [patchSize,patchSize,3]), display));
end

