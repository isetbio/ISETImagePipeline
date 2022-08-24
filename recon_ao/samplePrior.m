%% Create an estimator object
% Load prior and display
prior = load('sparsePrior.mat');
display = displayCreate('CRT12BitDisplay');

regPara = 1.0;
stride = 4;
imSize = [64, 64, 3];

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
    end
end

%% Adding noise: Langevin algorithm
nIter = 2500;
stepSize = 250;
nPlot = nIter / stepSize;

figure();

image = rand(prod(imSize), 1);
center = zeros(nIter, 1);
tau = 1e-5;
gamma = 2.5e-4;
for idx = 1:nIter
    % neg prior grad
    [~, grad] = estm.prior(reshape(image, imSize));

    % go uphill on the log-gradient with *added noise*
    % https://en.wikipedia.org/wiki/Metropolis-adjusted_Langevin_algorithm
    image = image + tau * (-grad) + sqrt(2 * gamma) * normrnd(0, 1, prod(imSize), 1);

    % record the value of a pixel
    pxIdx = 2080;
    center(idx) = image(pxIdx);
 
    % plot
    if mod(idx - 1, stepSize) == 0
        plotIdx = floor((idx - 1) / stepSize) + 1;
        subplot(1, nPlot, plotIdx);
        imshow(gammaCorrection(reshape(image, imSize), display));
        title(sprintf('n Iter = %d', idx));
    end
end

% plot ACF of the center pixel
[acf, lags] = autocorr(center, 'NumLags', 500);
figure();
plot(lags, acf, 'LineWidth', 2);
box off;

%% Sampler function
prior = load('sparsePrior.mat');
display = displayCreate('CRT12BitDisplay');

imSize = [100, 100, 3];
nSample = 6;

samples = lgvSampler(prior, nSample, imSize);

figure();
for idx = 1:nSample
    subplot(2, 3, idx);
    imshow(gammaCorrection(reshape(samples(idx, :, :, :), imSize), display));
end
