%% load the display setup
display = displayCreate('CRT12BitDisplay');

% Generate cone mosaic - 10,10 deg ecc
eccX = 10; eccY = 10;
retina = ConeResponsePeripheral(eccX, eccY, ...
    'fovealDegree', 1.0, 'display', display, 'pupilSize', 3.0);

retina.visualizeMosaic();

%% load an example image
imageSize = [128, 128, 3];
image = im2double(imread('images/2.jpeg'));
image = imresize(image, imageSize(1) / size(image, 1));

% normalize to insure pixels are in the range of [0, 1]
image = (image - min(image(:)));
image = image ./ max(image(:));

assert(min(image(:)) >= 0);
assert(max(image(:)) <= 1);

figure();
imshow(image, 'InitialMagnification', 500);

%% compute cone response example
[~, ~, linear, coneVec] = retina.compute(image);
retina.visualizeExcitation();
retina.visualizeOI();

%% render matrix
render = retina.forwardRender(imageSize);
render = double(render);

%% image reconstruction: cone response
% compute (noise-free) avearge cone response
response = render * linear(:);

% cone response with added Poisson noise
response = poissrnd(response);

%% reconstruction: running the estimator
% construct image estimator
prior = load('sparsePrior.mat');
regPara = 0.05; stride = 4; 

estimator = PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
reconImage = estimator.runEstimate(response, 'maxIter', 1000, 'display', 'iter');

%% show results
% compare the original and reconstructed images
figure();
subplot(1, 2, 1);
imshow(image, 'InitialMagnification', 500);
title('Original')

% reconstructed image is in linear pixel space
% convert it back to RGB space
recon = gammaCorrection(reconImage, display);

subplot(1, 2, 2);
imshow(recon, 'InitialMagnification', 500);
title('Reconstruction')
set(gcf,'position',[0, 0, 800, 400]);
