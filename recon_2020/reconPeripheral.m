%% Load the display setup and create mosaic object
display = displayCreate('CRT12BitDisplay');

% Generate cone mosaic - [eccX, eccY] deg ecc
% with the new cone mosaic (cMosaic class)
eccX = 18.0; eccY = 18.0;
retina = ConeResponseCmosaic...
    (eccX, eccY, 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 6);

%% Show mosaic and optical PSF
retina.visualizeMosaic();
retina.visualizePSF();

%% Load an example image
imageSize = [128, 128, 3];
image = im2double(imread('images/2.jpeg'));
image = imresize(image, imageSize(1) / size(image, 1));

% Normalize to insure pixels are in the range of [0, 1]
image = (image - min(image(:)));
image = image ./ max(image(:));

assert(min(image(:)) >= 0);
assert(max(image(:)) <= 1);

figure();
imshow(image, 'InitialMagnification', 500);

%% Test on example image
[allCone, linear] = retina.compute(image);

retina.visualizeOI();
retina.visualizeExcitation();

%% Render matrix
render = retina.forwardRender(imageSize, true, false);
render = double(render);

%% Image reconstruction: Compute cone response
% compute (noise-free) avearge cone response
response = render * linear(:);

% cone response with added Poisson noise
% response = poissrnd(response);

%% Reconstruction: Running the estimator
% Construct image estimator
prior = load('sparsePrior.mat');
regPara = 2e-3; stride = 4; 

estimator = PoissonSparseEstimator...
    (render, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
reconImage = estimator.runEstimate(response, 'maxIter', 500, 'display', 'iter');

%% Show results
% Compare the original and reconstructed images
figure();
subplot(1, 2, 1);
imshow(image, 'InitialMagnification', 500);
title('Original')

% Reconstructed image is in linear pixel space
% Convert it back to RGB space
recon = gammaCorrection(reconImage, display);

subplot(1, 2, 2);
imshow(recon, 'InitialMagnification', 500);
title('Reconstruction')
set(gcf,'position',[0, 0, 800, 400]);
