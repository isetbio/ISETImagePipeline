%% Create Retina object
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.1);

% load from file - treeshrew mosaic and PSF object
threwMosaic.noiseFlag = 'none';
retina.Mosaic = threwMosaic;
retina.PSF = threwPSF;
retina.FovealDegree = 5;

%% Compute mosaic response
imageSize = [128, 128, 3];
image   = imresize(im2double(imread('./input.jpeg')), 0.5);
patch = sampleImage(image, imageSize(1));

%% Show input image
patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 500);

%% Calculate response
retina.compute(patch);
retina.visualizeExcitation();
retina.visualizeOI();

%% Image reconstruction
% using small image size so the demo will run fast
imageSize = [32, 32, 3];

testImage = rand(imageSize);
[~, ~, linearImage, coneVec] = retina.compute(testImage);

%% Generate a "training set" for regression
nTrain = 4e3;
allConeVec   = zeros(nTrain, length(coneVec));
allLinearImg = zeros(nTrain, length(linearImage(:)));

parfor idx = 1:nTrain
    % In practice, we found that white noise images produce better results
    % and are easier to work with (since we can simply generate them on the
    % fly).
    input = rand(imageSize);
    [~, ~, linearImage, coneVec] = retina.compute(input);
    
    allConeVec(idx, :)   = coneVec;
    allLinearImg(idx, :) = linearImage(:);
end

%% Linear regression
regEstimator = RegressionEstimator(allLinearImg, allConeVec);
renderMatrix = regEstimator.W';

%% Test the render matrix
image   = imresize(im2double(imread('./input.jpeg')), 0.25);
patch = sampleImage(image, imageSize(1));

patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 500);

%% Calculation
[~, ~, linearImage, coneVec] = retina.compute(patch);
coneVecRender = renderMatrix * linearImage(:);

figure();
scatter(coneVec, coneVecRender);

%% Image reconstruction
estimator  = SparsePatchEstimator(renderMatrix, inv(regBasis), MU', 0.1, 2, imageSize);
reconImage = estimator.estimate(coneVec, 2.5e3, rand([prod(imageSize), 1]));

%% Compare reconstruction
figure();
subplot(1, 2, 1);
imshow(patch, 'InitialMagnification', 500);

subplot(1, 2, 2);
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 500);
