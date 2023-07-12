

%% Use ISETTreeshrew and create forward model for reconstruction
%
% Generate a generic cone response object
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.1);

% Get trew shrew optics and mosaic
fovDegs = 8;
shrewPSF = oiTreeShrewCreate('pupilDiameterMM', 2.0);
shrewMosaic = coneMosaicTreeShrewCreate(...
    shrewPSF.optics.micronsPerDegree, ...
    'fovDegs', fovDegs, ...
    'sConeMinDistanceFactor', 2, ...
    'integrationTime', 0.1);
shrewMosaic.noiseFlag = 'none';

% Write these into the cone response object
retina.Mosaic = shrewMosaic;
retina.PSF = shrewPSF;
retina.FovealDegree = fovDegs;

%% Read in an imaget to test with
%
% Use a standard Matlab image and reduce its size.
imageSize = [32, 32, 3];
[X,cmap] = imread('corn.tif');
theRawImage = ind2rgb(X,cmap);
image   = imresize(im2double(theRawImage), 0.5);
patch = sampleImage(image, imageSize(1));
patch(patch < 0) = 0;
patch(patch > 1) = 1;

% Show the test image
imshow(patch, 'InitialMagnification', 500);

%% Calculate cone responses to the test image
[~, ~, linearImage, coneVec] = retina.compute(patch);
retina.visualizeOI();
retina.visualizeExcitation();

%% Compute the render matrix for the tree shrew eye
%
% The render matrix is the linear matrix that maps between
% the image in vector form and the vector of mean cone excitations.
%
% The way we do this is to generate a set of "training" pairs
% between image vectors and cone excitation vectors, and then
% use linear regression to find the matrix that maps between them.
%
% The training set needs to be at least as larger as the number
% of pixels in the image times the number of color channels.
minNTrain = prod(imageSize);
nTrain = 1.3*minNTrain;

% Allocate space for training set
allConeVec   = zeros(nTrain, length(coneVec));
allLinearImg = zeros(nTrain, length(linearImage(:)));

% Calculate the training set 
%
% In practice, we found that white noise images produce better results
% and are easier to work with (since we can simply generate them on the
% fly).
parfor idx = 1:nTrain
   
    input = rand(imageSize);
    [~, ~, linearImage, coneVec] = retina.compute(input);
    
    allConeVec(idx, :)   = coneVec;
    allLinearImg(idx, :) = linearImage(:);
end

%% Linear regression
regEstimator = RegressionEstimator(allLinearImg, allConeVec);
renderMatrix = regEstimator.W';

%% Test the render matrix
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
