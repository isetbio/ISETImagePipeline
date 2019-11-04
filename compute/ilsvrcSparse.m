%% Load display setup
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Create the cone mosaic, test run
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.5, 'display', display.CRT12BitDisplay);

imageSize = [180, 180, 3];
testImage = rand(imageSize);
[~, ~, testLinearImage, testConeVec] = retina.compute(testImage);

inteTime = retina.Mosaic.integrationTime;

%% Visualization
showPlot = true;
if showPlot
    retina.visualizeMosaic();
    
    figure();
    imshow(testImage, 'InitialMagnification', 500);
    retina.visualizeOI();
    retina.visualizeExcitation();
end

%% Compute render matrix with basis function
retina.Mosaic.integrationTime = inteTime * 50;
render = zeros(length(testConeVec), length(testLinearImage(:)));

parfor idx = 1:length(testLinearImage(:))
    input = zeros(size(testLinearImage));
    input(idx) = 1.0;
    
    [~, ~, linear, coneVec] = retina.compute(input);
    render(:, idx) = coneVec;
end

retina.Mosaic.integrationTime = inteTime;
render = render ./ 50;

%% Zero thresholding / Scaling for numerical precision
testRender = render;
testRender(testRender < 1) = 0;

%% Load dataset
projectName  = 'ISETImagePipeline';
thisImageSet = 'ILSVRC_mat';
dataBaseDir  = getpref(projectName, 'dataDir');

imageSet = load(fullfile(dataBaseDir, thisImageSet, 'imageDataLinear.mat'));
imageSet = imageSet.imageDataLinear;

%% Whitening, SVD
[Z, U, SIG, MU] = whitening(imageSet, 'svd');

%% RICA analysis
nBasis = 16 * 16 * 3;
Mdl = rica(Z, nBasis, 'IterationLimit', 1e4, 'VerbosityLevel', 1, 'GradientTolerance', 1e-8, 'StepTolerance', 1e-8);
regBasis = U * diag(sqrt(SIG)) * Mdl.TransformWeights;

%% Visualization
[~] = visualizeBasis(regBasis, 16, size(regBasis, 2), false);

%% Option 1: Sample test image
projectName  = 'ISETImagePipeline';
thisImageSet = 'ILSVRC';
dataBaseDir  = getpref(projectName, 'dataDir');
imageName    = 'ILSVRC2017_test_00000021.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.5);

patch = sampleImage(image, imageSize(1));
patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 500);

%% Option 2: Generate Uniform Field
patch = zeros(imageSize);

patch(:, :, 1) = 1.0; % R channel
patch(:, :, 2) = 1.0; % G channel
patch(:, :, 3) = 1.0; % B channel

imshow(patch, 'InitialMagnification', 500);

%% Generate Cone Vec
[~, ~, patchLinear, patchConeVec] = retina.compute(patch);
showPlot = true;
if showPlot
    retina.visualizeMosaic();
    
    figure();
    imshow(patch, 'InitialMagnification', 500);
    retina.visualizeOI();
    retina.visualizeExcitation();
end

figure();
scatter(patchConeVec, renderMtx * patchLinear(:));
grid on; hold on;
refPoint = [-500, 4000];
plot(refPoint, refPoint);
axis square;
xlim(refPoint);
ylim(refPoint);

%% Reconstruction
estimator  = SparsePatchEstimator(renderMtx, inv(regBasis), MU', 0.05, 2, size(patch));
reconImage = estimator.estimate(renderMtx * patchLinear(:), 2.5e3, ones([numel(patch), 1]) * 0.5);

%% Show plot
figure();
subplot(1, 2, 1);
imshow(patch, 'InitialMagnification', 500); title('input');
subplot(1, 2, 2);
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 500); title('reconstruction');

%% Option 3: Single Cone Stimulation - Set Up
% regularization parameter [0.1, 1, 2, 10]
imageSize = [140, 140, 3];
estimator = SparsePatchEstimator(renderMtx, inv(regBasis), MU', 2, 1, imageSize);

%% Calculation
nIter    = 5;
coneType = 'L';
for idx = 1:nIter
    stimulus = 0.25 * ones(imageSize);
    [~, ~, linStim, coneRes]  = retina.compute(stimulus);
    coneResponse = renderMtx * linStim(:);
    % [coneResponse, coneCount] = retina.coneExcitationRnd(8, coneType);
    
    stimIdx = randi(length(coneResponse));
    coneResponse(stimIdx) = coneResponse(stimIdx) * 8;
    
    fig = figure(idx);
    % subplot(1, 2, 1);
    % retina.visualizeExcitation(true);
    
    reconLasso = estimator.estimate(coneResponse, 2.5e3, ones([prod(imageSize), 1]) * 0.25);
    reconLasso = invGammaCorrection(reshape(reconLasso, imageSize), display.CRT12BitDisplay);
    
    % subplot(1, 2, 2);
    imshow(reconLasso, 'InitialMagnification', 500);
    % title('Sparse Prior Reconstruction');
    % suptitle(sprintf('L: %d, M: %d, S: %d', coneCount(1), coneCount(2), coneCount(3)));
end

%% Two line stimulus
[~, allCone] = retina.computeWithOI(theOINoLca);

retina.visualizeOI();
retina.visualizeExcitation();

%% Reconstruction
% TODO? thresholding before reconstruction
imageSize = [140, 140, 3];
estimator  = SparsePatchEstimator(renderMtx, inv(regBasis), MU', 1, 1, imageSize);
reconImage = estimator.estimate(allCone, 2.5e3, ones([prod(imageSize), 1]) * 0.0);
reconImage = invGammaCorrection(reshape(reconImage, imageSize), display.CRT12BitDisplay);

%% Show Image
imshow(reconImage, 'InitialMagnification', 500);

%% Check cone excitation
retina.compute(reconImage);
retina.visualizeOI();
retina.visualizeExcitation();

%% Reconstruction for regular scence (with optics)
% Scene -> OI -> Cone Response -> Recon -> OI -> Cone Response
[~, allCone] = retina.computeWithScene(scene);

%% Visualization
retina.visualizeOI();
retina.visualizeExcitation();

%% Reconstruction with this stimulus
imageSize = [140, 140, 3];
estimator  = SparsePatchEstimator(renderMtx, inv(regBasis), MU', 1, 1, imageSize);
reconImageLCA = estimator.estimate(allCone, 2.5e3, ones([prod(imageSize), 1]) * 0.0);
reconImageLCA = invGammaCorrection(reshape(reconImageLCA, imageSize), display.CRT12BitDisplay);

%% Visualization
imshow(reconImageLCA, 'InitialMagnification', 500);

%% Check cone excitation
retina.compute(reconImageLCA);
retina.visualizeOI();
retina.visualizeExcitation();
