%% Load display setup
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Create the cone mosaic, test run
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display.CRT12BitDisplay);

testImage = rand([140, 140, 3]);
[~, ~, testLinearImage, testConeVec] = retina.compute(testImage);

%% Visualization
showPlot = false;
if showPlot
    retina.visualizeMosaic();
    
    figure();
    imshow(testImage, 'InitialMagnification', 500);
    retina.visualizeOI();
    retina.visualizeExcitation();
end

%% Compute render matrix with basis function
render = zeros(length(testConeVec), length(testLinearImage(:)));

parfor idx = 1:length(testLinearImage(:))
    input = zeros(size(testLinearImage));
    input(idx) = 1.0;
    
    [~, ~, linear, coneVec] = retina.compute(input);
    render(:, idx) = coneVec;
end

%% Zero thresholding for numerical precision
renderMtx = render;
renderMtx(renderMtx < 0.75) = 0;

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
image   = imresize(im2double(imread(fileDir)), 0.4);

patch = sampleImage(image, 140);
patch(patch < 0) = 0;
imshow(patch, 'InitialMagnification', 500);

%% Option 2: Generate Uniform Field
patch = zeros(140, 140, 3);

patch(:, :, 1) = 1; % R channel
patch(:, :, 2) = 0.6; % G channel
patch(:, :, 3) = 0.25; % B channel

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
