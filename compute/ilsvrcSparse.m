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
showPlot = true;
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
[Z, U, SIG, MU] = whitening(imageTr, 'svd');

%% RICA analysis
nBasis = 16 * 16;
Mdl = rica(Z, nBasis, 'IterationLimit', 5e3, 'VerbosityLevel', 1, 'GradientTolerance', 1e-4, 'StepTolerance', 1e-4);
regBasis = U * diag(sqrt(SIG)) * Mdl.TransformWeights;

%% Sample test image
projectName  = 'ISETImagePipeline';
thisImageSet = 'ILSVRC';
dataBaseDir  = getpref(projectName, 'dataDir');
imageName    = 'ILSVRC2017_test_00000021.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.4);

patch = sampleImage(image, 140);
patch(patch < 0) = 0;
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
estimator  = SparsePatchEstimator(renderMtx, inv(regBasis), mu', 0.1, 2, size(patch));
reconImage = estimator.estimate(renderMtx * patchLinear(:), 5e2, patchLinear(:));

%% Show plot
figure();
subplot(1, 2, 1); title('input');
imshow(patch, 'InitialMagnification', 500);
subplot(1, 2, 2); title('reconstruction');
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 500);