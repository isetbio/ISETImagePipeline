%% Load display setup
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Create the cone mosaic, test run
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.5, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

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
retina.Mosaic.integrationTime = inteTime * 10;
render = zeros(length(testConeVec), length(testLinearImage(:)));

parfor idx = 1:length(testLinearImage(:))
    input = zeros(size(testLinearImage));
    input(idx) = 1.0;
    
    [~, ~, linear, coneVec] = retina.compute(input);
    render(:, idx) = coneVec;
end

retina.Mosaic.integrationTime = inteTime;
render = render ./ 10;

%% Visualization
load('./sparsePrior.mat');
[~] = visualizeBasis(regBasis, 16, size(regBasis, 2), false);

%% Option 1: Sample test image
projectName  = 'ISETImagePipeline';
thisImageSet = 'ILSVRC';
dataBaseDir  = getpref(projectName, 'dataDir');
imageName    = 'ILSVRC2017_test_00000025.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.65);

patch = sampleImage(image, imageSize(1));
patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 400);

%% Generate Cone Vec
[~, ~, patchLinear, patchConeVec] = retina.compute(patch);
showPlot = true;
if showPlot        
    figure();
    imshow(patch, 'InitialMagnification', 400);
    retina.visualizeOI();
    retina.visualizeExcitation();
end

figure();
scatter(patchConeVec, render * patchLinear(:));
grid on; hold on;
refPoint = [-500, 4000];
plot(refPoint, refPoint);
axis square;
xlim(refPoint);
ylim(refPoint);

%% Reconstruction
estimator  = SparsePatchEstimator(render, inv(regBasis), MU', 0.1, 4, imageSize);

meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

reconImage = estimator.estimate(response, 4e3, rand([prod(imageSize), 1]));

%% Compare results
% Show plot
figure();

subplot(1, 2, 1);
imshow(patch, 'InitialMagnification', 400); 
title('input');

subplot(1, 2, 2);
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 400); 
title('reconstruction');
