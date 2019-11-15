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
imageName    = 'ILSVRC2017_test_00000173.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.3);

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
refPoint = [0, 15000];
plot(refPoint, refPoint);
axis square;
xlim(refPoint);
ylim(refPoint);

%% Reconstruction with Gaussian likelihood and sparse prior
estimator = SparsePatchEstimator(render, inv(regBasis), MU', 2e1, 4, imageSize);

meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

reconImage = estimator.estimate(response, 2.5e3, rand([prod(imageSize), 1]));

%% Results of reconstruction
% Show plot
figure();

subplot(1, 2, 1);
imshow(patch, 'InitialMagnification', 400); 
title('input');

subplot(1, 2, 2);
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 400); 
title('reconstruction');

%% Reconstruction with Poisson likelihood and sparse prior
estimator = PoissonSparseEstimator(render, inv(regBasis), MU', 1, 4, imageSize);

meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

reconImage = estimator.estimate(response, 2.5e3, rand([prod(imageSize), 1]));

%% Results of reconstruction
% Show plot
figure();

subplot(1, 2, 1);
imshow(patch, 'InitialMagnification', 400); 
title('input');

subplot(1, 2, 2);
imshow(invGammaCorrection(reconImage, display.CRT12BitDisplay), 'InitialMagnification', 400); 
title('reconstruction');

%% Effect of regularization - Gaussian distribution
meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

regPara = [10, 20, 40, 60, 80];
reconImageGauss = zeros([5, imageSize]);

for idx = 1:5
    estimator = SparsePatchEstimator(render, inv(regBasis), MU', regPara(idx), 4, imageSize);
    reconImageGauss(idx, :, :, :) = estimator.estimate(response, 2.5e3, rand([prod(imageSize), 1]));
end

%% Effect of regularization - Poisson distribution
meanResp = render * patchLinear(:);
response = poissrnd(meanResp);

regPara = [1e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1, 0.5, 1];
reconImagePoiss = zeros([3, imageSize]);

for idx = 1 : 3
    estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara(idx), 4, imageSize);
    reconImagePoiss(idx, :, :, :) = estimator.estimate(response, 3e3, rand([prod(imageSize), 1]), true);
end

%% Results of reconstruction
% Show plot
figure();
imshow(patch, 'InitialMagnification', 400); 
title('input');

mseRecon = zeros(1, 8);
for idx = 1 : 8
    figure();
    
    reconImage = reshape(reconImagePoiss(idx, :, :, :), imageSize);
    reconImage = invGammaCorrection(reconImage, display.CRT12BitDisplay);
    imshow(reconImage, 'InitialMagnification', 400); pause(0.5);
    title(strcat('reconstruction - lambda:', num2str(regPara(idx))));
    
    mseRecon(idx) = norm(patch(:) - reconImage(:));
end

figure();
plot(log(regPara), mseRecon, 'k', 'LineWidth', 2);
xlabel('Prior Weight');
ylabel('Mean Squared Error');
xticks(log(regPara)); 
xticklabels(arrayfun(@(x) num2str(x), regPara, 'UniformOutput', false));
