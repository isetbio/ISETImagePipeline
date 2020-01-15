%% Load display setup & constant
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

imageSize = [128, 128, 3];

%% Generate cone mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'pupilSize', 2.5);

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

%% Manipulation of M cone ratio
retina.resetCone();

ratio = [0.32, 0.2, 0.15, 0.10, 0.05, 0.01, 0.001, 0];
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('M Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignMCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize);
    renderArray(idx) = {renderMtx};
end

%% Manipulation of L cone ratio
retina.resetCone();

ratio = [0.6, 0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0];
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('L Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignLCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize);
    renderArray(idx) = {renderMtx};
end

%% Manipulate of S cone ratio
retina.resetCone();

ratio = [0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95];
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));
for idx = 1:length(ratio)
    retina.resetSCone();
    fprintf('L Cone Ratio: %.4f \n', ratio(idx));
    
    retina.reassignSCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize);
    renderArray(idx) = {renderMtx};
end

%% Input imageset for reconstruction
load('./inputImage_128.mat');
load('./sparsePrior.mat');

nImage = 12;
for idx = 1:nImage
    image = reshape(inputLinear(idx, :, :, :), imageSize);
    
    figure();
    imshow(invGammaCorrection(image, display.CRT12BitDisplay), 'InitialMagnification', 500);
end

%% Generate retina mosaic and cone response
retinaIdx = 1;

retina.Mosaic.pattern = mosaicArray{retinaIdx};
render = double(renderArray{retinaIdx});

imageIdx = 11;
image = invGammaCorrection(reshape(inputLinear(imageIdx, :, :, :), imageSize), display.CRT12BitDisplay);

[~, ~, imageLinear, imageCone] = retina.compute(image);
retina.visualizeExcitation();
retina.visualizeMosaic();

figure(); hold on;
scatter(imageCone, render * imageLinear(:));
plot(xlim, ylim, '--k', 'LineWidth', 2);
xlim([0, 1500]); ylim([0, 1500]);
axis square;

%% Reconstruction test: Image Dataset
inputTest = reshape(inputLinear(imageIdx, :, :, :), imageSize);

estimator = PoissonSparseEstimator(render, inv(regBasis), MU', 5e-4, 4, imageSize);
reconTest = estimator.estimate(render * inputTest(:), 5e2, rand([prod(imageSize), 1]), true);

figure(); subplot(1, 2, 1);
imshow(invGammaCorrection(inputTest, display.CRT12BitDisplay), 'InitialMagnification', 500);

subplot(1, 2, 2);
imshow(invGammaCorrection(reconTest, display.CRT12BitDisplay), 'InitialMagnification', 500);

%% Reconstruction with different retina
load('./sparsePrior.mat');
load('./inputImage.mat');
load('./retina_render.mat');

outputArray = cell(1, length(mosaicArray));
regPara = 5e-4;
parfor idx = 1:length(mosaicArray)
    fprintf('Reconstruction for Retina %d \n', idx);
    outputArray(idx) = {imageRecon(double(renderArray{idx}), inputLinear, regBasis, MU, regPara, imageSize)};
end

%% Show reconstruction
retinaIdx = 1;
output = outputArray{retinaIdx};

figure();
for idx = 1:size(output, 1)
    subplot(2, 5, idx);
    recon = reshape(output(idx, :, :, :), imageSize);
    imshow(invGammaCorrection(recon, display.CRT12BitDisplay), 'InitialMagnification', 500);
end

%% Helper function
function outputImage = imageRecon(render, inputImage, regBasis, MU, regPara, imageSize)

estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
outputImage = zeros(size(inputImage));
for idx = 1:size(inputImage, 1)
    
    input = inputImage(idx, :, :, :);
    coneVec = render * input(:);
    outputImage(idx, :, :, :) = estimator.estimate(coneVec, 1e3, rand([prod(imageSize), 1]), true);
    
end
end
