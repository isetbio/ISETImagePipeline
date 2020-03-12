% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

load('retinaRenderFov.mat');
render = double(renderFov);

regPara = 1e-3;
estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
outputImage = zeros(size(inputLinear));

parfor idx = 1:size(inputLinear, 1)
    
    input = inputLinear(idx, :, :, :);
    coneVec = render * input(:);
    
    outputImage(idx, :, :, :) = estimator.estimate(coneVec, 1.5e3, rand([prod(imageSize), 1]), true);
    
end

save('reconOutput.mat', 'outputImage', '-v7.3');

%% Analysis across retina: load data
load('sparsePrior.mat');

nMosaic = 5;
mosaicArray = cell(1, nMosaic);
renderArray = cell(1, nMosaic);

load('retinaRenderFov.mat');
mosaicArray(1) = {retinaFov};
renderArray(1) = {double(renderFov)};

load('retinaRender5.mat');
mosaicArray(2) = {retina5};
renderArray(2) = {double(render5)};

load('retinaRender10.mat');
mosaicArray(3) = {retina10};
renderArray(3) = {double(render10)};

load('retinaRender20.mat');
mosaicArray(4) = {retina20};
renderArray(4) = {double(render20)};

load('retinaRender40.mat');
mosaicArray(5) = {retina40};
renderArray(5) = {double(render40)};

%% Plot cone mosaic
% for idx = 1:nMosaic
%     retina = mosaicArray{idx};
%     retina.visualizeMosaic();
% end

%% Reconstruction of uniform patches: Stimulus
color1 = ([252, 123, 20] ./ 255);
color2 = ([20, 144, 252] ./ 255);
colors = [color1; color2];

% imshow(patchMaker(imageSize, color1));
% imshow(patchMaker(imageSize, color2));

%% Reconstruction of uniform patches: Reconstruction
regPara = 1e-3;
output = zeros([size(colors, 1), nMosaic, imageSize]);
for i = 1 : size(colors, 1)
    patch = patchMaker(imageSize, colors(i, :));
    [~, ~, linear] = retinaFov.compute(patch);
    for j = 1 : nMosaic
        render = renderArray{j};
        estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
        
        coneVec = render * patch(:);
        output(i, j, :, :, :) = estimator.estimate(coneVec, 20, rand([prod(imageSize), 1]), true);
    end
end

%% Visualization
figure();
plotAxis = tight_subplot(2, 5, [0.01, 0.01], 0.01, 0.01);
plotID   = 1;

dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));
imageSize = [128, 128, 3];

for i = 1 : size(colors, 1)
    for j = 1 : nMosaic
        axes(plotAxis(plotID));
        image = reshape(output(i, j, :, :), imageSize);
        imshow(invGammaCorrection(image, display.CRT12BitDisplay), 'InitialMagnification', 400);
        
        plotID = plotID + 1;
    end
end

%% helper function
function patch = patchMaker(imageSize, color)
nDim = imageSize(1);
patch = zeros(imageSize);
for idx = 1:3
    patch(:, :, idx) = color(idx) .* ones(nDim, nDim);
end
end
