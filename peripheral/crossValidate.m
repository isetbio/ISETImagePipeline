%% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

load('retinaRender10.mat');
render = double(render10);

regParas = [1e-5, 1e-4, 5e-4, 1e-3, 1e-2, 5e-2, 0.1, 1];
outputImage = zeros([length(regParas), size(inputLinear)]);

%% Run reconstruction
for regIdx = 1:length(regParas)
    regPara = regParas(regIdx);
    estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
    parfor idx = 1:size(inputLinear, 1)
        
        input = inputLinear(idx, :, :, :);
        coneVec = render * input(:);
        
        outputImage(regIdx, idx, :, :, :) = estimator.estimate(coneVec, 2e3, rand([prod(imageSize), 1]), true);
        
    end
end

save('cvOutput.mat', 'outputImage', '-v7.3');

%% Show results
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

figure();
plotAxis = tight_subplot(length(regParas), 12, [0.01, 0.01], 0.01, 0.01);

plotID = 1;
for i = 1:length(regParas)
    for j = 1:size(inputLinear, 1)
        axes(plotAxis(plotID));
        image = reshape(outputImage(i, j, :, :, :), imageSize);
        imshow(invGammaCorrection(image, display.CRT12BitDisplay), 'InitialMagnification', 400);
        
        plotID = plotID + 1;
    end
end