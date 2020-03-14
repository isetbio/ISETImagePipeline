%% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');
load('retinaRender10.mat');
render = double(render10);

dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

regPara = 1e-2;
estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
outputImage = zeros(size(inputLinear));

%% Reconstruction
parfor idx = 1:size(inputLinear, 1)
    image = reshape(inputLinear(idx, :, :, :), imageSize);
    
    fltSD = 5;
    image = imgaussfilt(image, fltSD);
    
    coneVec = render * image(:);
    outputImage(idx, :, :, :) = estimator.estimate(coneVec, 10, rand([prod(imageSize), 1]), true);
end

fprintf('Finish reconstruction, save results... Done! \n');
save('blurRecon.mat', 'outputImage', '-v7.3');

%% Show results
for idx = 1:size(inputLinear, 1)
    image = reshape(inputLinear(idx, :, :, :), imageSize);
    
    fltSD = 5;
    image = imgaussfilt(image, fltSD);
    image = invGammaCorrection(image, display.CRT12BitDisplay);
    
    figure();
    subplot(1, 2, 1);
    imshow(image, 'InitialMagnification', 400);
    
    subplot(1, 2, 2);
    output = reshape(outputImage(idx, :, :, :), imageSize);
    imshow(invGammaCorrection(output, display.CRT12BitDisplay), 'InitialMagnification', 400);
end