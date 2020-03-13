%% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');
load('retinaRender10.mat');
render = double(render10);

dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

regPara = 1e-2;
outputImage = zeros(size(inputLinear));

%% Reconstruction
for idx = 1:size(inputLinear, 1)
    image = reshape(inputLinear(idx, :, :, :), imageSize);
    
    fltSD = 5;    
    image = imgaussfilt(image, fltSD);
    
    % image = invGammaCorrection(image, display.CRT12BitDisplay);
    % retina10.compute(image);
    % retina10.visualizeOI();
        
    coneVec = render * image(:);
    outputImage(idx, :, :, :) = estimator.estimate(coneVec, 2e3, rand([prod(imageSize), 1]), true);
end