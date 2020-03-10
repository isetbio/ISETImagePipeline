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