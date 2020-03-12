% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

load('retinaRender10.mat');
render = double(render10);

regParas = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1];
outputImage = zeros([length(regPara), size(inputLinear)]);

for regIdx = 1:length(regParas)
    regPara = regParas(regIdx);
    estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
    parfor idx = 1:size(inputLinear, 1)
        
        input = inputLinear(idx, :, :, :);
        coneVec = render * input(:);
        
        outputImage(regIdx, idx, :, :, :) = estimator.estimate(coneVec, 1.0e3, rand([prod(imageSize), 1]), true);
        
    end
end

save('cvOutput.mat', 'outputImage', '-v7.3');