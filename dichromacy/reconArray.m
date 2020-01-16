function outputArray = reconArray(inputImage, renderArray, regBasis, MU, regPara, imageSize)

outputArray = cell(1, length(renderArray));

parfor idx = 1:length(renderArray)
    fprintf('Reconstruction for Retina %d \n', idx);
    outputArray(idx) = {imageRecon(double(renderArray{idx}), inputImage, regBasis, MU, regPara, imageSize)};
end

end

% Helper function
function outputImage = imageRecon(render, inputImage, regBasis, MU, regPara, imageSize)

estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
outputImage = zeros(size(inputImage));

for idx = 1:size(inputImage, 1)
    
    input = inputImage(idx, :, :, :);
    coneVec = render * input(:);
    outputImage(idx, :, :, :) = estimator.estimate(coneVec, 1e3, rand([prod(imageSize), 1]), true);
    
end
end