function outputArray = reconArray(inputImage, renderArray, regBasis, MU, regPara, imageSize)

outputArray = cell(1, length(renderArray));

% parfor for each render matrix
% parfor idx = 1:length(renderArray)
for idx = 1:length(renderArray)
    fprintf('Reconstruction for Retina %d \n', idx);
    outputArray(idx) = {imageRecon(double(renderArray{idx}), inputImage, regBasis, MU, regPara, imageSize)};
end

end

% Helper function
function outputImage = imageRecon(render, inputImage, regBasis, MU, regPara, imageSize)

estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
outputImage = zeros(size(inputImage));

% alternatively, parfor for each image
% for idx = 1:size(inputImage, 1)
parfor idx = 1:size(inputImage, 1)
    input = inputImage(idx, :, :, :);
    coneVec = render * input(:);
    outputImage(idx, :, :, :) = estimator.estimate(poissrnd(coneVec), 1e3, ...
                                rand([prod(imageSize), 1]), true, 1.0, 'final', false);
end
end