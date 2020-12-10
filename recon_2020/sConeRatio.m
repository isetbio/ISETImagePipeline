%% define constant
imageSize = [64, 64, 3];
display = load('display.mat');
prior   = load('sparsePrior.mat');

%% generate a cone mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display.CRT12BitDisplay, 'pupilSize', 2.5);

%% change the S cone proportion
% and generate the corresponding render matrix
retina.resetSCone();
ratio = [0, 0.01, 0.05, 0.1, 0.25, 0.50, 0.75, 0.9];

mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));
for idx = 1:length(ratio)
    % manipulate the number of S cone in the mosaic
    retina.resetSCone();
    retina.reassignSCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize, false, true, false);
    renderArray(idx) = {double(renderMtx)};
end

%% load images
nImage = 8;
input = zeros([nImage, imageSize]);

fileType = '.jpeg';
for idx = 1:nImage
    fileName = strcat(num2str(idx), fileType);
    filePath = fullfile('.', 'images', fileName);
    image = imresize(im2double(imread(filePath)), 0.25);
    
    image = sampleImage(image, imageSize(1));
    image = image - min(image(:));
    image = image ./ max(image(:));
    
    [~, ~, linearImage] = retina.compute(image);
    input(idx, :, :, :) = linearImage;
end

%% reconstruction
output = zeros([length(ratio), nImage, imageSize]);
parfor i = 1:length(ratio)
    render = renderArray{i};
    estimator = ...
        PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', 1e-4, 4, imageSize);
    
    for j = 1:nImage
        image = input(j, :, :, :);
        resp  = render * image(:);
        
        recon = estimator.estimate(resp, 500, rand([prod(imageSize), 1]), true, 1.0, 'final');
        output(i, j, :, :, :) = recon;
    end
end

%% RMSE and visualization

%% analysis without chromatic abberation

