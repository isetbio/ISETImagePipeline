%% Effect of S Cone Proportion on Image Reconstruction
% Analysis 1: With regular optics

%% define constant
imageSize = [64, 64, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

%% generate a cone mosaic
% analysis with normal optics
pupilSize = 2.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

%% load images
nImage = 10;
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

%% change the S cone proportion
% and generate the corresponding render matrix
ratio = [0, 0.01, 0.05, 0.1, 0.25, 0.40, 0.50, 0.60, 0.75, 0.9];
[~, renderArray] = computeRender(ratio, retina, imageSize);

%% reconstruction
regPara = 5e-3;
output = computeRecon(input, renderArray, prior, regPara, imageSize);

%% show results
plotResults(input, output, ratio, display, imageSize);

%% Analysis 2: Turn off LCA or use diffraction-limited optics
% compute optics with 'no lca' / 'diffraction-limit' flag
pupilSize = 2.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

testImage = gammaCorrection(reshape(input(1, :, :, :), imageSize), display);

% compute OI and excitation
retina.compute(testImage);
retina.visualizeOI();
retina.visualizeExcitation();

%% change optics, turn off LCA / use diffraction limited optics
option = 'diffraction-limit';
switch option
    case 'no-lca'
        retina.PSF = ConeResponse.psfNoLCA(pupilSize);
    case 'diffraction-limit'
        retina.PSF = ConeResponse.psfDiffLmt(pupilSize);
    otherwise
        error('optics option not valid');
end

% compute OI and excitation
retina.compute(testImage);
retina.visualizeOI();
retina.visualizeExcitation();

%% change the S cone proportion
% and generate the corresponding render matrix
ratio = [0, 0.01, 0.05, 0.1, 0.25, 0.40, 0.50, 0.60, 0.75, 0.9];
[~, renderArray] = computeRender(ratio, retina, imageSize);

%% reconstruction
regPara = 5e-3;
output = computeRecon(input, renderArray, prior, regPara, imageSize);

%% show results
plotResults(input, output, ratio, display, imageSize);

%% Analysis 3: Turn off LCA, Lens and Macular pigment
% define constant & cone mosaic & load images
imageSize = [64, 64, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

pupilSize = 2.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

nImage = 10;
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

%% Turn off optics
retina.PSF = ConeResponse.psfDiffLmt(pupilSize);

%% change the optical density of the lens
optics = retina.PSF;

lens0 = oiGet(optics, 'lens');
wls = lens0.wave;

lensUnitDensity1 = lens0.unitDensity;
lensUnitDensity1 = zeros(size(lensUnitDensity1));
lensPeakDensity1 = 1;

lens1 = Lens('wave', wls, ...
    'unitDensity', lensUnitDensity1, 'density', lensPeakDensity1);
retina.PSF = oiSet(optics, 'lens', lens1);

%% change the macular density
macular0 = retina.Mosaic.macular;
wls = macular0.wave;

macularUnitDensity1 = zeros(size(macular0.unitDensity));
macularDensity1 = macular0.density;

macular1 = Macular('wave', wls, 'unitDensity', macularUnitDensity1, 'density', macularDensity1);
retina.Mosaic.macular = macular1;

% change pigment
pigment = retina.Mosaic.pigment;
pigment.opticalDensity = [0.2, 0.2, 0.5];

%% change the S cone proportion
% and generate the corresponding render matrix
ratio = [0, 0.01, 0.05, 0.1, 0.25, 0.40, 0.50, 0.60, 0.75, 0.9];
[~, renderArray] = computeRender(ratio, retina, imageSize);

%% reconstruction
regPara = 5e-3;
output = computeRecon(input, renderArray, prior, regPara, imageSize);

%% show results
plotResults(input, output, ratio, display, imageSize);

%% Analysis 4: With two class of cone type
% define constant
imageSize = [64, 64, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

% generate a cone mosaic
% analysis with normal optics
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

%% Turn off optics
% retina.PSF = ConeResponse.psfDiffLmt(pupilSize);

%% load images
nImage = 10;
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

%% change the S cone proportion
% and generate the corresponding render matrix
ratio = [0.0, 0.01, 0.05, 0.1, 0.25, 0.40, 0.50, 0.60, 0.75, 0.9, 0.95];
[~, renderArray] = computeRenderDichroma(ratio, retina, imageSize, 'noMCone');

% reconstruction
regPara = 5e-3;
output = computeRecon(input, renderArray, prior, regPara, imageSize);

%% show results
plotResults(input, output, ratio, display, imageSize);

%% helper function for plotting
% compute render matrix for different cone mosaic
function [mosaicArray, renderArray] = computeRender(ratio, retina, imageSize)
retina.resetSCone();
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
end

function [mosaicArray, renderArray] = computeRenderDichroma(ratio, retina, imageSize, flag)
retina.resetSCone();
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    % manipulate the number of S cone in the mosaic
    retina.resetCone();
    retina.resetSCone();
    retina.reassignSCone(ratio(idx));
    
    switch (flag)
        case 'noLCone'
            retina.reassignLCone(0.0, false);
        case 'noMCone'
            retina.reassignMCone(0.0, false);
        otherwise
            error('Invalid Option');
    end    
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize, false, true, false);    
    renderArray(idx) = {double(renderMtx)};
end
end

% compute image reconstruction for each mosaic
function output = computeRecon(input, renderArray, prior, regPara, imageSize)
nImage = size(input, 1);
output = zeros([length(renderArray), nImage, imageSize]);
parfor i = 1:length(renderArray)
    render = renderArray{i};
    estimator = ...
        PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regPara, 4, imageSize);
    
    for j = 1:nImage
        image = input(j, :, :, :);
        resp  = render * image(:);
        
        nIter = 800;
        recon = estimator.estimate(resp, nIter, rand([prod(imageSize), 1]), true, 1.0, 'final');
        output(i, j, :, :, :) = recon;
    end
end
end

% compute reconstruct error and show reconstructed images
function plotResults(input, output, ratio, display, imageSize)

% compute RMSE
nImage = size(input, 1);
rmse = zeros([length(ratio), nImage]);
for i = 1:length(ratio)
    for j = 1:nImage
        inputImage  = input(j, :, :, :);
        outputImage = output(i, j, :, :, :);
        rmse(i, j) = norm(inputImage(:) - outputImage(:));
    end
end

% plot RMSE
figure();
errorbar(ratio, mean(rmse, 2), std(rmse, 0, 2) / sqrt(nImage), '--ok', 'LineWidth', 2);
grid off; box off; hold on;

% show original images
figure();
plotAxis = tight_subplot(length(ratio) + 1, nImage, [.01 .01], [.01 .01], [.01 .01]);
for idx = 1:nImage
    image = reshape(input(idx, :, :, :), imageSize);
    image = gammaCorrection(image, display);
    axes(plotAxis(idx));
    imshow(image, 'InitialMagnification', 200);
end

% show reconstructed images
for i = 1:length(ratio)
    for j = 1:nImage
        image = reshape(output(i, j, :, :, :), imageSize);
        image = gammaCorrection(image, display);
        axes(plotAxis(i * nImage + j));
        imshow(image, 'InitialMagnification', 200);
    end
end
end
