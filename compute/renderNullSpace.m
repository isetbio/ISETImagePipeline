%% Null space of render matrix
basisMtx = null(render);

%% Check
for n = 1:10
    idx = randi(size(basisMtx, 2), 1);
    basisVec = basisMtx(:, idx);
    
    proj = norm(render * basisVec);
    fprintf('Vector norm: %.2f \n', proj);
end

%% Visualization
[~] = visualizeBasis(basisMtx, 140, size(basisMtx, 2), false);

%% Rand image
% weight = 0.5 * rand(1, size(basisMtx, 2))';
% weight = normrnd(0.2, 0.1, [1, size(basisMtx, 2)])';
nDim = size(basisMtx, 2);
nIdx = floor(nDim / 4);
weight = [zeros(1, nDim - nIdx), rand(1, nIdx)];
weight = weight(randperm(length(weight)))';

image  = basisMtx * weight;

figure(); subplot(1, 2, 1);
imshow(reshape(image, [140, 140, 3]), 'InitialMagnification', 500);
fprintf('Vector norm: %.4f \n', norm(render * image));

% subplot(1, 3, 2);
% scatter(render*patchLinear(:), render*(image + patchLinear(:)));

subplot(1, 2, 2);
sumImage = reshape(patchLinear(:) + image, [140, 140, 3]);
imshow(invGammaCorrection(sumImage, display.CRT12BitDisplay), 'InitialMagnification', 500);

%% Analysis
% Projecting images/stimulus to Row Space and the Null Space of Render Matrix
% Sample test image
projectName  = 'ISETImagePipeline';
thisImageSet = 'ILSVRC';
dataBaseDir  = getpref(projectName, 'dataDir');
imageName    = 'ILSVRC2017_test_00000025.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.5);

image = sampleImage(image, 100);
image(image < 0) = 0;
image(image > 1) = 1;
imshow(image, 'InitialMagnification', 500);

imageSize = [140, 140, 3];

%% Analysis with Projection onto Null Space
[vectRow, vectNull] = linearDecomp(renderMatrix, basisMtx, image(:));

% Visualization
visualizeImage(vectRow, vectNull, imageSize);

%% Repeat the analysis with two line stimulus
addpath('./twoline/');
stimConfig = 1;

projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

for idx = 1:5
    [~, img] = generateTwoLineScene(display.CRT12BitDisplay, stimConfig, 2 ^ idx);
    [vectRow, vectNull] = linearDecomp(renderMatrix, basisMtx, img(:));
    visualizeImage(vectRow, vectNull, imageSize);
end

% Fisher information matrix analysis
%% Gaussian noise assumption
nVec = 121;
FIM = renderMatrix' * renderMatrix;

disp('Compute largest eigenvectors');
[Vlarge, Dlarge] = eigs(FIM, nVec, 'largestabs');

disp('Compute smallest eigenvectors');
[Vsmall, Dsmall] = eigs(FIM, nVec, 'smallestabs');

%% Show direction vector
for idx = 1:nVec
    figure(1);
    subplot(sqrt(nVec), sqrt(nVec), idx);
    imshow(reshape(rescaleBasis(Vlarge(:, idx)), imageSize), 'InitialMagnification', 500);
    
    figure(2);
    subplot(sqrt(nVec), sqrt(nVec), idx);
    imshow(reshape(rescaleBasis(Vsmall(:, idx)), imageSize), 'InitialMagnification', 500);
end

% Poisson noise assumption
%% Calculate FIM for poisson noise
FIM = poissonFIM(renderMatrix, image(:));

%% Eigendecompostion
nVec = 121;

disp('Compute largest eigenvectors');
[Vlarge, Dlarge] = eigs(FIM, nVec, 'largestabs');

disp('Compute smallest eigenvectors');
[Vsmall, Dsmall] = eigs(FIM, nVec, 'smallestabs');

%% Show direction vector
for idx = 1:nVec
    figure(1);
    subplot(sqrt(nVec), sqrt(nVec), idx);
    imshow(reshape(rescaleBasis(Vlarge(:, idx)), imageSize), 'InitialMagnification', 500);
    
    figure(2);
    subplot(sqrt(nVec), sqrt(nVec), idx);
    imshow(reshape(rescaleBasis(Vsmall(:, idx)), imageSize), 'InitialMagnification', 500);
end

%% Helper function
function [vectRow, vectNull] = linearDecomp(render, nullBasis, imageVec)
% Projection onto Null Space
coeffNull = nullBasis' * imageVec;
vectNull  = nullBasis * coeffNull;

% Orthogonal Complement
vectRow = imageVec - vectNull;

% Check
fprintf('Vector norm for row image: %.4f \n', norm(render * vectRow));
fprintf('Vector norm for null image: %.4f \n', norm(render * vectNull));
fprintf('Vector dot product: %.4f \n', vectRow' * vectNull);

end

function visualizeImage(vectRow, vectNull, imageSize)
figure();
subplot(1, 3, 1);
imshow(reshape(vectRow + vectNull, imageSize), 'InitialMagnification', 500);
title('Image');

subplot(1, 3, 2);
imshow(reshape(vectRow, imageSize), 'InitialMagnification', 500);
title('Row Space');

subplot(1, 3, 3);
imshow(reshape(vectNull, imageSize), 'InitialMagnification', 500);
title('Null Space');
end

function FIM = poissonFIM(render, image)
nResponse = size(render, 1);
FIM = zeros(length(image), length(image));
for idx = 1:nResponse
    subRender = render(idx, :);
    lambda = subRender * image;
    FIM = FIM + subRender' * subRender / lambda;
end
end

function basisVec = rescaleBasis(basisVec)
basisVec = (basisVec - min(basisVec));
basisVec = basisVec ./ max(basisVec);
end
