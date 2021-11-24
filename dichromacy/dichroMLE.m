%% Generate human cone mosaic with optics
imageSize = [100, 100, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display, 'pupilSize', 2.5);

%% Change all M cone to L cone
retina.reassignCone(0.0, retina.M_Cone_Idx, retina.L_Cone_Idx, false);

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

retina.visualizeMosaic();

%% Change all L cone to M cone
retina.reassignCone(0.0, retina.L_Cone_Idx, retina.M_Cone_Idx, false);

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

retina.visualizeMosaic();

%% No S cone condition
retina.resetSCone();

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

retina.visualizeMosaic();

%% Render matrix
render = retina.forwardRender(imageSize);
render = double(render);

%% Reconstruction
load('./dichromacy/inputImage_100.mat');

% Build an image reconstruction object with sparse prior
regConst = 0.0; stride = 4;
estimator = ...
    PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regConst, stride, imageSize);

% Run reconstruction on cone response to each images
% reconstructed images are in linear pixel space, need to
% gamma correct them before visulization
nIter = 800; optDisp = 'final';
output = zeros(size(inputLinear));

for idx = 1 : size(inputLinear, 1)
    input = reshape(inputLinear(idx, :, :, :), imageSize);
    coneResp = render * input(:);
    recon = estimator.runEstimate(coneResp, 'maxIter', nIter, 'display', optDisp);
    output(idx, :, :, :) = gammaCorrection(recon, display);
end

%% Show results
figure();
for idx = 1 : size(output, 1)
    subplot(2, 5, idx);
    reconImage = reshape(output(idx, :, :, :), imageSize);
    imshow(reconImage, 'InitialMagnification', 200);
end

%% Understand the solution
imageID = 10;

inTest = reshape(inputLinear(imageID, :, :, :), imageSize);
outTest = reshape(output(imageID, :, :, :), imageSize);

[~, ~, outTest, ~] = retina.compute(outTest);

%% Check cone excitation
figure();
subplot(1, 2, 1); imshow(inTest);
subplot(1, 2, 2); imshow(outTest);

inTest(inTest >= 1) = 1 - 1e-10;
inTest(inTest <= 0) = 0 + 1e-10;

figure();
scatter(render * inTest(:), render * outTest(:));
axis equal; xlim([0, 2000]); ylim([0, 2000]);

%% Run reconstruction with original as input
regConst = 0.0; stride = 4;
estimator = PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regConst, stride, imageSize);

coneResp = render * inTest(:);
recon = estimator.runEstimate(coneResp, 'maxIter', 100, 'display', 'iter', 'init', inTest(:));
recon = estimator.runEstimate(coneResp, 'maxIter', 100, 'display', 'iter', 'init', outTest(:));

%% Intuition with null space
basisMtx = null(render);

%% Rand image from null space
image = zeros(imageSize);
image(:, :, 1) = 1;

[vectRow, vectNull] = linearDecomp(render, basisMtx, image(:));
visualizeImage(vectRow, vectNull, imageSize);

image = zeros(imageSize);
image(:, :, 2) = 1;

[vectRow, vectNull] = linearDecomp(render, basisMtx, image(:));
visualizeImage(vectRow, vectNull, imageSize);

image = zeros(imageSize);
image(:, :, 3) = 1;

[vectRow, vectNull] = linearDecomp(render, basisMtx, image(:));
visualizeImage(vectRow, vectNull, imageSize);

%% Image Diff
figure();
subplot(1, 3, 1); imshow(inTest);
subplot(1, 3, 2); imshow(outTest);

diff = inTest - outTest;
subplot(1, 3, 3); imshow(diff + 0.5);

figure();
histogram(render * diff(:))

%% interpolation
for idx = 1:6
    interpolation = outTest + (idx - 1) * 0.2 * diff;

    figure(1);
    subplot(2, 3, idx);
    imshow(interpolation);

    figure(2);
    subplot(2, 3, idx);
    scatter(render * inTest(:), render * interpolation(:));
end

%% Optimization with an interpolation
intp = outTest + 0.5 * diff;
recon = estimator.runEstimate(render * inTest(:), 'maxIter', 100, ...
                                'display', 'iter', 'init', intp(:));

%% Show recon
figure();
subplot(1, 2, 1); imshow(intp);
subplot(1, 2, 2); imshow(recon);

%% Run til absolute convergence
recon = estimator.runEstimate(render * inTest(:), 'maxIter', 1e4, ...
                                'display', 'iter', 'init', outTest(:));

%% Error landscape analysis
r = 0:0.01:1; 
g = 0:0.01:1;
[R, G] = meshgrid(r, g);
llhdGrid = zeros(size(R));

input = zeros(imageSize);
input(:, :, 1) = 0.75; 
input(:, :, 2) = 0.25; 

coneVec = render * input(:);

for i = 1:101
    for j = 1:101
        input = zeros(imageSize);

        input(:, :, 1) = R(i, j);
        input(:, :, 2) = G(i, j);
    
        [nlogll, ~] = estimator.likelihood(coneVec, input(:));
        llhdGrid(i, j) = nlogll;
    end
end

recon_1 = estimator.runEstimate(coneVec, 'maxIter', 1e3, 'display', 'iter');
recon_2 = estimator.runEstimate(coneVec, 'maxIter', 1e3, 'display', 'iter');

%% Show plot 1
contour(R, G, llhdGrid, linspace(-3.48e+07, -3.45e+07, 20), 'ShowText', 'on'); 

%% Show recon
figure();
subplot(1, 3, 1);
input = zeros(imageSize);
input(:, :, 1) = 0.75; 
input(:, :, 2) = 0.25; 
imshow(input);

subplot(1, 3, 2);
imshow(recon_1);

subplot(1, 3, 3);
imshow(recon_2);

%% Show plot 2
contour(R, G, llhdGrid, linspace(-6.81e+07, -6.7e+07, 10), 'ShowText', 'on'); 

%% null space
[vectRow, vectNull] = linearDecomp(render, basisMtx, input(:));
visualizeImage(vectRow, vectNull, imageSize);

%% error landspace with the prior
% Build an image reconstruction object with sparse prior
regConst = 0.01; stride = 4;
estimator = ...
    PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regConst, stride, imageSize);

r = 0:0.025:1; 
g = 0:0.025:1;
[R, G] = meshgrid(r, g);
llhdGrid = zeros(size(R));

input = zeros(imageSize);
input(:, :, 1) = 0.75; 
input(:, :, 2) = 0.25; 

coneVec = render * input(:);

for i = 1:41
    for j = 1:41
        input = zeros(imageSize);

        input(:, :, 1) = R(i, j);
        input(:, :, 2) = G(i, j);
    
        [nlogll, ~] = estimator.reconObjective(coneVec, input(:));
        llhdGrid(i, j) = nlogll;
    end
end

recon_1 = estimator.runEstimate(coneVec, 'maxIter', 800, 'display', 'iter');
recon_2 = estimator.runEstimate(coneVec, 'maxIter', 800, 'display', 'iter');

%% adding prior
contour(R, G, llhdGrid, linspace(-6.8109e+07, -6.81e+07, 10)); 

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
