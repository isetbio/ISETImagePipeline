%% One-dimensional signal
dimension = 36;
rho = 0.9;

mu = 0.5 * ones(1, dimension);
mtx = covMtxMarkov(dimension, 0.3, rho);
sample = mvnrnd(mu, mtx, 1);

figure();
plot(sample, '-k');
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

%% Two-dimensional
spatial = kron(mtx, mtx);
mu = 0.5 * ones(1, dimension * dimension);
sample = mvnrnd(mu, spatial, 1);

figure();
imshow(reshape(sample, [dimension, dimension]), 'InitialMagnification', 800);

%% Color channel
mtx = covMtxMarkov(3, 0.25, rho);
chromatic = kron(mtx, spatial);

mu = 0.5 * ones(1, dimension * dimension * 3);
sample = mvnrnd(mu, chromatic, 1);

figure();
imshow(reshape(sample, [dimension, dimension, 3]), 'InitialMagnification', 800);

%% Load display file
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

%% Retinal mosaic
imageSize = [36, 36, 3];
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.4, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

%% Cone ratio manipulation
retina.resetCone();
retina.resetSCone();
retina.reassignSCone(0.05);
[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

%% Manipulation of M cone ratio
ratio = [0.3, 0.2, 0.15, 0.10, 0.05, 0.01, 0.001, 0];
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('M Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignMCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize);
    renderArray(idx) = {renderMtx};
end

%% Cone ratio manipulation
retina.resetCone();
retina.resetSCone();
retina.reassignSCone(0.05);
[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

%% Manipulation of L cone ratio
ratio = [0.6, 0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0];
mosaicArray_L = cell(1, length(ratio));
renderArray_L = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('L Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignLCone(ratio(idx));
    mosaicArray_L(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize);
    renderArray_L(idx) = {renderMtx};
end

%% All mosaic and render matrix
allMosaic = cat(2, flip(mosaicArray_L), mosaicArray);
allRender = cat(2, flip(renderArray_L), renderArray);

for idx = 1:length(allMosaic)
    retina.Mosaic.pattern = allMosaic{idx};
    retina.visualizeMosaic();
end

%% Gaussian model of signal
dimension = 36;
[mu, covMtx, regBasis] = colorSignal(dimension, 0.9, 0.9);

%% Sample from distribution
sample = mvnrnd(mu, covMtx, 1);

figure();
subplot(1, 2, 1);
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

sample = regBasis * normrnd(0, 1, [prod(imageSize), 1]) + 0.5;
subplot(1, 2, 2);
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

%% Test reconstruction
mosaicIdx = 8;
render = double(allRender{mosaicIdx});

estimator = RidgeGaussianEstimator(render, regBasis, mu');
estimator.setLambda(1e-3);

recon = (estimator.estimate((render * sample)'))';

figure();
subplot(1, 2, 1);
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

subplot(1, 2, 2);
imshow(reshape(recon, imageSize), 'InitialMagnification', 800);

%% Helper functions
function [mu, covMtx, regBasis] = colorSignal(dimension, rhoS, rhoC)

mtx = covMtxMarkov(dimension, 0.3, rhoS);
spatial = kron(mtx, mtx);

mtx = covMtxMarkov(3, 0.2, rhoC);
covMtx = kron(mtx, spatial);

mu = 0.5 * ones(1, dimension * dimension * 3);

[U, S, ~] = svd(covMtx);
visualizeBasis(U, dimension, dimension^2 * 3, false);

regBasis = U * diag(sqrt(diag(S)));
end

function mtx = covMtxMarkov(dimension, var, rho)

mtx = zeros(dimension, dimension);
for i = 1:dimension
    for j = 1:dimension
        mtx(i, j) = rho ^ abs(i - j);
    end
end
mtx = mtx * var;

end
