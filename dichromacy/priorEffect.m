%% One-dimensional signal
dimension = 36;
rho = 0.9;

mu = 0.5 * ones(1, dimension);
mtx = MarkovPrior.covMtxMarkov(dimension, 0.3, rho);
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
mtx = MarkovPrior.covMtxMarkov(3, 0.25, rho);
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
allRatio  = cat(2, flip([0.6, 0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0]), ...
    1 - [0.3, 0.2, 0.15, 0.10, 0.05, 0.01, 0.001, 0]);

for idx = 1:length(allMosaic)
    retina.Mosaic.pattern = allMosaic{idx};
    retina.visualizeMosaic();
end

%% Gaussian model of signal
dimension = 36;

corrSpatial = 0.9;
corrChromat = 0.9;
[mu, covMtx, regBasis] = MarkovPrior.colorSignal(dimension, corrSpatial, corrChromat);

%% Sample from distribution
sample = mvnrnd(mu, covMtx, 1);

figure();
subplot(1, 2, 1);
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

sample = regBasis * normrnd(0, 1, [prod(imageSize), 1]) + mu';
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

%% Run analysis
nDim = 36;
corrSpatial = 0.0;
corrChromat = 0.0;
nRecon = 20;

errorMtx = MarkovPrior.reconFunc(allRender, nDim, corrSpatial, corrChromat, nRecon, false);

%% Show signal
imageSize = [nDim, nDim, 3];
[mu, ~, regBasis] = MarkovPrior.colorSignal(nDim, corrSpatial, corrChromat, true);
sample = regBasis * normrnd(0, 1, [prod(imageSize), 1]) + mu';

figure();
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

%% Plot curve
errorMean = mean(errorMtx, 2);
errorSD   = std(errorMtx, 0, 2);

figure();
errorbar(allRatio, errorMean, errorSD, '-ok', 'LineWidth', 1);
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

xlabel('L Cone Ratio');
ylabel('Reconstruction RMSE');

%% Different parameter for spatial and chromatic correlation
spatial = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1];
chromat = spatial;

nDim = 36;
nRecon  = 20;

nCorr = length(spatial);
allError = zeros(nCorr, nCorr, length(allRender), nRecon);

for i = 1:length(spatial)
    for j = 1:length(chromat)
        corrSpatial = spatial(i);
        corrChromat = chromat(j);
        errorMtx = MarkovPrior.reconFunc(allRender, nDim, corrSpatial, corrChromat, nRecon, false);
        
        allError(i, j, :, :) = errorMtx;
        printf('Recon Set ID: %d, %d \n', i, j);
    end
end
