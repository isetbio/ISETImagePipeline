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
mtx = MarkovPrior.covMtxMarkov(3, 0.3, rho);
chromatic = kron(mtx, spatial);

mu = 0.5 * ones(1, dimension * dimension * 3);
sample = mvnrnd(mu, chromatic, 1);

figure();
imshow(reshape(sample, [dimension, dimension, 3]), 'InitialMagnification', 800);

%% Load display file & Generate retinal mosaic
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

imageSize = [36, 36, 3];
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

retina.visualizeMosaic();

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
imageSize = [36, 36, 3];

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

%% Run analysis (with different spatial and chromatic correlation)
nDim = 36;
corrSpatial = 0.9;
corrChromat = 0.9;
nRecon = 10;

errorMtx = MarkovPrior.reconFunc(allRender, nDim, corrSpatial, corrChromat, nRecon, false);

%% Plot curve
errorMean = mean(errorMtx, 2);
errorSD   = std(errorMtx, 0, 2);

figure();
errorbar(allRatio, errorMean, errorSD / sqrt(nRecon) * 2, '-ok', 'LineWidth', 1);
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

xlabel('L Cone Ratio');
ylabel('Reconstruction RMSE');

%% Different parameter for spatial and chromatic correlation
spatial = [0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1];
chromat = spatial;

nDim = 36;
nRecon  = 10;

nCorr = length(spatial);
allError = zeros(nCorr, nCorr, length(allRender), nRecon);

for i = 1:length(spatial)
    for j = 1:length(chromat)
        corrSpatial = spatial(i);
        corrChromat = chromat(j);
        errorMtx = MarkovPrior.reconFunc(allRender, nDim, corrSpatial, corrChromat, nRecon, false);

        allError(i, j, :, :) = errorMtx;
        fprintf('Recon Set ID: %d, %d; Total: %d * %d \n', i, j, nCorr, nCorr);
    end
end

%% Plot results
% load priorEffect.mat & retinaRenderTiny.mat
nCorr = length(spatial);
exclude = [3, 9];
plotAxis = tight_subplot(nCorr - length(exclude), nCorr - length(exclude), [.05 .05], 0.05, 0.05);
plotIdx = 1;
limits = [1.25, 4.5; 1.25, 4.5; 1, 4.5; 1, 4.5; 1, 4.5; 1, 4; 0.5, 3.75; 0.25, 3.25; 0.25, 3.25];

matchedY = false;
for i = 1:nCorr
    for j = 1:nCorr
        if(sum(i == exclude) == 0 && sum(j == exclude) == 0)
            [~, ~, nMosaic, nRecon] = size(allError);
            errorMtx = reshape(allError(i, j, :, :), [nMosaic, nRecon]);
            errorMean = mean(errorMtx, 2);
            errorSD   = std(errorMtx, 0, 2);

            axes(plotAxis(plotIdx));
            plotIdx = plotIdx + 1;

            minError = min(errorMean);
            marginIdx = allRatio((errorMean <= minError + 0.1));
            
            % errorbar(allRatio, errorMean, errorSD / sqrt(nRecon) * 2, '-k', 'LineWidth', 1.5);
            plot(allRatio, errorMean, '-k', 'LineWidth', 1.5); hold on;
            scatter(allRatio, errorMean, 40, 'k', 'filled', 'LineWidth', 1.0);

            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');

            labelIDX = [4, 6:8, 10, 13];
            xticks(allRatio(labelIDX));

            if matchedY
                % Matched Y-axis
                ylim([0, 4.5]);
            else
                % Scaled Y-axis
                yaxisLim = ylim();
                yaxisLim = [floor(yaxisLim(1) * 2) / 2, ceil(yaxisLim(2) * 2) / 2];
                ylim(yaxisLim);
                ytickPos = yaxisLim(1) : 0.5 : yaxisLim(2);
                if length(ytickPos) >= 5
                    ytickPos = yaxisLim(1) : 1.0 : yaxisLim(2);
                end
                yticks(ytickPos);

            end

            yaxisLim = ylim();
            shade = area([marginIdx(1), marginIdx(end)], [yaxisLim(end), yaxisLim(end)], ...
                'FaceColor', ones(1, 3) * 0.8, 'LineStyle','none'); hold on;
            uistack(shade,'bottom');

            set(gca, 'linewidth', 0.75)
            set(gca, 'TickLength', [0.03, 0.025])

            if(~(i == 8))
                xticklabels([]);
            else
                xticklabels(allRatio(labelIDX) * 100);
            end

        end
    end
end

%% Add noise
nDim = 36;
corrSpatial = 0.9;
corrChromat = 0.9;
nRecon = 10;

errorMtxNoise = MarkovPrior.reconFuncNoise(allRender, nDim, corrSpatial, corrChromat, nRecon, false, 1e3, 0.01);

%% Plotting function
figure();
plotError(errorMtx, allRatio, nRecon);
plotError(errorMtxNoise, allRatio, nRecon);

%% plotting helper function
function plotError(errorMtx, allRatio, nRecon)

errorMean = mean(errorMtx, 2);
errorSD   = std(errorMtx, 0, 2);

hold on;
errorbar(allRatio, errorMean, errorSD / sqrt(nRecon) * 2, '-ok', 'LineWidth', 1);
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

xlabel('L Cone Ratio');
ylabel('Reconstruction RMSE');

end
