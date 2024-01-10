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

%% 1. Example reconstruction with Gaussian signal
dimension = 36;
imageSize = [36, 36, 3];

corrSpatial = 0.9;
corrChromat = 0.9;
[mu, covMtx, regBasis] = MarkovPrior.colorSignal(dimension, corrSpatial, corrChromat);

%% Sample from distribution
figure();
sample = regBasis * normrnd(0, 1, [prod(imageSize), 1]) + mu';
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

%% Test reconstruction
mosaicIdx = 8;
render = double(allRender{mosaicIdx});

estimator = RidgeGaussianEstimator(render, regBasis, mu');
estimator.setLambda(1e-3);

recon = (estimator.estimate((render * sample)'))';

figure();
subplot(1, 2, 1);
title('Original');
imshow(reshape(sample, imageSize), 'InitialMagnification', 800);

subplot(1, 2, 2);
title('Reconstruction');
imshow(reshape(recon, imageSize), 'InitialMagnification', 800);

%% 2. Run analysis with a fixed spatial and chromatic correlation and different cone ratio
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

%% 3. Full Analysis with different parameter for spatial and chromatic correlation
% Note: This takes a long time to run
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
nCorr = length(spatial);
exclude = [3, 9];
plotAxis = tight_subplot(nCorr - length(exclude), nCorr - length(exclude), [.05 .05], 0.05, 0.05);
plotIdx = 1;
limits = [1.25, 4.5; 1.25, 4.5; 1, 4.5; 1, 4.5; 1, 4.5; 1, 4; 0.5, 3.75; 0.25, 3.25; 0.25, 3.25];
for i = 1:nCorr
    for j = 1:nCorr
        if(sum(i == exclude) == 0 && sum(j == exclude) == 0)
            [~, ~, nMosaic, nRecon] = size(allError);
            errorMtx = reshape(allError(i, j, :, :), [nMosaic, nRecon]);
            errorMean = mean(errorMtx, 2);
            errorSD   = std(errorMtx, 0, 2);
            
            axes(plotAxis(plotIdx));
            plotIdx = plotIdx + 1;
            
            % errorbar(allRatio, errorMean, errorSD / sqrt(nRecon) * 2, '-k', 'LineWidth', 1.5);
            plot(allRatio, errorMean, '-ok', 'LineWidth', 1.5);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            
            labelIDX = [4, 6:8, 10, 13];
            xticks(allRatio(labelIDX));
            
            % Matched Y-axis
            % ylim(limits(i, :));
            % yticks(ceil(limits(i, 1)) : 0.5 : floor(limits(i, 2)));
            
            % Scaled Y-axis
            yaxisLim = ylim();
            yaxisLim = [floor(yaxisLim(1) * 2) / 2, ceil(yaxisLim(2) * 2) / 2];            
            ylim(yaxisLim);
            ytickPos = yaxisLim(1) : 0.5 : yaxisLim(2);
            if length(ytickPos) >= 5
                ytickPos = yaxisLim(1) : 1.0 : yaxisLim(2);
            end
            yticks(ytickPos);
            
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
