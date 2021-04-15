%% Read scene from database
[manData, ~] = manchesterDB();
[harData, ~] = harvardDB();

%% Sample image patch
wave = manData.wave;
assert(sum(wave ~= harData.wave) == 0);

imageSize = [18, 18, length(wave)];
edge = imageSize(1);
nImage = 1e4;

count = 1;
samples = zeros([nImage, imageSize]);

for idx = 1 : length(manData.image)
    image = manData.image{idx};
    [h, w, d] = size(image);
    
    assert(d == length(wave));
    nSample = floor(sqrt(h * w / (edge ^ 2)));
    
    for idj = 1 : nSample
        hIdx = randi(h - edge + 1);
        wIdx = randi(w - edge + 1);
        samples(count, :, :, :) = ...
            image(hIdx : (hIdx + edge - 1), wIdx : (wIdx + edge - 1), :);
        
        count = count + 1;
    end
end

for idx = 1 : length(harData.image)
    image = harData.image{idx};
    [h, w, d] = size(image);
    
    assert(d == length(wave));
    nSample = floor(sqrt(h * w / (edge ^ 2)));
    
    for idj = 1 : nSample
        hIdx = randi(h - edge + 1);
        wIdx = randi(w - edge + 1);
        samples(count, :, :, :) = ...
            image(hIdx : (hIdx + edge - 1), wIdx : (wIdx + edge - 1), :);
        
        count = count + 1;
    end
end

samples = samples(1:count-1, :, :, :);
meanLevel = mean(samples(:)) * 2;

% Histogram equalization
[f, x] = ecdf(samples(:)); [x, ia] = unique(x); f = f(ia);

equalized = interp1(x, f, samples(:));
equalized = reshape(equalized, size(samples));

%% Show a few sample
scene = sceneCreate('whitenoise');
scene.spectrum.wave = wave;
figure();

nShow = 50;
for idx = 1 : nShow
    image = equalized(randi(count - 1), :, : , :);
    image = reshape(image, imageSize);
    
    scene.data.photons = image * meanLevel;
    imshow(sceneGet(scene, 'rgbimage'), 'InitialMagnification', 1e3);
end

%% PCA Analysis
imgData = reshape(equalized, [size(samples, 1), prod(imageSize)]);

[pcaBasis, ~, pcaVar, ~, ~, mu] = pca(imgData);
scaleMatrix = diag(sqrt(pcaVar));
regBasis    = pcaBasis * scaleMatrix;

%% I. Test run with basic reconstruction
% Create mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.25, 'display', displayCreate('CRT12BitDisplay'));

retina.visualizeMosaic();

%% Test Run
imageID = randi(count - 1);
allCone = retina.computeHyperspectral...
    (wave, meanLevel, reshape(equalized(imageID, :, :, :), imageSize));
retina.visualizeExcitation();

%% Compute render matrix
renderMtx = retina.hyperRender(imageSize, wave, meanLevel, true);
renderMtx = double(renderMtx);

%% Validate render matrix
input = equalized(imageID, :, :, :);
coneVec = renderMtx * input(:);

scatter(coneVec, allCone);
axis equal; axis square;

%% Reconstruction routine
estimator = RidgeGaussianEstimator(renderMtx, regBasis, mu');
estimator.setLambda(1e-3);

recon = estimator.estimate(coneVec');

%% Show results
scene = sceneCreate('whitenoise');
scene.spectrum.wave = wave;

figure(); subplot(1, 2, 1);
scene.data.photons = reshape(input, imageSize) * meanLevel;
imshow(sceneGet(scene, 'rgbimage'), 'InitialMagnification', 1e3);

scene.data.photons = reshape(input, imageSize) * meanLevel;

subplot(1, 2, 2);
scene.data.photons = reshape(recon, imageSize) * meanLevel;
imshow(sceneGet(scene, 'rgbimage'), 'InitialMagnification', 1e3);

%% Cross validation
nTest = 20;
imageID = randi(count - 1, [50, 1]);
testInput = equalized(imageID, :, :, :);

regPara = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10];
rss = zeros(size(regPara));

estimator = RidgeGaussianEstimator(renderMtx, regBasis, mu');
for idx = 1:length(regPara)
    estimator.setLambda(regPara(idx));
    
    for idj = 1:nTest
        input = testInput(idj, :);
        cone  = input * renderMtx';
        output = estimator.estimate(cone);
        
        rss(idx) = rss(idx) + norm(input - output);
    end
end

figure();
plot(log10(regPara), rss, '-ok');

%% II. Render matrix for different cone mosaic
% Create mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.25, 'display', displayCreate('CRT12BitDisplay'));

% Cone ratio manipulation
retina.resetCone();
retina.resetSCone();
retina.reassignSCone(0.05);
[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

% Manipulation of M cone ratio
ratio = [0.3, 0.2, 0.15, 0.10, 0.05, 0.01, 0.001, 0];
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('M Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignMCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.hyperRender(imageSize, wave, meanLevel, false);
    renderArray(idx) = {renderMtx};
end

% Reset cone ratio
retina.resetCone();
retina.resetSCone();
retina.reassignSCone(0.05);
[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

% Manipulation of L cone ratio
ratio = [0.6, 0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0];
mosaicArray_L = cell(1, length(ratio));
renderArray_L = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('L Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignLCone(ratio(idx));
    mosaicArray_L(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.hyperRender(imageSize, wave, meanLevel, false);
    renderArray_L(idx) = {renderMtx};
end

% All mosaic and render matrix
allMosaic = cat(2, flip(mosaicArray_L), mosaicArray);
allRender = cat(2, flip(renderArray_L), renderArray);
allRatio  = cat(2, flip([0.6, 0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0]), ...
    1 - [0.3, 0.2, 0.15, 0.10, 0.05, 0.01, 0.001, 0]);

for idx = 1:length(allMosaic)
    retina.Mosaic.pattern = allMosaic{idx};
    retina.visualizeMosaic();
end

%% III. RSS as function of cone ratio
nMosaic = length(allRatio);
nImage  = 100;

imageID = randi(count - 1, [nImage, 1]);
inputSet = equalized(imageID, :, :, :);

rss = zeros(nMosaic, nImage);

for idx = 1 : nMosaic
    fprintf('Run ID %d \n', idx);
    render = double(allRender{idx});
    
    estimator = RidgeGaussianEstimator(render, regBasis, mu');
    estimator.setLambda(1e-2);
    
    parfor idj = 1 : nImage
        input  = inputSet(idj, :);
        cone   = input * render';
        output = estimator.estimate(cone);
        
        rss(idx, idj) = norm(input - output);
    end
end

%% Plot RSS and SEM
% Figure format
try 
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
        'figureWidthInches', 15, ...
        'figureHeightInches', 10);
catch EXP
    fprintf('plotlab not available, use default MATLAB style \n');
end

%% Plot
plotIdx = [1, 5:10, 12, 16];

figure();
meanRSS = mean(rss, 2);
stdRSS  = std(rss, 0, 2);

errorbar(allRatio(plotIdx), meanRSS(plotIdx), ...
    stdRSS(plotIdx) / sqrt(nImage), '-ok', 'LineWidth', 2);

box off; grid off;
yticks(2 : 0.5 : 4);
ylim([2, 4]);

xlabel('L Cone Ratio'); ylabel('RSS, Hyperspectral');

%% IV. S Cone ratio
ratio = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95];
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.50, 'display', displayCreate('CRT12BitDisplay'));

retina.resetSCone();
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    % Manipulate the number of S cone in the mosaic
    retina.resetSCone();
    retina.reassignSCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.hyperRender(imageSize, wave, meanLevel, false);    
    renderArray(idx) = {double(renderMtx)};
end

%% Run reconstruction
nMosaic = length(ratio);
nImage  = 100;

imageID = randi(count - 1, [nImage, 1]);
inputSet = equalized(imageID, :, :, :);

rss = zeros(nMosaic, nImage);

for idx = 1 : nMosaic
    fprintf('Run ID %d \n', idx);
    render = double(renderArray{idx});
    
    estimator = RidgeGaussianEstimator(render, regBasis, mu');
    estimator.setLambda(1e-2);
    
    parfor idj = 1 : nImage
        input  = inputSet(idj, :);
        cone   = input * render';
        output = estimator.estimate(cone);
        
        rss(idx, idj) = norm(input - output);
    end
end

%% Plot RSS and SEM
% Figure format
try 
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
        'figureWidthInches', 15, ...
        'figureHeightInches', 10);
catch EXP
    fprintf('plotlab not available, use default MATLAB style \n');
end

%% Plot
plotIdx = 1 : nMosaic;

figure();
meanRSS = mean(rss, 2);
stdRSS  = std(rss, 0, 2);

errorbar(ratio(plotIdx), meanRSS(plotIdx), ...
    stdRSS(plotIdx) / sqrt(nImage), '-ok', 'LineWidth', 2);

box off; grid off;
yticks(2 : 0.5 : 4);
ylim([2, 4]);

xlabel('S Cone Ratio'); ylabel('RSS, Hyperspectral');
