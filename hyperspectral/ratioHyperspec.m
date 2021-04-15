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
meanLevel = mean(samples(:));

% Histogram equalization
[f, x] = ecdf(samples(:)); [x, ia] = unique(x); f = f(ia);

equalized = interp1(x, f, samples(:));
equalized = reshape(equalized, size(samples));

%% Show a few sample
scene = sceneCreate('whitenoise');
scene.spectrum.wave = wave;
figure();

nShow = 5e2;
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

%% Create mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', displayCreate('CRT12BitDisplay'));

retina.visualizeMosaic();

% Test Run
allCone = retina.computeHyperspectral...
    (wave, meanLevel, reshape(equalized(1, :, :, :), imageSize));
retina.visualizeExcitation();

%% Compute render matrix
renderMtx = retina.hyperRender(imageSize, wave, meanLevel, true);
renderMtx = double(renderMtx);

%% Validate render matrix
input = equalized(1, :, :, :);
coneVec = renderMtx * input(:);

scatter(coneVec, allCone);
axis equal; axis square;
