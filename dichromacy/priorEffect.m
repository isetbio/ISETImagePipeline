%% One-dimensional signal
dimension = 32;
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
mtx = covMtxMarkov(3, 0.25, 0.9);
chromatic = kron(mtx, spatial);

mu = 0.5 * ones(1, dimension * dimension * 3);
sample = mvnrnd(mu, chromatic, 1);

figure();
imshow(reshape(sample, [dimension, dimension, 3]), 'InitialMagnification', 800);

%% Load display file
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

%% Retinal mosaic
imageSize = [32, 32, 3];
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
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));

for idx = 1:length(ratio)
    fprintf('L Cone Ratio: %.4f \n', ratio(idx));
    retina.reassignLCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize);
    renderArray(idx) = {renderMtx};
end

%% Helper functions
function mtx = covMtxMarkov(dimension, var, rho)
    mtx = zeros(dimension, dimension);
    for i = 1:dimension
        for j = 1:dimension
            mtx(i, j) = rho ^ abs(i - j);
        end
    end
    mtx = mtx * var;
end
