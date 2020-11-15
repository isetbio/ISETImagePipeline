%% Setup a cone mosaic & generate its corresponding render matrix
fovDegs = 0.50;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', fovDegs, 'integrationTime', 1.0);
retina.visualizeMosaic();

imageSize = [50, 50, 3];
render = retina.forwardRender(imageSize);
render = double(render);

%% Setup parameter & create scene
spatialFreq = 4;

% rmsContrast = 0.1;
% chromaDir = [1.0, 1.0, 0.0]';

rmsContrast = 0.1;
chromaDir = [1.0, -1.0, 0.0]';

chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

gratingScene = createGratingScene(chromaDir, spatialFreq, 'fovDegs', fovDegs);

crst = 1.0;
[theSceneSequence, temporalSupport] = gratingScene.compute(crst);

figure();
gratingScene.visualizeStaticFrame(theSceneSequence);

% Compute response from scene
[~, allCone] = retina.computeWithScene(theSceneSequence{:});

%% Reconstruction
load('../sparsePrior.mat');
estimator = PoissonSparseEstimator(double(render), inv(regBasis), mu', 0.1, 2, imageSize);

%% Single reconstruction
reconTest = estimator.estimate(allCone, 250, rand([prod(imageSize), 1]), true);

figure();
imshow(reshape(reconTest, imageSize), 'InitialMagnification', 1000);

%% Reconstruct Images
images = zeros(10, 50, 50, 3);

% achromatic
figure();
spatialFreq = [2, 4, 8, 12, 15];

rmsContrast = 0.1;
chromaDir = [1.0, 1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

crst = 0.50;
for idx = 1:length(spatialFreq)
    subplot(2, 5, idx);
    gratingScene = createGratingScene(chromaDir, spatialFreq(idx), 'fovDegs', fovDegs);
    
    [theSceneSequence, temporalSupport] = gratingScene.compute(crst);
    
    [~, allCone] = retina.computeWithScene(theSceneSequence{:});
    
    reconTest = estimator.estimate(allCone, 200, rand([prod(imageSize), 1]), true);
    images(idx, :, :, :) = reconTest;
    
    imshow(reshape(reconTest, imageSize), 'InitialMagnification', 1000);
end

% red-green
rmsContrast = 0.1;
chromaDir = [1.0, -1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

for idx = 1:length(spatialFreq)
    subplot(2, 5, idx + 5);
    gratingScene = createGratingScene(chromaDir, spatialFreq(idx), 'fovDegs', fovDegs);
    
    [theSceneSequence, temporalSupport] = gratingScene.compute(crst);
    
    [~, allCone] = retina.computeWithScene(theSceneSequence{:});
    
    reconTest = estimator.estimate(allCone, 200, rand([prod(imageSize), 1]), true);
    images(idx + 5, :, :, :) = reconTest;
    
    imshow(reshape(reconTest, imageSize), 'InitialMagnification', 1000);
end

%% images
image1 = mean(reshape(images(5, :, :, :), imageSize), 3);
image2 = mean(reshape(images(10, :, :, :), imageSize), 3);

figure();
plot(1:50, image1(25, :), '-k'); hold on;
plot(1:50, image2(25, :), '--k'); hold on;
