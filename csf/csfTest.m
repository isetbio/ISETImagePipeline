%% Setup a cone mosaic & generate its corresponding render matrix
fovDegs = 0.50;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', fovDegs, 'integrationTime', 1.0);
retina.visualizeMosaic();

imageSize = [50, 50, 3];
render = retina.forwardRender(imageSize);
render = double(render);

%% Setup parameter & create scene
spatialFreq = 2;

rmsContrast = 0.1;
chromaDir = [1.0, 1.0, 1.0]';

% rmsContrast = 0.1;
% chromaDir = [1.0, -1.0, 0.0]';

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
estimator = PoissonSparseEstimator(double(render), inv(regBasis), mu', 0.06, 2, imageSize);
reconTest = estimator.estimate(poissrnd(allCone), 250, rand([prod(imageSize), 1]), true);

figure();
imshow(reshape(reconTest, imageSize), 'InitialMagnification', 1000);

%% Reconstruct Images
% achromatic
figure();
spatialFreq = [4, 8, 12, 16, 24];

rmsContrast = 0.1;
chromaDir = [1.0, 1.0, 1.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

crst = 0.50;
for idx = 1:length(spatialFreq)
    subplot(2, 5, idx);
    gratingScene = createGratingScene(chromaDir, spatialFreq(idx), 'fovDegs', fovDegs);
    
    [theSceneSequence, temporalSupport] = gratingScene.compute(crst);
    
    [~, allCone] = retina.computeWithScene(theSceneSequence{:});
    
    reconTest = estimator.estimate(allCone, 250, rand([prod(imageSize), 1]), true);
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
    
    reconTest = estimator.estimate(allCone, 250, rand([prod(imageSize), 1]), true);
    imshow(reshape(reconTest, imageSize), 'InitialMagnification', 1000);
end
