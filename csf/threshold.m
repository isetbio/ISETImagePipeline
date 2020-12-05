%% Scene engine for static grating
spatialFreq = 4;

rmsContrast = 0.1;
chromaDir = [1.0, -1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

fovDegs = 0.50;
gratingScene = createGratingScene(chromaDir, spatialFreq, 'fovDegs', fovDegs);

%% Setup a cone mosaic & generate its corresponding render matrix
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'integrationTime', 1.0);

imageSize = [50, 50, 3];
render = retina.forwardRender(imageSize, false);
render = double(render);

%% Neural engine with ConeResponse class
load('../sparsePrior.mat');
estimator = PoissonSparseEstimator(double(render), inv(regBasis), mu', 0.1, 1, imageSize);

computeFunction = @(neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin) ...
    ...
    reconNeuralEngine(estimator, retina, neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin{:});

neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralEngine = neuralResponseEngine(computeFunction, neuralParams);

%% Recon Classifier engine
classifierEngine = responseClassifierEngine(@reconClassifier, struct());

trainFlag = 'none'; testFlag = 'random';
nTrain = 1; nTest = 48;

%% Test single stimulus
nullContrast = 0.0;
[nullScene, temporalSupport] = gratingScene.compute(nullContrast);

testContrast = 1.0;
[testScene, ~] = gratingScene.compute(testContrast);
gratingScene.visualizeStaticFrame(testScene);

[prediction, ~, responses] = computePerformanceTAFC...
        (nullScene, testScene, temporalSupport, nTrain, nTest, neuralEngine, classifierEngine, trainFlag, testFlag, true, false);

%% Compute performance
nullContrast = 0.0;
[nullScene, temporalSupport] = gratingScene.compute(nullContrast);

testContrast = [0.0, 0.25, 0.50, 0.75, 1.0];
pCorrect = zeros(size(testContrast));
for idx = 1:length(testContrast)
    
    [testScene, ~] = gratingScene.compute(testContrast(idx));
    
    [prediction, ~] = computePerformanceTAFC...
        (nullScene, testScene, temporalSupport, nTrain, nTest, neuralEngine, classifierEngine, trainFlag, testFlag);
    
    pCorrect(idx) = mean(prediction);
    
end
