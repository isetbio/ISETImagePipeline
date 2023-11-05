% List of spatial frequencies to be tested.
spatialFreqs = [2, 4, 8, 16, 25];

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'red-green';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 0.0]';
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0]';
end

% Set the RMS cone contrast of the stimulus.
% A vector length contrast of 0.09 should be OK.
rmsContrast = 0.10;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Setup a cone mosaic & generate its corresponding render matrix
fovDegs = 0.50;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'integrationTime', 0.5);
retina.PSF = ConeResponse.psfDiffLmt();

imageSize = [100, 100, 3];
render = retina.forwardRender(imageSize, false, true, false);
render = double(render);
fprintf('Compute Render Matrix... Done \n');

%% Neural engine with ConeResponse class
load('../sparsePrior.mat');
estimator = PoissonSparseEstimator(double(render), inv(regBasis), mu', 0.01, 3, imageSize);

computeFunction = @(neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin) ...
    reconNeuralEngine(estimator, retina, 'off', neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin{:});

neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralEngine = neuralResponseEngine(computeFunction, neuralParams);

%% Recon Classifier engine
classifierEngine = responseClassifierEngine(@reconClassifier, struct());
classifierPara = struct('trainFlag', 'none', 'testFlag', 'random', 'nTrain', 1, 'nTest', 32);

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.
thresholdPara = struct('logThreshLimitLow', 2.5, ...
    'logThreshLimitHigh', 0.0, ...
    'logThreshLimitDelta', 0.0125, ...
    'slopeRangeLow', 1, ...
    'slopeRangeHigh', 200, ...
    'slopeDelta', 5.0);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 640, 'maxTrial', 640, ...
    'numEstimator', 1, 'stopCriterion', 0.05);

%% Compute threshold for each spatial frequency
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
logThreshold = zeros(1, length(spatialFreqs));
questObj = {};
responseObj = {};

fprintf('Start CSF Calculation \n');
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingScene = createGratingScene(chromaDir, spatialFreqs(idx), 'fovDegs', fovDegs);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(idx), questObj{end + 1}, responseObj{end + 1}] = ...
        computeThresholdRecon(gratingScene, neuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara);
    
    if showPlot
        % Plot stimulus
        figure(dataFig);
        subplot(3, 4, idx * 2 - 1);
        
        visualizationContrast = 1.0;
        [theSceneSequence] = gratingScene.compute(visualizationContrast);
        gratingScene.visualizeStaticFrame(theSceneSequence);
        
        % Plot data and psychometric curve
        % with a marker size of 5.0
        subplot(3, 4, idx * 2);
        questObj{idx}.plotMLE(5.0);
    end       
end

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Show result
dataFig = figure();

fprintf('Start CSF Calculation \n');
for idx = 1:length(spatialFreqs)    
    gratingScene = createGratingScene(chromaDir, spatialFreqs(idx), 'fovDegs', fovDegs);
    
    % Plot stimulus
    figure(dataFig);
    subplot(3, 4, idx * 2 - 1);

    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);

    % Plot data and psychometric curve
    % with a marker size of 5.0
    subplot(3, 4, idx * 2);
    questObj{idx}.plotMLE(10.0);
end

set(dataFig, 'Position',  [0, 0, 800, 800]);
