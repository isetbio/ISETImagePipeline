% List of spatial frequencies to be tested.
spatialFreqs = [1, 2, 4, 8, 12, 15];

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0]';
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0]';
end

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  Vector length contrast of 0.09 should be fine

rmsContrast = 0.10;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Create neural response engine
%
% This calculations isomerizations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
fovDegs = 0.50;
neuralParams = nreConeResponse;
neuralParams.coneMosaicParams.fovDegs = fovDegs;
neuralParams.coneMosaicParams.timeIntegrationSeconds = 0.01;
neuralParams.coneMosaicParams.eccBased = true;
theNeuralEngine = neuralResponseEngine(@nreConeResponse, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
classifierPara = struct('trainFlag', 'none', ...
    'testFlag', 'random', ...
    'nTrain', 1, 'nTest', 256);

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.
thresholdPara = struct('logThreshLimitLow', 3.0, ...
    'logThreshLimitHigh', 0.0, ...
    'logThreshLimitDelta', 0.0125, ...
    'slopeRangeLow', 1, ...
    'slopeRangeHigh', 100, ...
    'slopeDelta', 2.5);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 2560, 'maxTrial', 2560, ...
    'numEstimator', 1, 'stopCriterion', 0.05);

%% Compute threshold for each spatial frequency

% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
dataFig = figure();
logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingScene = createGratingScene(chromaDir, spatialFreqs(idx), 'fovDegs', fovDegs);
    
    [logThreshold(idx), questObj] = ...
     computeThresholdTAFC(gratingScene, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara);
    
    % Plot stimulus
    figure(dataFig);
    subplot(3, 4, idx * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);
    
    subplot(3, 4, idx * 2);
    questObj.plotMLE(7.5);
end
% set(dataFig, 'Position',  [50, 50, 400, 1050]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
% set(theCsfFig, 'Position',  [500, 50, 600, 650]);
