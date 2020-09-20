%% setup variables and constants
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'integrationTime', 0.05);
retina.visualizeMosaic();

% estimate on the log contrast domain
estDomain  = -5 : 0.025 : 0;
slopeRange = 0.1 : 0.5 : 100;

observer = PoissonTemplateObserver(retina, display.CRT12BitDisplay, 'L+M+S', 1);

%% Single QUEST+ object with fixed number of trials
estimator = QuestThresholdEstimator('minTrial', 128, 'maxTrial', 256, ...
                                    'estDomain', estDomain, 'numEstimator', 1, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    response = observer.singleTrial(stimCrst) + 1;
    
    [crst, flag] = estimator.singleTrial(crst, response);
    [threshold, ~] = estimator.thresholdEstimate();
    
    fprintf('Trial count: %d, threshold estimate: %.3f \n', estimator.nTrial, threshold);
end

% Show results
[threshold, stderr] = estimator.thresholdEstimate();
fprintf('%d trials recorded, (log) threshold estimate: %.2f +/- %.2f \n', estimator.nTrial, threshold, stderr);

[threshold, para] = estimator.thresholdMLE('showPlot', true);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Single QUEST+ object but when it's easy to run multiple trials per contrast
estimator = QuestThresholdEstimator('minTrial', 1024, 'maxTrial', 1024, ...
                                    'estDomain', estDomain, 'numEstimator', 1, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();

% 32 trial for each contrast level
nRepeat = 32; 
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    [~, response] = observer.multiTrial(stimCrst, nRepeat);
    response = response + 1;    
    
    [crst, flag] = estimator.multiTrial(crst * ones(1, nRepeat), response);
    [threshold, ~] = estimator.thresholdEstimate();
    
    fprintf('Trial count: %d, threshold estimate: %.3f \n', estimator.nTrial, threshold);
end

% Show results
[threshold, stderr] = estimator.thresholdEstimate();
fprintf('%d trials recorded, (log) threshold estimate: %.2f +/- %.2f \n', estimator.nTrial, threshold, stderr);

[threshold, para] = estimator.thresholdMLE('showPlot', true);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Multiple QUEST+ object for adaptive procedure
