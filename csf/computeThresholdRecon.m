function [threshold, questObj, allResponse] = computeThresholdRecon...
    (theSceneEngine, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara)

estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;

estimator = questThresholdEngine(...
    'minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, ...
    'numEstimator', questEnginePara.numEstimator, ...
    'stopCriterion', questEnginePara.stopCriterion);

% Generate the NULL stimulus (zero contrast)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

allResponse = {};
[logContrast, nextFlag] = estimator.nextStimulus();
while (nextFlag)
    testContrast = 10 ^ logContrast;
    [thTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    oneSided = true;
    [predictions, ~, responses] = computePerformanceRecon(...
        theNullSceneSequence, thTestSceneSequence, theSceneTemporalSupportSeconds, ...
        classifierPara.nTrain, classifierPara.nTest, theNeuralEngine, classifierEngine, ...
        classifierPara.trainFlag, classifierPara.testFlag, oneSided);
    
    responses.logContrast = logContrast;
    allResponse{end + 1} = responses;
    
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, classifierPara.nTest), predictions);
end

% Return threshold value
[threshold, para] = estimator.thresholdMLE('showPlot', false);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

% Return the quest+ object wrapper for plotting and/or access to data
questObj = estimator;

end
