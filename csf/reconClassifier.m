function dataOut = reconClassifier(obj, operationMode, ~, nullResponses, testResponses)

% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Image Reconstruction Observer');
    return;
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end

if (strcmp(operationMode, 'train'))
    nullRecon = mean(nullResponses, 2);
    testRecon = mean(testResponses, 2);
    
    dataOut.trainedClassifier = [];
    dataOut.preProcessingConstants = struct('nullTemplate', nullRecon(:), 'testTemplate', testRecon(:));
    
    return;
end

if (strcmp(operationMode, 'predict'))
    % Get template we tucked away at training time.
    nullTemplate = obj.preProcessingConstants.nullTemplate;
    testTemplate = obj.preProcessingConstants.testTemplate;    
    
    % Make sure number of null and test instances matches.
    nTrials = size(nullResponses, 2);
    assert(nTrials == size(testResponses, 2));
        
    response = zeros(1, nTrials);
    parfor idx = 1:nTrials
        testRecon = testResponses(:, idx);
        nullRecon = nullResponses(:, idx);
        
        distCr = norm([testRecon - testTemplate; nullRecon - nullTemplate]);
        distIr = norm([testRecon - nullTemplate; nullRecon - testTemplate]);
        
        % For DV really close to 0, do a coin flip
        threshold = 1e-5;
        if (abs(distCr - distIr) <= threshold)
            response(idx) = (rand() > 0.5);
        else
            response(idx) = ((distCr - distIr) < 0);
        end
    end
    
    % Set up return
    dataOut.trialPredictions = response;
    dataOut.pCorrect = mean(response);
    
    return;
end

end


