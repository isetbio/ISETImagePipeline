function dataOut = reconClassifier(reconObj, obj, operationMode, ~, nullResponses, testResponses)

% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Image Reconstruction Observer');
    return;
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end

imageLength = prod(reconObj.Size);
if (strcmp(operationMode, 'train'))
    nullConeVec = mean(nullResponses, 2);
    testConeVec = mean(testResponses, 2);
    
    imageSize = reconObj.Size;
    mask = makeMast(imageSize(1));
    
    templateInit = rand([imageLength, 1]);
    nIteration = 500;
    nullRecon = (reconObj.estimate(nullConeVec, nIteration, templateInit, true, 1.0, 'off')) .* mask;
    testRecon = (reconObj.estimate(testConeVec, nIteration, templateInit, true, 1.0, 'off')) .* mask;
    
    dataOut.trainedClassifier = [];
    dataOut.preProcessingConstants = struct('nullTemplate', nullRecon(:), 'testTemplate', testRecon(:), 'mask', mask);
    
    return;
end

if (strcmp(operationMode, 'predict'))
    % Get template we tucked away at training time.
    nullTemplate = obj.preProcessingConstants.nullTemplate;
    testTemplate = obj.preProcessingConstants.testTemplate;
    mask = obj.preProcessingConstants.mask;
    
    % Make sure number of null and test instances matches.
    nTrials = size(nullResponses, 2);
    assert(nTrials == size(testResponses, 2));
    
    nIteration = 200;
    response = zeros(1, nTrials);
    parfor idx = 1:nTrials
        testRecon = (reconObj.estimate(testResponses(:, idx), nIteration, rand([imageLength, 1]), true, 1.0, 'off')) .* mask;
        
        distCr = norm(testRecon(:) - testTemplate);
        distIr = norm(testRecon(:) - nullTemplate);
        
        % For DV extremely close to 0, do a coin flip
        threshold = 1e-6;
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

function maskImage = makeMast(sideLength)

maskImage = ones([sideLength, sideLength, 3]);

center = sideLength * 0.5;
radius = sideLength * 0.4;

for i = 1:sideLength
    for j = 1:sideLength
        dist = sqrt((i - center) ^ 2 + (j - center) ^ 2);
        if (dist > radius)
            maskImage(i, j, :) = 0;
        end
    end
end

end
