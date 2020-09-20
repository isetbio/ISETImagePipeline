function [percentCorrect, stdErr] = computeResponseDiscriminability(nullStimResponseInstances, ...
    testStimResponseInstances, varargin)
?
% Compute the probability with which a binary classifier discriminates noisy responses to a test stimulus
% from noisy responses to a null stimulus. There is currently just one classifier implemented (svm) and two
% dimensionality reduction methods: 'noiseFreeResponseBased' and 'principalComponentAnalysisBased'.
% When the user selects the 'noiseFreeResponseBased' dimensionality reduction method, the supplied noise-free responses 
% to the null & test stimuli are projected into the space of these templates, thereby reducing the response dimensionality to 1. 
% Note that the dimensionality of the noise-free responses has to match that of the the noisy-response instances, so if the 
% noisy-response instances represent spatiotemporal responses, the passed noisy-free response templates must also represent 
% a spatiotemporal signal of the same dimensionality.
%
% Inputs:
%   nullStimResponseInstances           -  [Nreps x Mdim] matrix of Ntrials of noisy response instances
%                                           to the NULL stimulus, each consisting of M-dim data points
%
%   testStimResponseInstances           -  [Nreps x Mdim] matrix of Ntrials of noisy response instances
%                                           to the TEST stimulus, each consisting of M-dim data points
%
% Optional key/value pairs
%   'taskIntervals'                     - Scalar, either 1 or 2. In the one-interval task, 
%                                         null and test responses are labeled as two different classes.
%                                         In the two-interval tesk, null and test response vectors are
%                                         concatenated with the [NULL TEST] order being one class and the
%                                         [TEST-NULL] order being the other class.
%   'nullStimNoiseFreeResponse'         - [1 x Mdim] vector of the noise-free response to the NULL stimulus
%   'testStimNoiseFreeResponse'         - [1 x Mdim] vector of the noise-free response to the TEST stimulus
%   'dimensionalityReductionMethod'     - String - Valid choices:{'noiseFreeResponseBased', 'principalComponentAnalysisBased'}
%
% Outputs:
%   percentCorrect                      - Performance, percent correct the responses were discriminated correctly
%   stdErr                              - Standard error of the performance
%
% Usage:
%{
   % Noisy response instances to a test stimulus
   nTrials = 500;   % number of instances
   nRows = 16;      % spatial dimension #1
   mCols = 16;      % spatial dimension #2
   nTimeBins = 10;  % temporal dimension
   rng(1);
?
   % Generate noisy response instances to the test stimulus
   testStimResponseInstances = randn(nTrials, nRows, mCols, nTimeBins);
   testStimResponseInstances(1:nTrials, 1:nRows, mCols/2+(-1:1), 1:nTimeBins) = ...
        testStimResponseInstances(1:nTrials, 1:nRows, mCols/2+(-1:1), 1:nTimeBins) + 0.12;
?
   % Generate noisy response instances to the null stimulus
   nullStimResponseInstances = randn(nTrials, nRows, mCols, nTimeBins);
?
   % Reshape the response instances in the format expected by computeResponseDiscriminability(): nTrials x mDimensions
   mDims = nRows * mCols * nTimeBins;
   nullStimResponseInstances = reshape(nullStimResponseInstances, [nTrials mDims]);
   testStimResponseInstances = reshape(testStimResponseInstances, [nTrials mDims]);
   
   % Compute performance for an (svm-based) classifier that uses PCA with the first 5 principal components for dimensionality reduction
   percentCorrectForPCAbasedDimensionalityReduction = computeResponseDiscriminability(nullStimResponseInstances, testStimResponseInstances, ...
       'taskIntervals', 2, ...
       'PCAComponentsNum', 5, ...
       'dimensionalityReductionMethod', 'principalComponentAnalysisBased')
?
   % Compute performance for an (svm-based) classifier that uses noise-free response-based templates for dimensionality reduction
   percentCorrectForTemplateBasedDimensionalityReduction = computeResponseDiscriminability(nullStimResponseInstances, testStimResponseInstances, ...
       'taskIntervals', 2, ...
       'nullStimNoiseFreeResponse', mean(nullStimResponseInstances,1), ...
       'testStimNoiseFreeResponse', mean(testStimResponseInstances,1), ...
       'dimensionalityReductionMethod', 'noiseFreeResponseBase')
?
   % Compute performance for an (svm-based) classifier that uses no dimensionality reduction
   percentCorrectForNoDimensionalityReduction = computeResponseDiscriminability(nullStimResponseInstances, testStimResponseInstances, ...
       'taskIntervals', 2, ...
       'dimensionalityReductionMethod', 'none')
?
%}
%
%
% History
%
% 7/21/2020  NPC  Wrote it
%
?
?
    % Parse input
    p = inputParser;
    p.addRequired('nullStimResponseInstances', @isnumeric);
    p.addRequired('testStimResponseInstances', @isnumeric);
    p.addParameter('taskIntervals', 1, @(x)(isscalar(x)&&((x==1)||(x==2))));
    p.addParameter('dimensionalityReductionMethod', 'noiseFreeResponseBased', @(x)(ismember(x,{'noiseFreeResponseBased', 'principalComponentAnalysisBased', 'none'})));
    p.addParameter('nullStimNoiseFreeResponse', [], @isnumeric);
    p.addParameter('testStimNoiseFreeResponse', [], @isnumeric);
    p.addParameter('PCAComponentsNum', 100, @isscalar);
    p.addParameter('classifierType', 'svm', @(x)(ismember(x, {'svm'})));
    p.addParameter('svmClassifierParam_KernelFunction', 'linear', @(x)(ismember(x, {'linear', 'rbf'})));
    p.addParameter('crossValidationFoldsNum', 10, @isscalar);
    p.parse(nullStimResponseInstances,testStimResponseInstances,varargin{:});
    
    nullStimResponseInstances = p.Results.nullStimResponseInstances;
    testStimResponseInstances = p.Results.testStimResponseInstances;
    nullStimNoiseFreeResponse = p.Results.nullStimNoiseFreeResponse;
    testStimNoiseFreeResponse = p.Results.testStimNoiseFreeResponse;
    dimensionalityReductionMethod = p.Results.dimensionalityReductionMethod;
    PCAComponentsNum = p.Results.PCAComponentsNum;
    taskIntervals = p.Results.taskIntervals;
    classifierType = p.Results.classifierType;
    svmClassifierParam_KernelFunction = p.Results.svmClassifierParam_KernelFunction;
    crossValidationFoldsNum = p.Results.crossValidationFoldsNum;
    
    % Test that the NoiseFreeResponse inputs are either both empty, or both non empty
    assert(~xor( isempty(nullStimResponseInstances), isempty(testStimResponseInstances)), ...
        'null/test StimResponseInstances must either be both empty or both non-empty');
    
    % Test that the response dimensionalities match
    if (~isempty(nullStimNoiseFreeResponse))
        assert(all([...
            size(testStimNoiseFreeResponse,2) == size(nullStimNoiseFreeResponse,2) ...
            size(nullStimResponseInstances,2) == size(testStimResponseInstances,2) ...
            size(testStimNoiseFreeResponse,2) == size(testStimResponseInstances,2) ...
            ]), 'response dimensionalities do not match');
    else
        assert(size(nullStimResponseInstances,2) == size(testStimResponseInstances,2), ...
            'response dimensionalities do not match');
    end
    
    % Test that response repetitions numbers match
    assert(size(nullStimResponseInstances,1) == size(testStimResponseInstances,1), ...
        'response trials numbers do not match match');
    
    % Dimensionality reduction
    switch (dimensionalityReductionMethod)
        case 'noiseFreeResponseBased'
            if (~isempty(nullStimNoiseFreeResponse))
                % Project noisy response instances to their templates
                nullStimResponseInstances = nullStimResponseInstances * nullStimNoiseFreeResponse';
                testStimResponseInstances = testStimResponseInstances * testStimNoiseFreeResponse';
            else
                error('The user must pass the noise-free response to the null/test stimuli when the dimensionalityReductionMethod is set to ''noiseFreeResponseBased''.');
            end
?
        case 'principalComponentAnalysisBased'
            % Concatenate all responses
            allResponses = cat(1, nullStimResponseInstances, testStimResponseInstances);
        
            % Standardize responses
            m = mean(allResponses,1);
            s = std(allResponses,1,1);
            allResponses = (allResponses - repmat(m,size(allResponses,1),1)) ./ repmat(s,size(allResponses,1),1);
?
            % Find principal components
            principalComponents = pca(allResponses,'NumComponents',PCAComponentsNum);
?
            % Project all responses to the space defined by the selected PCs    
            allResponses = allResponses * principalComponents;
?
            % Extract the projections of the null/test stimulus response instances
            nullStimResponseInstances = allResponses(1:size(nullStimResponseInstances,1),:);
            testStimResponseInstances = allResponses(size(nullStimResponseInstances,1)+(1:size(testStimResponseInstances,1)),:);
        case 'none'
             % Do nothing. No dimensionality reduction
    end % Switch
    
    % Assemble the classification data
    nTrials = size(nullStimResponseInstances,1);
    responseSize = size(nullStimResponseInstances,2);
    
    if (taskIntervals == 1)
        classificationData = zeros(2*nTrials, responseSize);
        classLabels = zeros(2*nTrials, 1);
        
        % Class 1 data and labels
        classificationData(1:nTrials,:) = nullStimResponseInstances;
        classLabels(1:nTrials,1) = 0;
        
        % Class 2 data and labels
        classificationData(nTrials+(1:nTrials),:) = testStimResponseInstances;
        classLabels(nTrials+(1:nTrials),1) = 1;
    else
        classificationData = zeros(nTrials, 2*responseSize);
        classLabels = zeros(nTrials, 1);
        halfTrials = floor(nTrials/2);
        
        % Class 1 data and labels
        classificationData(1:halfTrials,:) = cat(2, nullStimResponseInstances(1:halfTrials,:), testStimResponseInstances(1:halfTrials,:));
        classLabels(1:halfTrials,1) = 0;
        
        % Class 2 data and labels
        classificationData(halfTrials+(1:halfTrials),:) = cat(2, testStimResponseInstances(halfTrials+(1:halfTrials),:), nullStimResponseInstances(halfTrials+(1:halfTrials),:));
        classLabels(halfTrials+(1:halfTrials),1) = 1;
    end
    
    % Classify
    switch (classifierType)
        case 'svm'
            svm = fitcsvm(classificationData, classLabels, ...
                'KernelFunction', svmClassifierParam_KernelFunction ...
                );
            CVSVM = crossval(svm,'KFold',crossValidationFoldsNum);
            percentCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
            stdErr = std(percentCorrect)/sqrt(crossValidationFoldsNum);
            percentCorrect = mean(percentCorrect);
    end % switch
end
