classdef QuestThresholdEstimator < ContrastThresholdEstimator
    %QuestThresholdEstimator  Adaptive contrast threshold estimation procedure
    % based on the QUEST+ routine. The class maintain an array of questData
    % ojbect and loop through them sequentially. A stop criterion can be
    % triggered if the standard error among the set of questData object
    % drops below a pre-defined a threshold and we have the minimal number of
    % trials, or maximum number of trials is reached. 
    %
    % Usage:
    %   See t_questAdaptiveEstimator.m
    %   Also see base class QuestThresholdEstimator
    %
    %
    % ContrastThresholdEstimator Properties:
    %   estimators     - The array of questData object.
    %   numEstimator   - The number questData object to maintain
    %   stopCriterion  - Stop criterion for deciding if we have enough trials
    %   slopeRange     - Range of slopes for psychometric curve
    %   guessRate      - Range of guess rate for psychometric curve
    %   lapseRate      - Range of lapseRate rate for psychometric curve
    %   estIdx         - Current questData object being used
    %   
    %   Also see base class QuestThresholdEstimator
    %
    % QuestThresholdEstimator Methods:
    %   thresholdEstimate    - Current running estimate of threshold and
    %                          its standard error across questData object
    %   combineData          - Return all stimulus - response data we have
    %                          recorded so far
    %
    %   thresholdMLE         - Run a MLE on all data we have recorded so
    %                          far, plot data and the psychometric curve
    %
    %   Also see base class QuestThresholdEstimator 
    %
    % Inputs:
    %   None.
    %
    % Outputs:
    %   QuestThresholdEstimator Object.
    %
    % Optional key/value pairs:
    %   'numEstimator'    - Int. Number of questData object to run
    %
    %   'stopCriterion'   - Double. Should be between [0, 1]. Stop the
    %                       procedure if the standard error drops below 
    %                       stopCriterion * thresholdEstimate
    %
    %   'slopeRange'       - Array. An array of all possible slope for the
    %                        psychometric curve
    %    
    %   'guessRate'        - Array. An array of all possible guess rate for
    %                        the psychometric curve
    %  
    %   'lapseRate'        - Array. An array of all possible lapse rate for
    %                        the psychometric curve
    % 
    %    Also see base class QuestThresholdEstimator 
    
    
    % Class properties
    properties %(Access = private)
        
        estimators;
        numEstimator;
        stopCriterion;
        
        slopeRange;
        guessRate;
        lapseRate;
        
        estIdx;
        
    end
    
    methods
        
        % Constructor method
        function this = QuestThresholdEstimator(varargin)
            
            this@ContrastThresholdEstimator(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('numEstimator', 1, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('stopCriterion', 0.05, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('slopeRange', 0.1 : 0.5 : 50);
            p.addParameter('guessRate', 0.5);
            p.addParameter('lapseRate', 0.0);
            
            
            parse(p, varargin{:});
            this.numEstimator  = p.Results.numEstimator;
            this.stopCriterion = p.Results.stopCriterion;
            this.slopeRange = p.Results.slopeRange;
            this.guessRate = p.Results.guessRate;
            this.lapseRate = p.Results.lapseRate;
            
            % Initialize QUEST+ objects specified by 'numEstimator'
            this.estimators = cell(this.numEstimator, 1);
            for idx = 1 : this.numEstimator
                this.estimators{idx} = ...
                    qpInitialize('stimParamsDomainList', {this.estDomain}, ...
                    'psiParamsDomainList',  {this.estDomain, this.slopeRange, this.guessRate, this.lapseRate});
            end
            
            % Set the current estimator to #1
            this.estIdx = 1;
            
            this.nTrial = 0;
            this.nextFlag = true;
            this.testCrst = qpQuery(this.estimators{this.estIdx});
            
        end
        
        % Implement single trial and multi trial protocol
        function [nextCrst, nextFlag] = singleTrial(this, stim, response)
            
            this.estimators{this.estIdx} = qpUpdate(this.estimators{this.estIdx}, stim, response);
            
            this.nTrial = this.nTrial + 1;
            this.estIdx = mod(this.estIdx, this.numEstimator) + 1;            
            
            if this.nTrial >= this.maxTrial
                this.nextFlag = false;
            end
            
            [threshold, stderr] = thresholdEstimate(this);
            % Stop only if  stderr drop below criterion and we have at least minTrial # of trials
            if ((stderr / abs(threshold)) < this.stopCriterion) && this.nTrial >= this.minTrial
                this.nextFlag = false;
            end
            
            this.testCrst = qpQuery(this.estimators{this.estIdx});
            
            [nextCrst, nextFlag] = this.nextStimulus();
            
        end
        
        function [nextCrst, nextFlag] = multiTrial(this, stimVec, responseVec)
            
            for idx = 1:length(stimVec)
                this.singleTrial(stimVec(idx), responseVec(idx));
            end
            
            [nextCrst, nextFlag] = this.nextStimulus();
            
        end
        
        function estimates = parameterEstimate(this)
            
            % For now, hard code the number of parameter to be 4
            estimates = zeros(this.numEstimator, 4);
            
            for idx = 1 : this.numEstimator
                estimator = this.estimators{idx};
                psiParamsIndex = qpListMaxArg(estimator.posterior);
                estimates(idx, :) = estimator.psiParamsDomain(psiParamsIndex, :);
            end
            
        end
        
        function [threshold, stderr] = thresholdEstimate(this)
            
            estimates = this.parameterEstimate();
            
            % Assume threshold parameter is the 1st one
            threshold = mean(estimates(:, 1));
            stderr = std(estimates(:, 1)) / sqrt(this.numEstimator);
            
        end
        
        function [stimVec, responseVec, structVec] = combineData(this)
            
            structVec = [];
            for idx = 1 : this.numEstimator
                questData = this.estimators{idx};
                structVec = [structVec; questData.trialData];
            end
            
            nTotal = length(structVec);
            
            stimVec = zeros(1, nTotal);
            responseVec = zeros(1, nTotal);
            
            for idx = 1:nTotal
                stimVec(idx) = structVec(idx).stim;
                responseVec(idx) = structVec(idx).outcome - 1;
            end
            
        end
        
        % Combine all data to make a final MLE
        function [threshold, para] = thresholdMLE(this, varargin)
            
            [stimVec, responseVec, structVec] = this.combineData();
            
            questData = this.estimators{1};
            psiParamsIndex = qpListMaxArg(questData.posterior);
            psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
            
            para = qpFit(structVec, questData.qpPF, psiParamsQuest, questData.nOutcomes, ...
                'lowerBounds', [min(this.estDomain) min(this.slopeRange) min(this.guessRate) min(this.lapseRate)], ...
                'upperBounds', [max(this.estDomain) max(this.slopeRange) max(this.guessRate) max(this.lapseRate)]);
            
            threshold = para(1);
            
            p = inputParser;
            p.addParameter('showPlot', false);
            parse(p, varargin{:});
            
            if p.Results.showPlot
                figure();
                stimVal = unique(stimVec);
                for idx = 1:length(stimVal)
                    prop = responseVec(stimVec == stimVal(idx));
                    
                    scatter(stimVal(idx), sum(prop) / length(prop), 2.5e3 / length(stimVec) * length(prop), ...
                        'MarkerEdgeColor', zeros(1, 3), 'MarkerFaceColor', ones(1, 3) * 0.5, 'MarkerFaceAlpha', 0.5);
                    hold on;
                end
                                
                stimSpace = this.estDomain;
                
                fitCurve = qpPFWeibull(stimSpace', para);
                plot(stimSpace, fitCurve(:, 2), '-', 'Color', [1.0 0.2 0.0], 'LineWidth', 3);
            end
            
        end
        
    end
end

