classdef ContrastThresholdEstimator < handle
    %ContrastThresholdEstimator  Base class for adaptive contrast threshold estimation procedure.
    %
    % Usage:
    %   See instance class QuestThresholdEstimator
    %
    %
    % ContrastThresholdEstimator Properties:
    %   minTrial       - The minimal number of trials to run
    %   maxTrial       - The maximum number of trials to run
    %   estDomain      - The domain for which contrast is being estimated
    %   nTrial         - Number of trials ran so far
    %   testCrst       - Next test contrast value
    %   nextFlag       - Indicate if new trials are still needed
    %
    % ContrastThresholdEstimator Methods:
    %   nextStimulus         - What contrast should be the next trial and
    %                          if new trials are still required
    %
    %   singleTrial          - Submit the result of a single trial, and return
    %                          'nextStimulus'
    %
    %   multiTrial           - Submit the result of multiple trials, and return
    %                          'nextStimulus'
    %
    % Inputs:
    %   None.
    %
    % Outputs:
    %   ContrastThresholdEstimator Object.
    %
    % Optional key/value pairs:
    %   'minTrial'        - Int. The minimal number of trials to run
    %
    %   'maxTrial'        - Int. The maximum number of trials to run
    %
    %   'estDomain'       - Array. The domain (range) for which contrast is
    %                       being estimated. For examples, linear contrast
    %                       [0:0.01:1] or log contrast [-10:0.1:0]. User
    %                       should handle the conversion externally
    
    
    properties (Access = public)
        
        minTrial;
        maxTrial;
        
        estDomain;
        
        nTrial;
        testCrst;
        nextFlag;
        
    end              
    
    methods (Access = public)
        
        function this = ContrastThresholdEstimator(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('minTrial', 64, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('maxTrial', 256, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('estDomain', 0 : 0.01 : 1);
            
            parse(p, varargin{:});
            this.minTrial = p.Results.minTrial;
            this.maxTrial = p.Results.maxTrial;
            this.estDomain = p.Results.estDomain;
            
        end
        
        function [crst, flag] = nextStimulus(this)
            
            crst = this.testCrst;
            flag = this.nextFlag;
            
        end
        
    end
    
    methods (Abstract)
        
        [nextCrst, nextFlag] = singleTrial(this, stim, response)
        
        [nextCrst, nextFlag] = multiTrial(this, stimVec, responseVec)
        
    end
    
end

