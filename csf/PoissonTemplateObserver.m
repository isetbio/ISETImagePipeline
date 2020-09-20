classdef PoissonTemplateObserver < Observer
    
    properties
        retina;
        display;
        
        bgResponse;
    end
    
    methods (Access = public)
        
        function this = PoissonTemplateObserver(retina, display, stimType, stimFreq)
            this@Observer(stimType, stimFreq);
            
            this.retina = retina;
            this.display = display;
            
            % use 0-contrast stim as background, compute its response
            [~, ~, ~, rateB] = this.retina.compute(this.getStimulus(0.0));
            this.bgResponse = rateB;
            
        end
                
        % simulate a single trial response
        % @Override
        function response = singleTrial(this, stimCrst)
            targetStim = this.getStimulus(stimCrst);
            [~, ~, ~, tgResponse] = this.retina.compute(targetStim);
            sampleT = poissrnd(tgResponse);
            
            llB = sum(log(poisspdf(sampleT, this.bgResponse)));
            llT = sum(log(poisspdf(sampleT, tgResponse)));
            response = (llT - llB) > 0;
        end
        
        % simulate mutiple trials
        function [pCorr, choice] = multiTrial(this, stimCrst, nTrial)
                                                  
            targetStim = this.getStimulus(stimCrst);
            [~, ~, ~, tgResponse] = this.retina.compute(targetStim);                        
            
            sampleT = poissrnd(repmat(tgResponse, 1, nTrial));
            llB = sum(log(poisspdf(sampleT, repmat(this.bgResponse, 1, nTrial))), 1);
            llT = sum(log(poisspdf(sampleT, repmat(tgResponse, 1, nTrial))), 1);
            logRatio = llT - llB;
            
            % choice generation
            choice = (logRatio > 0);
            pCorr = (sum(choice) / nTrial);
        end
        
    end
end

