classdef ReconTemplateObserver < Observer
    
    properties
        retina;
        display;
        stimSize;
        render;
        estimator;
        
        bgRecon;
    end
    
    methods (Access = public)
        
        function this = ReconTemplateObserver(retina, display, stimType, stimFreq, stimSize, reconPara)
            this@Observer(stimType, stimFreq);
            
            this.retina = retina;            
            this.display = display;
            this.stimSize = stimSize;
            
            this.render = this.retina.forwardRender(this.stimSize);
            this.estimator = PoissonSparseEstimator(this.render, reconPara.basis, reconPara.mu, reconPara.reg, reconPara.stride, this.stimSize);
            
            % use 0-contrast stim as background
            bgStim = this.getStimulus(0.0);
            this.bgRecon = this.estimator(this.render * bgStim(:), 1e3, rand([prod(this.stimSize), 1]), true);
            
        end
        
        % simulate a single trial response
        % @Override
        function response = singleTrial(this, stimCrst)
            tgStim = this.getStimulus(stimCrst);
            tgResp = this.render * tgStim(:);
            tgRecon = this.estimator(tgResp, 1e3, rand([prod(this.stimSize), 1]), true);
            
            sampleRecon = this.estimator(poissrnd(tgResp), 1e3, rand([prod(this.stimSize), 1]), true);
            
            bgDist = norm(this.bgRecon(:) - sampleRecon(:));
            tgDist = norm(tgRecon(:) - sampleRecon(:));
            response = (tgDist < bgDist);
        end
        
        
    end
end

