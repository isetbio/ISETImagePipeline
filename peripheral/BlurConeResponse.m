classdef BlurConeResponse < ConeResponsePeripheral
    
    properties (Access = public)
        blurSD;
    end
    
    methods(Access = public)
        function this = BlurConeResponse(eccX, eccY, blurSD, varargin)
            this@ConeResponsePeripheral(eccX, eccY, varargin{:});
            this.blurSD = blurSD;
        end
        
        % Override
        function [excitation, theOI, linearizedImage, allCone, L, M, S] = compute(this, image)
            meanLuminanceCdPerM2 = [];
            [~, ~, linearizedImage] = sceneFromFile(image, 'rgb', ...
                meanLuminanceCdPerM2, this.Display);
            blurInput = imgaussfilt(linearizedImage, this.blurSD, 'Padding', 0, 'FilterDomain', 'spatial');            
            zeroIdx = (blurInput == 0);
            
            blurInput = invGammaCorrection(blurInput, this.Display);
            blurInput(zeroIdx) = 0;
            [excitation, theOI, ~, allCone, L, M, S] = compute@ConeResponse(this, blurInput);
        end
        
    end
end

