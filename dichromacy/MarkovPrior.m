classdef MarkovPrior
    
    methods(Static)
        function mtx = covMtxMarkov(dimension, var, rho)
            
            mtx = zeros(dimension, dimension);
            for i = 1:dimension
                for j = 1:dimension
                    mtx(i, j) = rho ^ abs(i - j);
                end
            end
            mtx = mtx * var;
            
        end
        
        function [mu, covMtx, regBasis] = colorSignal(dimension, rhoS, rhoC, showPlot)
            if(~exist('showPlot', 'var'))
                showPlot = true;
            end
            
            mtx = MarkovPrior.covMtxMarkov(dimension, 0.2, rhoS);
            spatial = kron(mtx, mtx);
            
            mtx = MarkovPrior.covMtxMarkov(3, 0.2, rhoC);
            covMtx = kron(mtx, spatial);
            
            mu = 0.5 * ones(1, dimension * dimension * 3);
            
            if nargout > 2
                [U, S, ~] = svd(covMtx);
                if showPlot
                    visualizeBasis(U, dimension, dimension^2 * 3, false);
                end
                regBasis = U * diag(sqrt(diag(S)));
            end
        end
        
        function errorMtx = reconFunc(renderArray, nDim, corrSpatial, corrChromat, nRecon, showPlot)
            [mu, ~, regBasis] = MarkovPrior.colorSignal(nDim, corrSpatial, corrChromat, false);
            nMosaic = length(renderArray);
            imageSize = [nDim, nDim, 3];
            
            errorMtx = zeros(nMosaic, nRecon);
            
            parfor reconID = 1:nRecon
                sample = regBasis * normrnd(0, 1, [prod(imageSize), 1]) + mu';
                for mosaicID = 1:nMosaic
                    render = double(renderArray{mosaicID});
                    
                    estimator = RidgeGaussianEstimator(render, regBasis, mu');
                    estimator.setLambda(1e-4);
                    
                    recon = (estimator.estimate((render * sample)'))';
                    errorMtx(mosaicID, reconID) = norm(sample(:) - recon(:));
                    
                    if showPlot
                        figure();
                        subplot(1, 2, 1);
                        imshow(reshape(sample, imageSize), 'InitialMagnification', 800);
                        
                        subplot(1, 2, 2);
                        imshow(reshape(recon, imageSize), 'InitialMagnification', 800);
                    end
                end
            end
        end
    end
    
end

