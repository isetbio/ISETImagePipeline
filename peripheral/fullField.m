%% Compute reconstruction from a LARGE field
eccX = 0.5 : 1 : 14.5;
eccY = 14.5 : -1 : 0.5;

imgEdge = 100;
imageSize = [imgEdge, imgEdge, 3];

numX = length(eccX);
numY = length(eccY);

display = displayCreate('CRT12BitDisplay');
prior = load('./sparsePrior.mat');

input = load('./largeLinear.mat');
input = input.inputLinear;

output = zeros(size(input));
for idx = 1 : numX
    renderArray = cell(1, numY);
    outputArray = cell(1, numY);
    
    for idy = 1 : numY
        retina = ConeResponseCmosaic...
            (eccX(idx), eccY(idy), 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);
        
        render = retina.forwardRender(imageSize, false, false);
        renderArray{idy} = double(render);
    end
    
    fprintf('Finished render matrix calculation \n');
    
    parfor idy = 1 : numY        
        regPara = 1.5e-3; stride = 4;
        estimator = PoissonSparseEstimator...
            (renderArray{idy}, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
        
        inputPatch = input(cvtIdx(idy, imgEdge), cvtIdx(idx, imgEdge), :);
        
        response = render * inputPatch(:);
        recon = estimator.runEstimate(response, 'maxIter', 500, 'display', 'off');
        
        outputArray{idy} = gammaCorrection(recon, display);
    end
    
    for idy = 1 : numY
        output(cvtIdx(idy, imgEdge), cvtIdx(idx, imgEdge), :) = outputArray{idy};
    end
    
    fprintf('x: %d / %d \n', idx, numX);
end

%% Helper function
function index = cvtIdx(index, step)

startIdx = (index - 1) * step + 1;
endIdx = startIdx + step - 1;

index = startIdx : 1 : endIdx;

end
