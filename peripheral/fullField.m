%% Compute reconstruction from a LARGE field
eccX = -2.5 : 1 : 14.5;
eccY = 14.5 : -1 : -2.5;

nStep = length(eccX);
assert(nStep == length(eccY));

imgEdge = 64;
imageSize = [64, 64, 3];

display = displayCreate('CRT12BitDisplay');
prior = load('../sparsePrior.mat');

input = im2double(imread('./imgLarge.jpeg'));
input = input(:, floor(size(input, 2) * 0.2) : ceil(size(input, 2) * 0.8), :);
input = sampleImage(input, nStep * imageSize(1)); 

output = zeros(size(input));
for idx = 1 : nStep
    for idy = 1 : nStep
        retina = ConeResponseCmosaic...
            (eccX(idx), eccY(idy), 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);
        
        render = retina.forwardRender(imageSize, false, false);        
        render = double(render);
        
        regPara = 2e-3; stride = 4;
        estimator = PoissonSparseEstimator...
            (render, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
        
        inputPatch = input(cvtIdx(idx, imgEdge), cvtIdx(idy, imgEdge), :);
        
        response = render * inputPatch(:);
        recon = estimator.runEstimate(response, 'maxIter', 500, 'display', 'off');        
        
        output(cvtIdx(idx, imgEdge), cvtIdx(idy, imgEdge), :) = ...
            gammaCorrection(recon, display);
        
        fprintf('x: %d, y: %d \n', idx, idy);
    end
end

%% Helper function
function index = cvtIdx(index, step)

startIdx = (index - 1) * step + 1;
endIdx = startIdx + step - 1;

index = startIdx : 1 : endIdx;

end