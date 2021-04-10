%% Compute reconstruction from a LARGE field
eccX = -1.5 : 1 : 15.5;
eccY = 9.5 : -1 : -1.5;

imgEdge = 100;
imageSize = [imgEdge, imgEdge, 3];

numX = length(eccX);
numY = length(eccY);

display = displayCreate('CRT12BitDisplay');
prior = load('../sparsePrior.mat');

input = imresize(im2double(imread('./imgLarge.jpg')), 1/3);
input = input(1:numY * imgEdge, 1:numX * imgEdge, :);

output = zeros(size(input));

for idx = 1 : numX
    for idy = 1 : numY       
        inputPatch = input(cvtIdx(idy, imgEdge), cvtIdx(idx, imgEdge), :);
        imshow(inputPatch, 'InitialMagnification', 600);
        pause(1);
    end
end

% output = zeros(size(input));
% for idx = 1 : nStep
%     for idy = 1 : nStep
%         retina = ConeResponseCmosaic...
%             (eccX(idx), eccY(idy), 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);
%         
%         render = retina.forwardRender(imageSize, false, false);        
%         render = double(render);
%         
%         regPara = 2e-3; stride = 4;
%         estimator = PoissonSparseEstimator...
%             (render, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
%         
%         inputPatch = input(cvtIdx(idx, imgEdge), cvtIdx(idy, imgEdge), :);
%         
%         response = render * inputPatch(:);
%         recon = estimator.runEstimate(response, 'maxIter', 500, 'display', 'off');        
%         
%         output(cvtIdx(idx, imgEdge), cvtIdx(idy, imgEdge), :) = ...
%             gammaCorrection(recon, display);
%         
%         fprintf('x: %d, y: %d \n', idx, idy);
%     end
% end

%% Helper function
function index = cvtIdx(index, step)

startIdx = (index - 1) * step + 1;
endIdx = startIdx + step - 1;

index = startIdx : 1 : endIdx;

end