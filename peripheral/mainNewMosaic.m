%% Load the display setup and create mosaic object
display = displayCreate('CRT12BitDisplay');
prior = load('sparsePrior.mat');
imageSize = [100, 100, 3];

retinaEccs = [1, 0; 5, 0; 10, 0; 10, 10; 18, 0; 18, 18];
nRetina = size(retinaEccs, 1);

output = zeros([nRetina, size(inputLinear, 1), imageSize]);

%% Run reconstructions
for retinaID = 1 : nRetina
    
    % Generate cone mosaic - [eccX, eccY] deg ecc
    eccX = retinaEccs(retinaID, 1);
    eccY = retinaEccs(retinaID, 2);
    
    retina = ConeResponseCmosaic...
        (eccX, eccY, 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);
    
    render = retina.forwardRender(imageSize, false, false);    
    render = double(render);
    
    regPara = 2e-3; stride = 4;
    estimator = PoissonSparseEstimator...
        (render, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
    
    parfor idx = 1 : size(inputLinear, 1)
        input = reshape(inputLinear(idx, :, :, :), [128, 128, 3]);
        input = imresize(input, imageSize(1) / size(input, 1));
        
        response = render * input(:);
        recon = estimator.runEstimate(response, 'maxIter', 500, 'display', 'iter');        
        
        output(retinaID, idx, :, :, :) = recon;
    end
end
