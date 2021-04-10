%% Compute reconstruction from a LARGE field
% Setup constants
tbUseProject('ISETImagePipeline');
cd ./peripheral;

cond = 2;
switch cond
    case 1
        eccX = 0.5 : 1 : 14.5;
        eccY = 14.5 : -1 : 0.5;
        fovDegs = 1.0;
    case 2
        eccX = 5 : 2 : 15;
        eccY = 15 : -2 : 5;
        fovDegs = 2.0;
    otherwise
        error('Invalid Condition');
end

imgEdge = 100;
imageSize = [imgEdge, imgEdge, 3];

numX = length(eccX);
numY = length(eccY);

display = displayCreate('CRT12BitDisplay');
prior = load('./sparsePrior.mat');

input = load('./largeLinear0.mat');
input = input.input;
input = imresize(input, imgEdge * numX / size(input, 1));

%% Run computation
output = zeros(size(input));
for idx = 1 : numX
    renderArray = cell(1, numY);
    outputArray = cell(1, numY);
    
    startX = (idx - 1) * imgEdge + 1;
    endX = startX + imgEdge - 1;
    
    fileName = strcat('./renderArray/', sprintf('rnd_x%d', idx), '.mat');
    
    if isfile(fileName)
        savedRender = load(fileName);
        renderArray = savedRender.renderArray;
    else
        for idy = 1 : numY
            retina = ConeResponseCmosaic...
                (eccX(idx), eccY(idy), 'fovealDegree', fovDegs, 'pupilSize', 3.0, 'subjectID', 9);
            
            render = retina.forwardRender(imageSize, false, false);
            renderArray{idy} = double(render);
        end
        save(fileName, 'renderArray', '-v7.3');
    end
    
    fprintf('Finished render matrix calculation \n');
    
    for idy = 1 : numY
        regPara = 1.5e-3; stride = 4; render = renderArray{idy};
        estimator = PoissonSparseEstimator...
            (render, inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
        
        startY = (idy - 1) * imgEdge + 1;
        endY = startY + imgEdge - 1;
        
        inputPatch = input(startY:endY, startX:endX, :);
        
        response = render * inputPatch(:);
        recon = estimator.runEstimate(response, 'maxIter', 500, 'display', 'off');
        
        outputArray{idy} = gammaCorrection(recon, display);
    end
    
    for idy = 1 : numY
        startY = (idy - 1) * imgEdge + 1;
        endY = startY + imgEdge - 1;
        
        output(startY:endY, startX:endX, :) = outputArray{idy};
    end
    fprintf('Current x: %d / %d \n', idx, numX);
end

save result.mat input output;
