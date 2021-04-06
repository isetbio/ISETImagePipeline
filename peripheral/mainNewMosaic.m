%% Load the display setup and create mosaic object
load('inputLinear.mat');

display = displayCreate('CRT12BitDisplay');
prior = load('../sparsePrior.mat');
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

%% Show results
plotAxis = tight_subplot(nRetina, size(inputLinear, 1), [.01 .01], [.01 .01], [.01 .01]);

rmsePixel = zeros(nRetina, size(inputLinear, 1));
rmseGray = zeros(nRetina, size(inputLinear, 1));

plotID = 1;
for rID = 1 : nRetina
    for idx = 1 : size(inputLinear, 1)
        input = reshape(inputLinear(idx, :, :, :), [128, 128, 3]);
        input = imresize(input, imageSize(1) / size(input, 1));
        input = gammaCorrection(input, display);
        
        recon = reshape(output(rID, idx, :, :, :), imageSize);
        recon = gammaCorrection(recon, display);
        
        inputGray = rgb2gray(input); reconGray = rgb2gray(recon);
        rmsePixel(rID, idx) = norm(input(:) - recon(:));
        rmseGray(rID, idx) = norm(inputGray(:) - reconGray(:));
        
        axes(plotAxis(plotID));
        imshow(recon, 'initialMagnification', 500);
        
        plotID = plotID + 1;
    end
end

% Plot RMSE
figure();
yyaxis left
errorbar(1:size(rmsePixel, 1), mean(rmsePixel, 2), ...
    std(rmsePixel, 0, 2)/sqrt(size(rmsePixel, 2)), '-o', 'LineWidth', 2);
ylabel('RMSE RGB');

yyaxis right
errorbar(1:size(rmseGray, 1), mean(rmseGray, 2), ...
    std(rmseGray, 0, 2)/sqrt(size(rmseGray, 2)), '-o', 'LineWidth', 2);
ylabel('RMSE Gray');

box off;

xlim([0.5, 6.5]);
xticks(1:size(rmsePixel, 1));
xticklabels({'1, 0', '5, 0', '10, 0', '10, 10', '18, 0', '18, 18'});
xlabel('Visual Eccentricity');
