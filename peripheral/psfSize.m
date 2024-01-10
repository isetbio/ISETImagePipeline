%% Load the display setup and create mosaic object
load('inputLinear.mat');

display = displayCreate('CRT12BitDisplay');
prior = load('../sparsePrior.mat');
imageSize = [100, 100, 3];

retinaEccs = [1, 0; 5, 0; 10, 0; 10, 10; 18, 0; 18, 18];
nRetina = size(retinaEccs, 1);

output = zeros([nRetina, size(inputLinear, 1), imageSize]);

ecc = 10;
retina = ConeResponseCmosaic...
    (ecc, ecc, 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);

%% Run reconstructions
for retinaID = 1 : nRetina
    fprintf('Run ID: %d \n', retinaID);
    
    % Generate cone mosaic - [eccX, eccY] deg ecc
    eccX = retinaEccs(retinaID, 1);
    eccY = retinaEccs(retinaID, 2);
    
    retina_ = ConeResponseCmosaic...
        (eccX, eccY, 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);
    
    retina.PSF = retina_.PSF;
    retina.psfData = retina_.psfData;
    
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

plotID = 1;
for rID = 1 : nRetina
    for idx = 1 : size(inputLinear, 1)
        recon = reshape(output(rID, idx, :, :, :), imageSize);
        recon = gammaCorrection(recon, display);
        
        axes(plotAxis(plotID));
        imshow(recon, 'initialMagnification', 500);
        
        plotID = plotID + 1;
    end
end

%% Error metric
inputSize = size(inputLinear);
inputAll = reshape(inputLinear, [prod(inputSize(1:3)), 3]);

loadings = pca(inputAll);

rmseImage = zeros(nRetina, size(inputLinear, 1));  
rmseSpatial = zeros(nRetina, size(inputLinear, 1));
rmseChromat = zeros(nRetina, size(inputLinear, 1));

for rID = 1 : nRetina
    for idx = 1 : size(inputLinear, 1)
        input = reshape(inputLinear(idx, :, :, :), [128, 128, 3]);
        input = imresize(input, imageSize(1) / size(input, 1));
        recon = reshape(output(rID, idx, :, :, :), imageSize);
                
        rmseImage(rID, idx) = norm(input(:) - recon(:));
        
        input = reshape(input, [100 * 100, 3]) * loadings;
        recon = reshape(recon, [100 * 100, 3]) * loadings;
        
        rmseSpatial(rID, idx) = norm(input(:, 1) - recon(:, 1));
        
        diff = input(:, [2, 3]) - recon(:, [2, 3]);
        rmseChromat(rID, idx) = norm(diff(:));
    end
end

yyaxis left
rmseSpatial = rmseSpatial - rmseSpatial(4, :);
errorbar(1:size(rmseSpatial, 1), mean(rmseSpatial, 2), ...
    std(rmseSpatial, 0, 2) / sqrt(size(rmseSpatial, 2)), '-o', 'LineWidth', 2);
ylabel('Relative RSS, Spatial');
ylim([-0.05, 3]);

yyaxis right
rmseChromat = rmseChromat - rmseChromat(6, :);
errorbar(1:size(rmseChromat, 1), mean(rmseChromat, 2), ...
    std(rmseChromat, 0, 2) / sqrt(size(rmseChromat, 2)), '-o', 'LineWidth', 2);
ylabel('Relative RSS, Chromatic'); 

box off;

xlim([0.5, 6.5]);
ylim([-0.05, 3]);
xticks(1:size(rmsePixel, 1));
xticklabels({'1, 0', '5, 0', '10, 0', '10, 10', '18, 0', '18, 18'});
xlabel('PSF at Eccentricity');

%% Plot PSF
ecc = 10;

for retinaID = 1 : nRetina
        % Generate cone mosaic - [eccX, eccY] deg ecc
    eccX = retinaEccs(retinaID, 1);
    eccY = retinaEccs(retinaID, 2);
    
    retina_ = ConeResponseCmosaic...
        (eccX, eccY, 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 9);
    
    retina.PSF = retina_.PSF;
    retina.psfData = retina_.psfData;
    
    mosaic = ConeResponsePeripheral(ecc, ecc, 'fovealDegree', 0.20);
    cmStruct = mosaic.Mosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    coneData.conePositionsArcMin = cmStruct.coneLocs * 60;
    coneData.coneAperturesArcMin = cmStruct.coneApertures * 60;
    
    [~, wIdx] = min(abs(retina.psfData.supportWavelength-550));
    wavePSF = squeeze(retina.psfData.data(:,:,wIdx));
    zLevels = 0.1:0.2:0.9;
    xyRangeArcMin = 8 * [-1, 1];
    
    figure();
    PolansOptics.renderPSF(gca(), ...
        retina.psfData.supportX, retina.psfData.supportY, wavePSF/max(wavePSF(:)), ...
        xyRangeArcMin, zLevels,  gray(1024), [0.5 0.5 0.5], 'withConeData', coneData);
    
    pause(1); grid off;
    xlim(xlim() * 0.60); ylim(ylim() * 0.60);
    xlabel(''); ylabel('');
    xticks(''); yticks('');
end
