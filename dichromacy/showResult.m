%% Load array of retina and input image, and reconstruction
nRecon = size(outputArray_M1{1}, 1);
nRetina = length(mosaicArray);
inputImage = inputLinear;
outputArray = outputArray_M1;

%% Show results
plotAxis = tight_subplot(nRetina + 1, nRecon, [.01 .03], [.1 .01], [.01 .01]);
for idx = 1:nRecon
    axes(plotAxis(idx));
    rgbImage = reshape(inputImage(idx, :, :, :), imageSize);
    imshow(invGammaCorrection(rgbImage, display.CRT12BitDisplay), 'InitialMagnification', 400);
end

for idr = 1:nRetina
    output = outputArray{idr};
    baseIdx = idr * nRecon;
    
    for idi = 1:nRecon
        rgbImage = reshape(output(idi, :, :, :), imageSize);
        
        plotIdx = baseIdx + idi;
        axes(plotAxis(plotIdx));
        imshow(invGammaCorrection(rgbImage, display.CRT12BitDisplay), 'InitialMagnification', 400);
    end
end

%% Compare two visualization method
tbUseProject('ISETImagePipeline');
plotAxis = tight_subplot(4, nRecon , [.01 .03], [.1 .01], [.01 .01]);
for idx = 1:nRecon
    axes(plotAxis(idx));
    rgbImage = reshape(inputImage(idx, :, :, :), imageSize);
    imshow(invGammaCorrection(rgbImage, display.CRT12BitDisplay), 'InitialMagnification', 400);
end

dichrArray = 8;
for idr = 1
    output = outputArray{dichrArray(idr)};
    baseIdx = idr * nRecon;
    
    for idi = 1:nRecon
        rgbImage = reshape(output(idi, :, :, :), imageSize);
        
        plotIdx = baseIdx + idi;
        axes(plotAxis(plotIdx));
        imshow(invGammaCorrection(rgbImage, display.CRT12BitDisplay), 'InitialMagnification', 400);
    end
end

% Brettel
baseIdx = nRecon * 2;
for idx = 1:nRecon
    axes(plotAxis(baseIdx + idx));
    rgbImage = reshape(inputImage(idx, :, :, :), imageSize);
    imshow(lmsVis(rgbImage, 2, 'brettel'), 'InitialMagnification', 400);
end

% Linear Transformation
baseIdx = nRecon * 3;
for idx = 1:nRecon
    axes(plotAxis(baseIdx + idx));
    rgbImage = reshape(inputImage(idx, :, :, :), imageSize);
    imshow(lmsVis(rgbImage, 2, 'linear'), 'InitialMagnification', 400);
end

%% Difference between original image and reconstruction with different preceptual metric
tbUse({'Psychtoolbox-3','isetcam', 'ISETPipelineToolbox'});

errorMSE = zeros(nRecon, nRetina);
errorLAB = zeros(nRecon, nRetina);
errorSLAB = zeros(nRecon, nRetina);

for idx = 1:nRetina
    [errorMSE(:, idx), errorLAB(:, idx), errorSLAB(:, idx)] = reconError(inputImage, outputArray{idx}, imageSize, nRecon, display);
end

%% Show plots
plotError(ratio, errorMSE, nRecon, '-ok', 'MSE Distance');
plotError(ratio, errorLAB, nRecon, '-or', 'LAB Distance');
plotError(ratio, errorSLAB, nRecon, '-ob', 'SLAB Distance');

%% Show plot: S cone ratio
figure(); subplot(1, 2, 1);
errorbar(ratio, mean(errorMSE), std(errorMSE) / sqrt(nRecon), '-ok', 'LineWidth', 1.5);
xlabel('Ratio of L cone');
ylabel('RMSE');

subplot(1, 2, 2);
errorbar(ratio, mean(errorLAB), std(errorLAB) / sqrt(nRecon), '-or', 'LineWidth', 1.5);
xlabel('Ratio of L cone');
ylabel('LAB');

%% Helper function: preceptual error metric
function [errorMSE, errorLAB, errorSLAB] = reconError(original, recon, imageSize, nRecon, display)

errorMSE  = zeros(nRecon, 1);
errorLAB  = zeros(nRecon, 1);
errorSLAB = zeros(nRecon, 1);

for idx = 1:nRecon
    
    rgbOriginal = invGammaCorrection(reshape(original(idx, :, :, :), imageSize), display.CRT12BitDisplay);
    rgbReconst  = invGammaCorrection(reshape(recon(idx, :, :, :), imageSize), display.CRT12BitDisplay);
    
    errorMSE(idx) = norm(rgbOriginal(:) - rgbReconst(:));
    [errorLAB(idx), errorSLAB(idx)] = labDistance(rgbOriginal, rgbReconst);
    
end

end

%% Plot error metric as a function of L cone ratio
function plotError(ratio, error, nRecon, formatStr, ylabelText)

figure(); subplot(1, 2, 1);
errorbar(ratio, mean(error), std(error) / sqrt(nRecon), formatStr, 'LineWidth', 1);
set(gca, 'XDir','reverse');
xlabel('M cone ratio');
ylabel(ylabelText);

set(gca, 'box', 'off')
set(gca, 'TickDir', 'out');

subplot(1, 2, 2);
plot(log(ratio), log(mean(error)), formatStr, 'LineWidth', 1); hold on;

infPoint = 1e-3;
plot(log(infPoint), log(mean(error(:, end))), formatStr, 'LineWidth', 1);

xlabel(strcat('log M cone ratio'));
ylabel(strcat('log', ylabelText));

xticks(log(flip(ratio)));
xticklabels(arrayfun(@(x) num2str(x), flip(ratio), 'UniformOutput', false));

error = mean(error);
yticks(log(error));
yticklabels(arrayfun(@(x) num2str(x), error, 'UniformOutput', false));

set(gca, 'box', 'off')
set(gca, 'TickDir', 'out');
set(gca, 'XDir','reverse');
set(gcf,'Position',[0, 0, 900, 300])
end

