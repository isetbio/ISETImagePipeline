%% Setup a cone mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

fovDegs = 0.20;
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

retina.resetSCone();
retina.reassignSCone(0.1);
retina.visualizeMosaic();

%% Regular optics
retina.PSF = oiCreate('human', pupilSize);

render = retina.forwardRender(imageSize, true, true, false);
render = double(render);

%% Run image reconstruction
% construct image estimator
regPara = 1e-3; stride = 4;
useGPU = true;
try
    estimator = PoissonSparseEstimator(gpuArray(single(render)), ...
        inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
catch NoGPU
    estimator = PoissonSparseEstimator(render, ...
        inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
    useGPU = false;
end

figure();
inputSize = [1024, 1024, 3];
outputOptics = zeros([6, imageSize]);
for idx = 1:size(inputImage, 1)
    input = reshape(inputImage(idx, :, :, :), inputSize);
    
    subplot(2, 6, idx);
    imshow(input, 'InitialMagnification', 500);
    
    [~, ~, ~, allCone] = retina.compute(input);
    reconImage = estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    outputOptics(idx, :, :, :) = reconImage;
    
    subplot(2, 6, 6 + idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);    
end

%% Turn off optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

%% Run image reconstruction
figure();
inputSize = [1024, 1024, 3];
outputDiflmt = zeros([6, imageSize]);
for idx = 1:size(inputImage, 1)
    input = reshape(inputImage(idx, :, :, :), inputSize);
    
    subplot(2, 6, idx);
    imshow(input, 'InitialMagnification', 500);
    
    [~, ~, ~, allCone] = retina.compute(input);
    reconImage = estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    outputDiflmt(idx, :, :, :) = reconImage;
    
    subplot(2, 6, 6 + idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);    
end
