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

%% construct image estimator
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

%% Run image reconstruction
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

%% Generate special cone mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

fovDegs = 1.0;
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 0.1);

retina.resetSCone();
retina.reassignSCone(0.1);
retina.visualizeMosaic();

eccDeg = 10.0;
[~, psf] = PeripheralModel.eyeModelEcc(eccDeg, eccDeg, fovDegs, pupilSize);

retina.PSF = psf;

%% Generate render matrix
render = retina.forwardRender(imageSize, true, true, false);
render = double(render);

%% Test with regular optics
spatialFreq = 6;

rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 1.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

gratingScene = ...
    createGratingScene(chromaDir, spatialFreq, 'fovDegs', fovDegs, 'pixelsNum', 1024);

crst = 0.50;
[theSceneSequence, ~] = gratingScene.compute(crst);

% Compute response from scene
[~, allCone] = retina.computeWithScene(theSceneSequence{:});
visualizeOpticalImage(retina.LastOI, 'displayRadianceMaps', true, ...
    'displayRetinalContrastProfiles', true);

retina.visualizeExcitation();

%% construct image estimator
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

%% Reconstruction
% estimate images
reconImage = ...
    estimator.runEstimate(allCone, 'maxIter', 150, ...
    'display', 'iter', 'gpu', useGPU);

figure();
imshow(gammaCorrection(reconImage, display), ...
    'InitialMagnification', 500);

%% Run reconstruction
stimFreq = [2, 8, 16, 32, 64, 128];
rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 1.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

outputOptics = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    gratingScene = ...
        createGratingScene(chromaDir, stimFreq(idx), 'fovDegs', fovDegs, 'pixelsNum', 2048);
    
    crst = 1.0;
    [theSceneSequence, ~] = gratingScene.compute(crst);
    
    % Compute response from scene
    [~, allCone] = retina.computeWithScene(theSceneSequence{:});
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputOptics(idx, :, :, :) = reconImage;
    
    subplot(2, 3, idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
    
end

%% Diffraction-limited Optics
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

%% Run reconstruction
stimFreq = [2, 6, 12, 16, 20, 25, 32, 64, 96, 128];
rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 1.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

outputDiflmt = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    gratingScene = ...
        createGratingScene(chromaDir, stimFreq(idx), 'fovDegs', fovDegs, 'pixelsNum', 2048);
    
    crst = 0.50;
    [theSceneSequence, ~] = gratingScene.compute(crst);
    
    % Compute response from scene
    [~, allCone] = retina.computeWithScene(theSceneSequence{:});
    % visualizeOpticalImage(retina.LastOI, 'displayRadianceMaps', true, ...
    %    'displayRetinalContrastProfiles', true);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputDiflmt(idx, :, :, :) = reconImage;
    
    subplot(2, 4, idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
    
end
