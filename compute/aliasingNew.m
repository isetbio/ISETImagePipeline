%% PART 1: foveal mosaic

%% Setup a cone mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

fovDegs = 0.20;
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

retina.visualizeMosaic();

%% Regular optics
retina.PSF = oiCreate('human', pupilSize);

render = retina.forwardRender(imageSize, true, true, false);
render = double(render);

%% Construct image estimator
regPara = 5e-3; stride = 4;
useGPU = true;
try
    estimator = PoissonSparseEstimator(gpuArray(single(render)), ...
        inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
catch NoGPU
    estimator = PoissonSparseEstimator(render, ...
        inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
    useGPU = false;
end

%% Run reconstruction with special scene
stimFreq = [8, 16, 32, 64, 96, 110, 120, 128];
rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

outputOptics = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    gratingScene = ...
        createGratingScene(chromaDir, stimFreq(idx), 'fovDegs', fovDegs, 'pixelsNum', 2048);
    
    crst = 1.0;
    [theSceneSequence, ~] = gratingScene.compute(crst);
    scene = theSceneSequence{:};
    
    % Special scene with narrow bandwith
    scene.data.photons(:, :, ...
        (scene.spectrum.wave ~= 600 & scene.spectrum.wave ~= 620 & scene.spectrum.wave ~= 640)) = 0.0;
    
    % Compute response from scene
    [~, allCone] = retina.computeWithScene(scene);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 200, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputOptics(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% Turn off optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

%% Run reconstruction with special scene
stimFreq = [8, 16, 32, 64, 96, 110, 120, 128];
rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

outputDiflmt = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    gratingScene = ...
        createGratingScene(chromaDir, stimFreq(idx), 'fovDegs', fovDegs, 'pixelsNum', 2048);
    
    crst = 1.0;
    [theSceneSequence, ~] = gratingScene.compute(crst);
    scene = theSceneSequence{:};
    
    % Special scene with narrow bandwith
    scene.data.photons(:, :, ...
        (scene.spectrum.wave ~= 600 & scene.spectrum.wave ~= 620 & scene.spectrum.wave ~= 640)) = 0.0;
    
    % Compute response from scene
    [~, allCone] = retina.computeWithScene(scene);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 200, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputDiflmt(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% PART 2: new mosaic
% Calculation with the new mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');
fovDegs = 1.0;

% Generate cone mosaic - [eccX, eccY] deg ecc
eccX = 18.0; eccY = 18.0;
retina = ConeResponseCmosaic...
    (eccX, eccY, 'fovealDegree', fovDegs, 'pupilSize', 3.0, 'subjectID', 6);

render = retina.forwardRender(imageSize, true, false);
render = double(render);

%% Construct image estimator
regPara = 5e-3; stride = 4;
useGPU = true;
try
    estimator = PoissonSparseEstimator(gpuArray(single(render)), ...
        inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
catch NoGPU
    estimator = PoissonSparseEstimator(render, ...
        inv(prior.regBasis), prior.mu', regPara, stride, imageSize);
    useGPU = false;
end

%% Run reconstruction with special scene
stimFreq = [6, 10, 16, 20, 25, 32];
rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

outputOptics = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    gratingScene = ...
        createGratingScene(chromaDir, stimFreq(idx), 'fovDegs', fovDegs, 'pixelsNum', 2048);
    
    crst = 1.0;
    [theSceneSequence, ~] = gratingScene.compute(crst);
    scene = theSceneSequence{:};
    
    % Special scene with narrow bandwith
    scene.data.photons(:, :, ...
        (scene.spectrum.wave ~= 600 & scene.spectrum.wave ~= 620 & scene.spectrum.wave ~= 640)) = 0.0;
    
    % Compute response from scene
    allCone = retina.computeWithScene(scene);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 200, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputOptics(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% Turn off optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

%% Run reconstruction with special scene
stimFreq = [6, 10, 16, 20, 25, 32];
rmsContrast = 1.0;
chromaDir = [1.0, 1.0, 0.0]';
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;

outputDiflmt = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    gratingScene = ...
        createGratingScene(chromaDir, stimFreq(idx), 'fovDegs', fovDegs, 'pixelsNum', 2048);
    
    crst = 1.0;
    [theSceneSequence, ~] = gratingScene.compute(crst);
    scene = theSceneSequence{:};
    
    % Special scene with narrow bandwith
    scene.data.photons(:, :, ...
        (scene.spectrum.wave ~= 600 & scene.spectrum.wave ~= 620 & scene.spectrum.wave ~= 640)) = 0.0;
    
    % Compute response from scene
    allCone = retina.computeWithScene(scene);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 200, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputDiflmt(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end
