%% Show test stimulus
display = displayCreate('CRT12BitDisplay');

stimPixel = 1024;
spatialFreq = 6;
ratio = 1;
nCycle = spatialFreq / ratio;

average = 0.5;
amplitude = 0.5;

stim = createStim(nCycle, average, amplitude, stimPixel, display);
imshow(stim, 'InitialMagnification', 100);

%% Generate cone mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

fovDegs = 0.20; ratio = 5;
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

% retina.resetSCone();
% retina.reassignSCone(0.1);
retina.visualizeMosaic();

%% Generate render matrix
retina.PSF = oiCreate('human', pupilSize);

render = retina.forwardRender(imageSize, true, true, false);
render = double(render);

%% Construct estimator
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

%% Run reconstruction
retina.PSF = oiCreate('human', pupilSize);

stimFreq = [0, 16, 32, 64, 96, 110, 128, 150];
stimPixel = 2048;

average = 0.5;
amplitude = 0.5;

outputOptics = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    nCycle = stimFreq(idx) / ratio;
    stim = createStim(nCycle, average, amplitude, stimPixel, display);
    
    % Compute response from scene
    [~, ~, ~, allCone] = retina.compute(stim);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputOptics(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% Turn off optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

% Run reconstruction
stimFreq = [0, 16, 32, 64, 96, 110, 128, 150];
stimPixel = 2048;

average = 0.05;
amplitude = 0.04;

outputDiflmt = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    nCycle = stimFreq(idx) / ratio;
    stim = createStim(nCycle, average, amplitude, stimPixel, display);
    
    % Compute response from scene
    [~, ~, ~, allCone] = retina.compute(stim);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputDiflmt(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% Create mosaic
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

%% Constant
display = displayCreate('CRT12BitDisplay');

stimPixel = 2048;
ratio = 1;

average = 0.25;
amplitude = 0.25;

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

%% Run reconstruction
stimFreq = [5, 25];
outputOptics = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    nCycle = stimFreq(idx) / ratio;
    stim = createStim(nCycle, average, amplitude, stimPixel, display);
    
    % Compute response from scene
    allCone = retina.compute(stim);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputOptics(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% Turn off optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

% Run reconstruction
stimFreq = [5, 15, 20, 22, 25];

average = 0.03;
amplitude = 0.015;

outputDiflmt = zeros([length(stimFreq), imageSize]);

figure();
for idx = 1:length(stimFreq)
    nCycle = stimFreq(idx) / ratio;
    stim = createStim(nCycle, average, amplitude, stimPixel, display);
    
    % Compute response from scene
    allCone = retina.compute(stim);
    
    reconImage = ...
        estimator.runEstimate(allCone, 'maxIter', 150, ...
        'display', 'iter', 'gpu', useGPU);
    
    outputDiflmt(idx, :, :, :) = reconImage;
    
    subplot(2, ceil(length(stimFreq)/2), idx);
    imshow(gammaCorrection(reconImage, display), 'InitialMagnification', 500);
end

%% Helper function
function stim = createStim(nCycle, average, amplitude, stimPixel, display)

signal = cos(linspace(0, 2 * pi, stimPixel) * nCycle) * amplitude + average;
% stim = repmat(repmat(signal, [stimPixel, 1]), [1, 1, 3]);

stim = zeros([stimPixel, stimPixel, 3]);
stim(:, :, 1) = repmat(signal, [stimPixel, 1]);

stim = gammaCorrection(stim, display);

end
