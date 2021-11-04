%% Generate human cone mosaic with optics
imageSize = [100, 100, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display, 'pupilSize', 2.5);

%% Change all M cone to L cone
if cond == 0
    retina.reassignCone(0.0, retina.M_Cone_Idx, retina.L_Cone_Idx, false);

    [L, M, S] = retina.coneCount();
    fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

    retina.visualizeMosaic();
end

%% No S cone condition
if cond == 1
    retina.resetSCone();

    [L, M, S] = retina.coneCount();
    fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

    retina.visualizeMosaic();
end

%% Render matrix
render = retina.forwardRender(imageSize);

%% Reconstruction
load('./dichromacy/inputImage_100.mat');

% Build an image reconstruction object with sparse prior
regConst = 0.0; stride = 4;
estimator = ...
    PoissonSparseEstimator(render, inv(prior.regBasis), prior.mu', regConst, stride, imageSize);

% Run reconstruction on cone response to each images
% reconstructed images are in linear pixel space, need to 
% gamma correct them before visulization
nIter = 800; optDisp = 'final';
output = zeros(size(inputLinear));
for idx = 1:nImage
    input = reshape(inputLinear(idx, :, :, :), imageSize);
    coneResp = render * input(:);    
    recon = estimator.runEstimate(coneResp, 'maxIter', nIter, 'display', optDisp);
    output(idx, :, :, :) = gammaCorrection(recon, display);
end
