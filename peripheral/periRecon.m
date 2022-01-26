% Reconstruction with a set of peripheral mosaic
% for comparison with the denoiser prior reconstruction

% LQZ, Jan 25, 2022

%% setup
imageSize = [100, 100, 3];
display = displayCreate('CRT12BitDisplay');
prior = load('sparsePrior.mat');

renderMtx = load('periRender.mat');
renderMtx = renderMtx.renderMtx;

load('inputLinear.mat');

%% run reconstruction
allOutput = cell(5, 1);
for rid = 1:5
    fprintf('rid %d / 5 \n', rid);
    render = double(renderMtx{rid});

    % Build an image reconstruction object with sparse prior
    regConst = 2.5e-3; stride = 4;
    estimator = PoissonSparseEstimator(render, inv(prior.regBasis), ...
        prior.mu', regConst, stride, imageSize);

    % Run reconstruction on cone response to each images
    % reconstructed images are in linear pixel space, need to
    % gamma correct them before visulization
    nIter = 1e3; optDisp = 'off';
    output = zeros(size(inputLinear));

    parfor idx = 1 : size(inputLinear, 1)
        input = reshape(inputLinear(idx, :, :, :), imageSize);
        coneResp = render * input(:);
        recon = estimator.runEstimate(coneResp, 'maxIter', nIter, 'display', optDisp);
        output(idx, :, :, :) = gammaCorrection(recon, display);
    end

    allOutput{rid} = output;
end

save('reconOut.mat', 'imageSize', 'allOutput', '-v7.3');
