load('sparsePrior.mat');
load('ishihara_9_linear.mat');

load('render_normal.mat');
load('render_anomalous.mat');
load('render_deuteranopia.mat');

renders = {5 * renderNormal, 5 * renderAnomalous, 5 * renderDeuteranopia};
results = zeros(3, 128, 128, 3);

regPara = 0.04;
parfor idx = 1 : 3
    fprintf('Reconstruction for Render %d \n', idx);
    
    render  = renders{idx};
    estimator = PoissonSparseEstimator(render, inv(regBasis), MU', regPara, 4, imageSize);
    
    coneVec = render * patchLinear(:);
    results(idx, :, :, :) = estimator.estimate(poissrnd(coneVec), 1.5e3, rand([prod(imageSize), 1]), true);
end

save('ishihara_9.mat', 'results', '-v7.3');