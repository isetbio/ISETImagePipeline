% set up the variables
imageSize = [128, 128, 3];

display = displayCreate('CRT12BitDisplay');

prior = load('sparsePrior.mat');

dataset = load('input_cone_ratio.mat');
dataset = double(dataset.input_cone_ratio);

% run reconstruction for all matrix
nRow = 12;
nCol = 15;

allOutput = cell(nRow, nCol);

basePath = '/mnt/home/lzhang/ceph/render/allMtx';
for idx = 1:nRow
    for idy = 1:nCol
        fprintf('%d, %d \n', idx, idy);

        % load the matrix
        fileName = sprintf('mtx_%d_%d.mat', idx, idy);
        fullName = fullfile(basePath, fileName);
        mtx = load(fullName);
        mtx = double(mtx.renderMtx);

        % construct the optimizer        
        regConst = 1e-2; stride = 4;
        estimator = PoissonSparseEstimator(mtx, inv(prior.regBasis), ...
                            prior.mu', regConst, stride, imageSize);

        % parameters for the reconstruction
        nIter = 6e2; optDisp = 'off';
        output = zeros(size(dataset));

        % run reconstruction
        parfor iid = 1 : size(dataset, 1)
            input = reshape(dataset(iid, :, :, :), imageSize);

            coneResp = mtx * input(:);
            recon = estimator.runEstimate(coneResp, 'maxIter', nIter, 'display', optDisp);

            output(iid, :, :, :) = gammaCorrection(recon, display);
        end

        % save the result
        allOutput{idx, idy} = output;
    end
end

save('ConeRatioRecon.mat', 'allOutput', '-v7.3');
