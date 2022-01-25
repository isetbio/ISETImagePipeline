%% Set up retinal mosaic
display = displayCreate('CRT12BitDisplay');

% Create a 1.0 deg foveal mosaic with a pupil size of 2.5 mm
% Note that to make the tutorial run faster, reduce the size of the mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display, 'pupilSize', 2.5, 'integrationTime', 0.1);
retina.visualizeMosaic();

imageSize = [128, 128, 3];

%% Cone ratio analysis
% generate mtxs
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0];

% figure();
for i = 1:length(s)
    for j = 1:length(l)
        lR = l(j) * (1 - s(i));
        mR = (1 - l(j)) * (1 - s(i));
        % scatter(lR, mR, '*r'); hold on;

        retina.setConeRatio(lR, mR);
        renderMtx = retina.forwardRender(imageSize, false, true, false);

        fileName = sprintf('mtx_%d_%d.mat', i, j);
        save(fullfile('.', 'allMtx', fileName), 'renderMtx', '-v7.3');

    end
end