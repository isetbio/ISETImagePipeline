%% generate a cone mosaic
imageSize = [64, 64, 3];
display = load('display.mat');
prior   = load('sparsePrior.mat');

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display.CRT12BitDisplay, 'pupilSize', 2.5);

%% change the S cone proportion
% and generate the corresponding render matrix
retina.resetSCone();
ratio = [0, 0.01, 0.05, 0.1, 0.25, 0.50, 0.75, 0.9];

mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));
for idx = 1:length(ratio)
    % manipulate the number of S cone in the mosaic
    retina.resetSCone();
    retina.reassignSCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize, false, true, false);
    renderArray(idx) = {double(renderMtx)};
end

%% analysis without chromatic abberation
