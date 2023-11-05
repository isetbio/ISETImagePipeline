%% Generate cone mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display, 'pupilSize', 2.5);

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

%% Manipulate of S cone ratio
retina.resetCone();

ratio = [0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95];
mosaicArray = cell(1, length(ratio));
renderArray = cell(1, length(ratio));
for idx = 1:length(ratio)
    retina.resetSCone();
    fprintf('S Cone Ratio: %.4f \n', ratio(idx));
    
    retina.reassignSCone(ratio(idx));
    mosaicArray(idx) = {retina.Mosaic.pattern};
    
    % Generate render matrix for each cone mosaic
    renderMtx = retina.forwardRender(imageSize, false, true, false);
    renderArray(idx) = {renderMtx};
end

% Save Matrix
save('sMtxArray.mat', 'ratio', 'renderArray', '-v7.3');
