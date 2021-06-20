% Generate render matrix across eccentricity

% Parameters
eccX = -2.5 : 1 : 14.5;
eccY = 14.5 : -1 : -2.5;
fovDegs = 1.0;

imgEdge = 128;
imageSize = [imgEdge, imgEdge, 3];

display = displayCreate('CRT12BitDisplay');

numX = length(eccX);
numY = length(eccY);
renderArray = cell(numX, numY);

% Loop through cone mosaics
for idx = 1:numX
    for idy = 1:numY
        fprintf("%d, %d \n", idx, idy);
        retina = ConeResponseCmosaic...
            (eccX(idx), eccY(idy), 'fovealDegree', fovDegs, 'pupilSize', 3.0, 'subjectID', 6);
        
        render = retina.forwardRender(imageSize, false, false);
        [~, ~, v] = svd(render, 'econ');
        
        renderArray{idx, idy} = single(v);
    end
end

% Save results
save('renderArray.mat', 'renderArray', ...
    'eccX', 'eccY', 'imageSize', '-v7.3');