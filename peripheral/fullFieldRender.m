% Generate render matrix across eccentricity

% Parameters
eccX = -1.5 : 1 : 9.5;
eccY = flip(eccX);
fovDegs = 1.0;

imgEdge = 128;
imageSize = [imgEdge, imgEdge, 3];

display = displayCreate('CRT12BitDisplay');

numX = length(eccX);
numY = length(eccY);
renderArray = cell(numX, numY);

showPSF = false; runRender = true;
if showPSF
    plotAxis = tight_subplot(numX, numY, ...
        [.01 .01], [.01 .01], [.01 .01]);
end

% Loop through cone mosaics
for idx = 1:numX
    for idy = 1:numY
        fprintf("%d, %d \n", idx, idy);
        
        % Use subject 1, 2, 4, 7 for more reliable measurement
        retina = ConeResponseCmosaic...
            (eccX(idx), eccY(idy), 'fovealDegree', fovDegs, 'pupilSize', 3.0, 'subjectID', 6);
        
        if showPSF
            axes(plotAxis((idy - 1) * numY + idx));
            retina.visualizePSF();
            
            xlabel(''); ylabel('');
            xticklabels([]); yticklabels([]);
        end
        
        if runRender
            render = retina.forwardRender(imageSize, false, false);
            [~, ~, v] = svd(render, 'econ');
            
            renderArray{idx, idy} = single(v);
        end
    end
end

% Save results
save('renderArray.mat', 'renderArray', ...
    'eccX', 'eccY', 'imageSize', '-v7.3');