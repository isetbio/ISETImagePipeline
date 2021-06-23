function fullFieldRender(ecc, fovDegs, saveName, runRender, showPSF)
% Generate render matrix across eccentricity

% Parameters
eccX = ecc;
eccY = flip(eccX);

imgEdge = 128;
imageSize = [imgEdge, imgEdge, 3];

numX = length(eccX);
numY = length(eccY);
renderArray = cell(numX, numY);

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
            (eccX(idx), eccY(idy), 'fovealDegree', fovDegs, 'pupilSize', 3.0, 'subjectID', 9);
        
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
save(saveName, 'renderArray', ...
    'eccX', 'eccY', 'imageSize', '-v7.3');
end
