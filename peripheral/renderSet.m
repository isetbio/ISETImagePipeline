% Produce a set of render matrix
% with randomized cone types and location

nIter = 10;
retina = cell(nIter, 1);
render = cell(nIter, 1);

imageSize = [128, 128, 3];
for idx = 1:nIter
    
    cm = ConeResponseCmosaic(1.0, 1.0, ...
        'fovealDegree', 0.25, 'pupilSize', 3.0, 'subjectID', 9, 'randomMesh', true);
    
    mtx = retina.forwardRender(imageSize, false, false);
    
    retina{idx} = cm;
    render{idx} = mtx;
end