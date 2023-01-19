function pr = prFromBase(prBase,index,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
    stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
    forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor)

% Create pr structure from base structure plus filling in parameters in the
% passed arrays,according to index

pr = prBase;
pr.stimSizeDegs = stimSizeDegs(index);
pr.stimRVal = stimRVal(index);
pr.stimGVal = stimGVal(index);
pr.stimBVal = stimBVal(index);
pr.stimCenter = stimCenter(:,index);
pr.forwardDefocusDiopters = forwardDefocusDiopters(index);
pr.reconDefocusDiopters = reconDefocusDiopters(index);
pr.regPara = regPara(index);
pr.forwardChrom = forwardChrom(index);
pr.reconChrom = reconChrom(index);
pr.forwardPupilDiamMM = forwardPupilDiamMM(index);
pr.reconPupilDiamMM = reconPupilDiamMM(index);
pr.displayScaleFactor = displayScaleFactor(index);

% Make sure some new fields have a value if not specified
% by the caller.
if ~isfield(pr,'inputImageScaleFactor')
    pr.inputImageScaleFactor = 1;
end


end