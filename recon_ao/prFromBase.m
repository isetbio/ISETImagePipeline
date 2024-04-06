function pr = prFromBase(prBase,index,stimSizeDegs,stimrVal,stimgVal,stimbVal, ...
    stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
    forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor,focalRegion, ...
    focalPropL, focalVariant)

% Create pr structure from base structure plus filling in parameters in the
% passed arrays,according to index

pr = prBase;
pr.stimSizeDegs = stimSizeDegs(index);
pr.stimrVal = stimrVal(index);
pr.stimgVal = stimgVal(index);
pr.stimbVal = stimbVal(index);
pr.stimCenter = stimCenter(:,index);
pr.forwardDefocusDiopters = forwardDefocusDiopters(index);
pr.reconDefocusDiopters = reconDefocusDiopters(index);
pr.regPara = regPara(index);
pr.forwardPupilDiamMM = forwardPupilDiamMM(index);
pr.reconPupilDiamMM = reconPupilDiamMM(index);
pr.displayScaleFactor = displayScaleFactor(index);
pr.focalRegion = focalRegion(index);
pr.focalPropL = focalPropL(index);
pr.focalVariant = focalVariant(index); 

% Make sure some new fields have a value if not specified
% by the caller.
if ~isfield(pr,'inputImageScaleFactor')
    pr.inputImageScaleFactor = 1;
end


end