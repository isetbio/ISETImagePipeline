function pr = prFromBase(prBase,index,stimSizeDegs,mosaicStimSizeDegs,stimSizePixels,stimrVal,stimgVal,stimbVal, ...
    stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
    forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor,focalRegion, ...
    focalPropL,focalVariant)

% Create pr structure from base structure plus filling in parameters in the
% passed arrays,according to index

pr = prBase;
pr.stimSizeDegs = stimSizeDegs(index);
pr.stimSizePixels = stimSizePixels(index);
pr.stimrVal = stimrVal(index);
pr.stimgVal = stimgVal(index);
pr.stimbVal = stimbVal(index);
pr.stimCenter = stimCenter(:,index);
pr.forwardDefocusDiopters = forwardDefocusDiopters(index);
pr.reconDefocusDiopters = reconDefocusDiopters(index);
pr.regPara = regPara(index);
pr.forwardPupilDiamMM = forwardPupilDiamMM(index);
pr.reconPupilDiamMM = reconPupilDiamMM(index);
pr.displayScaleFactor = displayScaleFactor{index};
pr.focalRegion = focalRegion(index);
pr.focalPropL = focalPropL(index);
pr.focalVariant = focalVariant(index);
pr.mosaicStimSizeDegs = mosaicStimSizeDegs(index);

% Set propS field based on stimSizeDegs. This
% is done to allow us to drive the number of S
% cones in small focal stimulated regions.
if (pr.stimSizeDegs < prBase.targetSizeSPropThresholdDegs)
    pr.propS = prBase.propSSmallTarget;
else
    pr.propS = prBase.propSLargeTarget;
end

% Make sure some new fields have a value if not specified
% by the caller.
if ~isfield(pr,'inputImageScaleFactor')
    pr.inputImageScaleFactor = 1;
end


end