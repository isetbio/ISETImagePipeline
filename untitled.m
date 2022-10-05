function pr = prFromBase(prBase,index,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
    stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
    forwardChrom,reconChrom)
% Create pr structure from base structure plus filling in parameters in the
% passed arrays

pr = prBase;
pr.stimSizeDegs = stimSizeDegs(pp);
pr.stimRVal = stimRVal(pp);
pr.stimGVal = stimGVal(pp);
pr.stimBVal = stimBVal(pp);
pr.stimCenter = stimCenter(:,pp);
pr.forwardDefocusDiopters = forwardDefocusDiopters(pp);
pr.reconDefocusDiopters = reconDefocusDiopters(pp);
pr.regPara = regPara(pp);
pr.forwardChrom = forwardChrom(pp);
pr.reconChrom = reconChrom(pp);

end