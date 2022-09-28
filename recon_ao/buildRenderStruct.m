function renderStructure = buildRenderStruct(aoReconDir, helpDir, ...
    eccXDegs, eccYDegs, fieldSizeDegs, nPixels, pupilDiamMM, aoRender, ...
    defocusDiopters, overwriteDisplayGamma, displayName, displayFieldName, ...
    displayGammaBits, displayGammaGamma, randSeed, replaceCones, startCones, ...
    newCones, eccVars)
% Synopsis:
%    Build render matrix if desired/needed
%
% Description:
%    Run this function if you would like to rebuild a new mosaic and 
%    render matrix.  This also gets run if there is no cached file corresponding
%    to the desired parameters. Once built, this file can be loaded from cache
%    for quicker running.
%
% See also: aoStimRecon, aoStimReconRunMany

% History:
%   08/16/22  chr  Made it a callable function
%   08/25/22  chr  Included portion for dichromacy
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options


% Get display
theDisplayLoad = load(fullfile(aoReconDir, helpDir, [displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
if (overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end
clear theDisplayLoad;

% Spline underlying display wavelength down to the wavelength sampling we
% will eventually use in the calculations.
wls = (400:10:700)';
theDisplay = displaySet(theDisplay,'wave',wls);

% Create and setup cone mosaic
% 
% For AO, we create a dummy object with 3 mm pupil and then adjust
% pupil and make the OI diffraction limited with no LCA.  The field
% name PSF is not optimal, because it is actually OI. We need the dummy
% 3 mm pupil because the code that pulls out the initial Polens optics
% checks that the desired pupil size is smaller than the 4 mm at which
% those data were measured.
if (aoRender)
    if (eccVars)
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'wave', wls, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', randSeed);
        theConeMosaic.PSF = ConeResponse.psfDiffLmt(pupilDiamMM,'defocusDiopters', defocusDiopters);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', wls, ... 
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true);
        theConeMosaic.PSF = ConeResponse.psfDiffLmt(pupilDiamMM,'defocusDiopters', defocusDiopters);
    end

% We build a normal optics structure. Allow specified defocus.
else
    if (eccVars)
        % Build normal optics structure. 
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'wave', wls, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters);
    else
       theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', wls, ... 
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true);
    end
end


% Option to replace cones in mosaic with another kind to simulate
% dichromacy. 
if (replaceCones)
    for i=1:length(startCones)
        coneInd = find(theConeMosaic.Mosaic.coneTypes == startCones(i));
        theConeMosaic.Mosaic.reassignTypeOfCones(coneInd, newCones);
    end
end

theConeMosaic.Display = theDisplay;

% Generate render matrix
renderMatrix = theConeMosaic.forwardRender([nPixels nPixels 3], ...
    true, true, 'useDoublePrecision', true);
renderMatrix = double(renderMatrix);

% Push new info back into structure and save
renderStructure.theDisplay = theDisplay;
renderStructure.renderMatrix = renderMatrix;
renderStructure.theConeMosaic = theConeMosaic;
renderStructure.fieldSizeDegs = fieldSizeDegs;
renderStructure.eccX = eccXDegs;
renderStructure.eccY = eccYDegs;
renderStructure.nPixels = nPixels;
renderStructure.pupilDiamMM = pupilDiamMM;
renderStructure.AORender = aoRender;
renderStructure.defocusDiopters = defocusDiopters;
end
    
