function renderStructure = buildRenderStruct_chr(aoReconDir, ...
    eccXDegs, eccYDegs, fieldSizeDegs, nPixels, pupilDiamMM, aoRender, ...
    defocusDiopters, overwriteDisplayGamma, displayName, displayFieldName, ...
    displayGammaBits, displayGammaGamma, randSeed)
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


% Get display
theDisplayLoad = load(fullfile(aoReconDir,[displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
if (overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end
clear theDisplayLoad;

% Create and setup cone mosaic
%
% For AO, we create a dummy object with 3 mm pupil and then adjust
% pupil and make the OI diffraction limited with no LCA.  The field
% name PSF is not optimal, because it is actually OI. We need the dummy
% 3 mm pupil because the code that pulls out the initial Polens optics
% checks that the desired pupil size is smaller than the 4 mm at which
% those data were measured.
if (aoRender)
    % We build a normal optics structure, and then overwrite the PSF
    % property.  Allow specified defocus.
    theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', randSeed);
    theConeMosaic.PSF = ConeResponse.psfDiffLmt(pupilDiamMM,'defocusDiopters', defocusDiopters);
else
    % Build normal optics structure. 
    theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
        'defocusDiopters',defocusDiopters);
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
    
