function cnv = computeConvenienceParams(pr)
% Compute convenience parameter struct from parameter struct
%
% Syntax:
%     cnv = computeConvenienceParams(pr)
%
% Description:
%    Compute convenience parameters from the core parameters
%    The convenience parameters massage the parameters into
%    a more directly useful form for some parts of the code.

% Render directory
cnv.renderDir  = fullfile(pr.aoReconDir , pr.versEditor);
if (~exist(cnv.renderDir ,'dir'))
    mkdir(cnv.renderDir );
end

% Determine which method will be used for the reconstruction: ISETBIO or
% Render Matrix
if (pr.reconstructfromRenderMatrix)
    cnv.exciteSource = 'renderMatrix';
else
    cnv.exciteSource = 'isetbio';
end

% Field size in degrees from minutes
cnv.fieldSizeDegs = pr.fieldSizeMinutes/60;

% Determine AO/pupil rendering string
cnv.forwardPupilDiamMM = pr.forwardPupilDiamMM;
if (pr.forwardAORender)
    cnv.forwardAOStr = ['AO' num2str(cnv.forwardPupilDiamMM)];
    cnv.forwardAOStrPlt = ['AO' num2str(cnv.forwardPupilDiamMM)];
else
    cnv.forwardAOStr = ['NOAO' num2str(cnv.forwardPupilDiamMM) '_' pr.forwardZernikeDataBase '_' num2str(pr.forwardSubjectID)];
    cnv.forwardAOStrPlt = ['NOAO' num2str(cnv.forwardPupilDiamMM) ' ' pr.forwardZernikeDataBase ' ' num2str(pr.forwardSubjectID)];
end

cnv.reconPupilDiamMM = pr.reconPupilDiamMM;
if (pr.reconAORender)
    cnv.reconAOStr = ['AO' num2str(cnv.reconPupilDiamMM)];
    cnv.reconAOStrPlt = ['AO' num2str(cnv.reconPupilDiamMM)];
else
    cnv.reconAOStr = ['NOAO' num2str(cnv.reconPupilDiamMM) '_' pr.reconZernikeDataBase '_' num2str(pr.reconSubjectID)];
    cnv.reconAOStrPlt = ['NOAO' num2str(cnv.reconPupilDiamMM) ' ' pr.reconZernikeDataBase ' ' num2str(pr.reconSubjectID)];
end

% Set the display name
switch (pr.displayName)
    case 'conventional'
        cnv.displayFieldName = 'CRT12BitDisplay';
        cnv.overwriteDisplayGamma = false;
    case 'mono'
        cnv.displayFieldName = 'monoDisplay';
        cnv.overwriteDisplayGamma = true;
    otherwise
        error('Unknown display specified');
end

% Establish if a random seed will be used when building the cone mosaic for
% the render matrices
if (pr.forwardRandSeed)
    cnv.forwardSeedStr = 'rand';
else
    cnv.forwardSeedStr = 'noRand';
end
if (pr.reconRandSeed)
    cnv.reconSeedStr = 'rand';
else
    cnv.reconSeedStr = 'noRand';
end

%% TO-DO: UPDATE the Render structure name based on new approaches, maybe
% also consider nested directories with shorter names here. 
if (pr.forwardAORender)
    cnv.forwardRenderStructDir = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%s_%d_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM), ...
        pr.forwardDefocusDiopters,cnv.forwardSeedStr, pr.forwardEccVars, pr.forwardNoLCA);
else
    cnv.forwardRenderStructDir = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%s_%d_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM),...
        pr.forwardDefocusDiopters,pr.forwardZernikeDataBase,pr.forwardSubjectID, cnv.forwardSeedStr,...
        pr.forwardEccVars, pr.forwardNoLCA);
end
if (pr.reconAORender)
    cnv.reconRenderStructDir = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%s_%d_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM), ...
        pr.reconDefocusDiopters, cnv.reconSeedStr, pr.reconEccVars, pr.reconNoLCA);
else
    cnv.reconRenderStructDir = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%s_%d_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM),...
        pr.reconDefocusDiopters, pr.reconZernikeDataBase,pr.reconSubjectID, cnv.reconSeedStr,...
        pr.reconEccVars, pr.reconNoLCA);
end

% Naming sequence for the render matrix subdirectory, which as a
% simplification will be a static thing for our purposes between
% forward/recon conditions in the small_quads project
cnv.reconRenderStructName = sprintf('%s_');


%% Determine Poisson Noise string
if (pr.addPoissonNoise)
    noiseStr = 'noise';
else
    noiseStr = 'nonoise';
end

% Determine LCA string
if (pr.forwardNoLCA)
    cnv.forwardLCAStr = 'NoLCA';
else
    cnv.forwardLCAStr = 'LCA';
end

if (pr.reconNoLCA)
    cnv.reconLCAStr = 'NoLCA';
else
    cnv.reconLCAStr = 'LCA';
end

%% Updates
% Incorporating some portion to account for the new naming scheme and
% approach to running mosaics. Also do some general function formating
%
% Below is a hodgepodge directory name comprised of all the components that
% are expected to stay constant in the current small_quads state. For the
% sake of making output directories more immediately informative, will no
% longer be including all of this in the name but rather store it and group
% under some umbrella generalConditions_version1. Can then compare
% input parameters (which again should be unchanged) against the stored
% version and if a match can proceed down to make the actual output
% directories. If not, sends a flag that something has been altered and a
% new generalConditions_version should be established. Will become relevant
% if we decide to expand our parameter search again or boost pixel/FOV size
cnv.generalConditions_verion1 = sprintf(['%s_%s_%0.2f_%0.2f_%d_%d_%s_%0.2f_' ...
    '%s_%0.1f_%0.1f_%0.6f_%d_%s_%d_%d_%d_%d_%s_%d'], ...
    cnv.forwardAOStr,cnv.reconAOStr,pr.forwardDefocusDiopters, ...
    pr.reconDefocusDiopters,pr.nPixels,pr.fieldSizeMinutes,pr.displayName, ...
    pr.displayScaleFactor,pr.sparsePriorStr,pr.eccXDegs,pr.eccYDegs, ...
    pr.regPara,pr.stride, cnv.exciteSource,pr.forwardEccVars, ...
    pr.forwardNoLCA, pr.reconEccVars, pr.reconNoLCA, ...
    noiseStr,pr.boundedSearch);

% Throw error if see that some parameter that should have been constant is
% changed. 
error(["General parameter conditions unrecognized, adjust ouptut directories"]);


% TO-DO: Include a new constant param accounting for same mosaic structure on
% forward encoding and recon decoding since no longer have QS 
mosaicForwardReconSame = true; 


% TO-DO: Include a chunk to call out the region variant based on the focal variant
% of choice and the focal prop of choice

% Nested directory chain corresponding to the things most pertinent to the
% small quads routine, organized in a way that makes post-processing most
% straightforward. This organization may not be the best suited for other
% projects, but if so can expand from here. Maybe make it versEditor
% dependent. 
cnv.firstDir = sprintf(['stimSize_%0.1fArcmin_focalRegion_%s_stimPosition_' ...
    '%d_%d'], ...
    60*pr.stimSizeDegs,pr.focalRegion,pr.stimCenter(1),pr.stimCenter(2));
cnv.secondDir = sprintf(['regionProportions_%0.2fL_%0.2fL_%0.2fL_' ...
    '%0.2fS_%0.2fS_%0.2fS'],pr.propL(1),pr.propL(2),pr.propL(3), ...
    pr.propS(1),pr.propS(2),pr.propS(3));
cnv.thirdDir = sprintf(['regionVariant_%d_%d_%d'], ...
    pr.regionVariant(1),pr.regionVariant(2),pr.regionVariant(3));
cnv.fourthDir = sprintf(['stimColor_%0.4f_%0.4f_%0.4f_%0.4f'], ...
    pr.stimBgVal(1),pr.stimrVal,pr.stimgVal,pr.stimbVal); 

cnv.outputDir = fullfile(pr.aoReconDir, pr.versEditor, pr.system, ...
    cnv.firstDir, cnv.secondDir, cnv.thirdDir, cnv.fourthDir);

end