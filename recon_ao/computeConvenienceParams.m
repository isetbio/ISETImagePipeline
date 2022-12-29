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
else
    cnv.forwardAOStr = ['NOAO' num2str(cnv.forwardPupilDiamMM) '_' pr.forwardZernikeDataBase '_' num2str(pr.forwardSubjectID)];
end

cnv.reconPupilDiamMM = pr.reconPupilDiamMM;
if (pr.reconAORender)
    cnv.reconAOStr = ['AO' num2str(cnv.reconPupilDiamMM)];
else
    cnv.reconAOStr = ['NOAO' num2str(cnv.reconPupilDiamMM) '_' pr.reconZernikeDataBase '_' num2str(pr.reconSubjectID)];
end

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

% Process mosaic chromatic type (trichromatic, deuteranopic, etc.) to be used.
% This provides variables that are used to process the mosaic below according
% to the specified chromatic type. See routine assignCones.
[cnv.replaceForwardCones, cnv.forwardStartCones, ...
    cnv.forwardNewCones] = assignCones(pr.forwardChrom);
[cnv.replaceReconCones, cnv.reconStartCones, ...
    cnv.reconNewCones] = assignCones(pr.reconChrom);

% Render structure name
if (pr.forwardAORender)
    cnv.forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%s_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM),pr.forwardDefocusDiopters,cnv.forwardSeedStr, pr.forwardChrom, pr.forwardEccVars);
else
    cnv.forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%s_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM),pr.forwardDefocusDiopters,pr.forwardZernikeDataBase,pr.forwardSubjectID, ...
        cnv.forwardSeedStr, pr.forwardChrom, pr.forwardEccVars);
end
if (pr.reconAORender)
    cnv.reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%s_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM),pr.reconDefocusDiopters, cnv.reconSeedStr, pr.reconChrom, pr.reconEccVars);
else
    cnv.reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%s_%d.mat', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM),pr.reconDefocusDiopters, pr.reconZernikeDataBase,pr.reconSubjectID, ...
        cnv.reconSeedStr, pr.reconChrom, pr.reconEccVars);
end

% Output directory name
if (pr.addPoissonNoise)
    noiseStr = 'noise';
else
    noiseStr = 'nonoise';
end

% Determine which quadrants are stimulated for small quads
if isfield(pr,'quadSelect')
    stimQuads = (find(pr.quadSelect' ~= 0));
    stimQuadsName = sprintf('%d',stimQuads');
else
    stimQuads = [];
    stimQuadsName = '';
end

cnv.outputMainName = sprintf('%s_%s_%s_%0.2f_%0.2f_%d_%d_%s_%s', ...
    pr.versEditor,cnv.forwardAOStr,cnv.reconAOStr,pr.forwardDefocusDiopters,pr.reconDefocusDiopters,pr.nPixels,pr.fieldSizeMinutes,pr.displayName,pr.sparsePriorStr);
if (length(pr.stimBgVal) > 1)
    cnv.outputSubName = sprintf('%0.1f_%0.6f_%d_%s_%s_%s_%d_%s_%d_%d_%d_%s_%d_%s', ...
        60*pr.stimSizeDegs, pr.regPara,pr.stride,pr.imageName,cnv.exciteSource, pr.forwardChrom, pr.forwardEccVars, pr.reconChrom, pr.reconEccVars, ...
        pr.stimCenter(1),pr.stimCenter(2),noiseStr,pr.boundedSearch, stimQuadsName);
else
    cnv.outputSubName = sprintf('%0.1f_%0.6f_%d_%0.4f_%0.4f_%0.4f_%0.4f_%s_%s_%d_%s_%d_%d_%d_%s_%d_%s', ...
        60*pr.stimSizeDegs, pr.regPara,pr.stride,pr.stimBgVal,pr.stimRVal,pr.stimGVal,pr.stimBVal, cnv.exciteSource, pr.forwardChrom, pr.forwardEccVars, pr.reconChrom, pr.reconEccVars, ...
        pr.stimCenter(1),pr.stimCenter(2),noiseStr,pr.boundedSearch, stimQuadsName);
end
cnv.outputDir = fullfile(pr.aoReconDir ,cnv.outputMainName,cnv.outputSubName);

end