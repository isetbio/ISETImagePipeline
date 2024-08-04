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
%
% This switch sorts out what the name of the display variable
% in the read in file is.  In some ideal world it would always
% be yoked to the filename, or always be the same, but neither
% of these is the case.
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

% Determine Poisson Noise string
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

% Initiate vectors and override default values with desired region
% proportion and variant 
regionVariant = pr.regionVariant; 
propL = pr.propL;
propS = pr.propS;

switch pr.focalRegion
    case 'center'
        propL(1) = pr.focalPropL;
        regionVariant(1) = pr.focalVariant;
    case 'nearSurround'
        propL(2) = pr.focalPropL;
        regionVariant(2) = pr.focalVariant;
    case 'distantSurround'
        propL(3) = pr.focalPropL;
        regionVariant(3) = pr.focalVariant;
    case 'multiple'
        if length(pr.focalPropL) ~= 3 || length(pr.focalVariant) ~=3
            error(['Must provide proportions and variants for all three regions'])
        end
        propL = pr.focalPropL;
        regionVariant = pr.focalVariant;
    case 'global'% This needs to be a list of 3 or the name will break
        propL = pr.focalPropL * ones(1,3);
        regionVariant = pr.focalVariant * ones(1,3);
    otherwise
        error(['Unrecognized focal region entered'])
end

%% Build render structure directories
if (pr.forwardAORender)
    cnv.forwardRenderDirFirst = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM), ...
        pr.forwardDefocusDiopters,cnv.forwardSeedStr, pr.forwardEccVars, pr.forwardNoLCA);
else
    cnv.forwardRenderDirFirst = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM),...
        pr.forwardDefocusDiopters,pr.forwardZernikeDataBase,pr.forwardSubjectID, cnv.forwardSeedStr,...
        pr.forwardEccVars, pr.forwardNoLCA);
end
if (pr.reconAORender)
    cnv.reconRenderDirFirst = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM), ...
        pr.reconDefocusDiopters, cnv.reconSeedStr, pr.reconEccVars, pr.reconNoLCA);
else
    cnv.reconRenderDirFirst = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM),...
        pr.reconDefocusDiopters, pr.reconZernikeDataBase,pr.reconSubjectID, cnv.reconSeedStr,...
        pr.reconEccVars, pr.reconNoLCA);
end

% Subdirectory naming based on aspects relevant to the small_quads project,
% for now this is constant between forward/recon conditions. 
cnv.renderDirSecond = sprintf('%0.1fArcmin', ...
    60*pr.stimSizeDegs);
cnv.renderDirThird = sprintf('regionVariant_v%d_v%d_v%d_%0.2f', ...
    regionVariant(1),regionVariant(2),regionVariant(3),60*pr.forwardOpticalBlurStimSizeExpansionDegs);

% Build the nested render directories for forward and recon conditions
cnv.forwardRenderDirFull = fullfile(pr.aoReconDir,pr.versEditor,'xRenderStructures', ...
    cnv.forwardRenderDirFirst,cnv.renderDirSecond,cnv.renderDirThird);
if (~exist(cnv.forwardRenderDirFull,'dir'))
    mkdir(cnv.forwardRenderDirFull);
end

cnv.reconRenderDirFull = fullfile(pr.aoReconDir,pr.versEditor,'xRenderStructures', ...
    cnv.reconRenderDirFirst,cnv.renderDirSecond,cnv.renderDirThird);
if (~exist(cnv.reconRenderDirFull,'dir'))
    mkdir(cnv.reconRenderDirFull);
end

% The actual file name is set to be the proportions since this is the most 
% pertinent information when dealing with small_quads render structures.
cnv.renderName = sprintf(['regionProps_%0.2fL_%0.2fL_%0.2fL_' ...
    '%0.2fS_%0.2fS_%0.2fS.mat'],propL(1),propL(2),propL(3), ...
    propS(1),propS(2),propS(3));


%% Build mosaic montage directories
% Build the nested render directories for forward and recon conditions
cnv.forwardMontageDirFull = fullfile(pr.aoReconDir,pr.versEditor,'xRenderStructures', ...
    cnv.forwardRenderDirFirst,'mosaicMontages');
if (~exist(cnv.forwardMontageDirFull,'dir'))
    mkdir(cnv.forwardMontageDirFull);
end

cnv.reconMontageDirFull = fullfile(pr.aoReconDir,pr.versEditor,'xRenderStructures', ...
    cnv.reconRenderDirFirst,'mosaicMontages');
if (~exist(cnv.reconMontageDirFull,'dir'))
    mkdir(cnv.reconMontageDirFull);
end

%% Build output directories
%
% Nested directory chain going from general to most pertinent things to the
% small quads routine, organized in a way that makes post-processing most
% straightforward. This organization may not be the best suited for other
% projects, but if so can expand from here (maybe make it versEditor
% dependent).

% General conditions are those that were found during the initial parameter
% search and have since been constant throughout subsequent simulations
cnv.generalConditions = sprintf(['%s_%s_%0.2f_%0.2f_%d_%d_%s_%0.2f_' ...
    '%s_%0.1f_%0.1f_%0.6f_%d_%s_%d_%d_%d_%d_%s_%d'], ...
    cnv.forwardAOStr,cnv.reconAOStr,pr.forwardDefocusDiopters, ...
    pr.reconDefocusDiopters,pr.nPixels,pr.fieldSizeMinutes,pr.displayName, ...
    pr.displayScaleFactor,pr.sparsePriorStr,pr.eccXDegs,pr.eccYDegs, ...
    pr.regPara,pr.stride,cnv.exciteSource,pr.forwardEccVars, ...
    pr.reconEccVars,pr.forwardNoLCA,pr.reconNoLCA,noiseStr,pr.boundedSearch);

cnv.outputDirFirst = sprintf(['%0.1fArcmin_%s'], ...
    60*pr.stimSizeDegs,pr.focalRegion);
cnv.outputDirSecond = sprintf(['regionVariant_v%d_v%d_v%d_%0.2f'], ...
    regionVariant(1),regionVariant(2),regionVariant(3),60*pr.forwardOpticalBlurStimSizeExpansionDegs);
cnv.outputDirThird = sprintf(['regionProps_%0.2fL_%0.2fL_%0.2fL_' ...
    '%0.2fS_%0.2fS_%0.2fS'],propL(1),propL(2),propL(3), ...
    propS(1),propS(2),propS(3));
cnv.outputDirFourth = sprintf(['stimColor_%0.4f_%0.4f_%0.4f_%0.4f_' ...
    'stimPosition_%d_%d'], ...
    pr.stimBgVal(1),pr.stimrVal,pr.stimgVal,pr.stimbVal, ...
    pr.stimCenter(1), pr.stimCenter(2)); 

cnv.outputDirFull = fullfile(pr.aoReconDir, pr.versEditor, cnv.generalConditions,...
    cnv.outputDirFirst, cnv.outputDirSecond, cnv.outputDirThird, cnv.outputDirFourth);

% Add a component to make output directory for summary figures as well
cnv.outputDirSummaryFigs = fullfile(pr.aoReconDir, pr.versEditor, cnv.generalConditions,...
    cnv.outputDirFirst, cnv.outputDirSecond, 'summaryFigs');
if (~exist(cnv.outputDirSummaryFigs,'dir'))
    mkdir(cnv.outputDirSummaryFigs);
end
end