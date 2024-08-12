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

%% Set up some base variables
%
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

%% Name the output directories
%
% Nested directory chain going from general to most pertinent things to the
% small quads routine, organized in a way that makes post-processing most
% straightforward. This organization may not be the best suited for other
% projects, but if so can expand from here (maybe make it versEditor
% dependent).
%
% The code below is fancy in that it lets you define a series of mosaic
% based subdirectory names of form outputSubdirMosaicX, X = 1:N and
% a series of stimulus based subdirectory names of form outputSubdirStimX,
% X = 1:M, and then chains these all together with appropriate general
% conditions directories for recon output, forward render matrices and
% recon render matrices in a way that one could in principle add more or
% have fewer of the subdirs with various names and it will all still work.
% Does this through the magic of Matlab string processing. 

% General conditions directory for output files, largely constant
cnv.outputDirGeneral  = sprintf(['%s_%s_%0.2f_%0.2f_%d_%d_%s_%0.2f_' ...
    '%s_%0.1f_%0.1f_%0.6f_%d_%s_%d_%d_%d_%d_%s_%d'], ...
    cnv.forwardAOStr,cnv.reconAOStr,pr.forwardDefocusDiopters, ...
    pr.reconDefocusDiopters,pr.nPixels,pr.fieldSizeMinutes,pr.displayName, ...
    pr.displayScaleFactor,pr.sparsePriorStr,pr.eccXDegs,pr.eccYDegs, ...
    pr.regPara,pr.stride,cnv.exciteSource,pr.forwardEccVars, ...
    pr.reconEccVars,pr.forwardNoLCA,pr.reconNoLCA,noiseStr,pr.boundedSearch);

% General conditions directory for forward render matrices
if (pr.forwardAORender)
    cnv.forwardRenderDirGeneral = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM), ...
        pr.forwardDefocusDiopters,cnv.forwardSeedStr, pr.forwardEccVars, pr.forwardNoLCA);
else
    cnv.forwardRenderDirGeneral = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.forwardPupilDiamMM),...
        pr.forwardDefocusDiopters,pr.forwardZernikeDataBase,pr.forwardSubjectID, cnv.forwardSeedStr,...
        pr.forwardEccVars, pr.forwardNoLCA);
end

% General conditions directory for recon render matrices
if (pr.reconAORender)
    cnv.reconRenderDirGeneral = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_AO_%0.2f_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM), ...
        pr.reconDefocusDiopters, cnv.reconSeedStr, pr.reconEccVars, pr.reconNoLCA);
else
    cnv.reconRenderDirGeneral = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%s_NOAO_%0.2f_%s_%d_%s_%d_%d', ...
        pr.displayName,pr.fieldSizeMinutes,pr.eccXDegs,pr.eccYDegs,pr.nPixels,num2str(cnv.reconPupilDiamMM),...
        pr.reconDefocusDiopters, pr.reconZernikeDataBase,pr.reconSubjectID, cnv.reconSeedStr,...
        pr.reconEccVars, pr.reconNoLCA);
end

% Output sublevels pertaining to mosaic properties.  See comment below that
% notes that there is one mosaic property that is indicated lower down
% amidst the stimulus subdirectories.
outputSubdirMosaic1 = sprintf(['%0.1fArcmin_%sRender'], ...
    60*pr.stimSizeDegs, pr.focalRegion);

outputSubdirMosaic2 = sprintf(['regionProps_%0.2fL_%0.2fL_' ...
    '%0.2fS_%0.2fS_%0.2fS'],propL(2),propL(3), ...
    propS(1),propS(2),propS(3));

outputSubdirMosaic3 = sprintf(['regionVariant_v%d_v%d_v%d_%0.2f'], ...
    regionVariant(1),regionVariant(2),regionVariant(3), ...
    60*pr.forwardOpticalBlurStimSizeExpansionDegs);

% Output sublevels pertaining to stimulus properties.  Strictly speaking
% these are not all determined by stimulus properties, as we include
% the centerPropsL variation down in this chain. But it just snuck in
% there; most of these are in fact stimulus determined.  The reason it
% snuck in is that our main summarization code is interested in the cross
% between focal region propL and stimulus in the reconstructions.
outputSubdirStim1 = sprintf(['bgColor_%0.4f_%0.4f_%0.4f_' ...
    'stimSeriesVariant_%d_' ...
    'stimPosition_%d_%d'], ...
    pr.stimBgVal(1), pr.stimBgVal(2), pr.stimBgVal(3), ...
    pr.stimSeriesVariant, ...
    pr.stimCenter(1), pr.stimCenter(2));

outputSubdirStim2 = sprintf(['%0.1fArcmin_%sStim'], ...
    60*pr.stimSizeDegs, pr.focalRegion);

outputSubdirStim3 = sprintf(['centerPropsL_%0.2fL'], ...
    propL(1));

outputSubdirStim4 = sprintf(['stimColor_%0.4f_%0.4f_%0.4f'], ...
    pr.stimrVal,pr.stimgVal,pr.stimbVal);

% Output sublevel strings pertaining to summary figs and render matrices
summaryFigsString = 'summaryFigs'; 
xRenderString = 'xRenderStructures';

% Render matrix .mat file name
cnv.renderName = sprintf(['regionProps_%0.2fL_%0.2fL_%0.2fL_' ...
    '%0.2fS_%0.2fS_%0.2fS.mat'],propL(1),propL(2),propL(3), ...
    propS(1),propS(2),propS(3));

%% Build recon computation output directory name
%
% Find all instances of subdirectory naming from above and capture. More
% effort in creating this framework but also more resilient to changes in
% directory order / additional layers
%
% This will break because the outputSubdir values will need to be fed into
% the cnv structure, should replace with a contains(fieldnames(cnv), 'output')
vars = who();
indMosaicSubdirsAll = contains(vars, 'outputSubdirMosaic');
indStimSubdirsAll = contains(vars, 'outputSubdirStim');

indMosaicSubdirs = find(indMosaicSubdirsAll == 1);
indStimSubdirs = find(indStimSubdirsAll == 1);

% Fullfile name for the simulation outputs, with identification of how many
% layers of the full option set you want to go down 
levelMosaicSubdirs = length(indMosaicSubdirs); 
levelStimSubdirs = length(indStimSubdirs);

nameMosaicSubdirs = strjoin(vars(indMosaicSubdirs(1:levelMosaicSubdirs)), ',');
nameStimSubdirs = strjoin(vars(indStimSubdirs(1:levelStimSubdirs)), ',');

eval(['cnv.outputDirFull = fullfile(pr.aoReconDir, pr.versEditor,' ...
    'cnv.outputDirGeneral,' nameMosaicSubdirs ',' nameStimSubdirs ');']);

%% Build summary figs routine information
%
% Need to tell it where to loop over conditions that it will summarize
% as well as where it should write the summary figure output.
%
% This directory is in same output cascade and ends two levels up from
% where the recons get written for each mosaic/stim combination that we
% want to summarize. We find this this by loping off 2 of the "stim"
% subdirs, which in fact correspond to the mosaic focal region propL and
% the stimulus set being run.
levelMosaicSubdirs = length(indMosaicSubdirs); 
levelStimSubdirs = length(indStimSubdirs)-2;
nameMosaicSubdirs = strjoin(vars(indMosaicSubdirs(1:levelMosaicSubdirs)), ',');
nameStimSubdirs = strjoin(vars(indStimSubdirs(1:levelStimSubdirs)), ',');
eval(['cnv.outputSubdirSummaryFigs = fullfile(pr.aoReconDir, pr.versEditor,' ...
    'cnv.outputDirGeneral,' nameMosaicSubdirs ',' nameStimSubdirs ...
    ', summaryFigsString);']);
if (~exist(cnv.outputSubdirSummaryFigs,'dir'))
    mkdir(cnv.outputSubdirSummaryFigs);
end

%% Build the nested render matrix directories for forward and recon
levelMosaicSubdirs = length(indMosaicSubdirs); 
nameMosaicSubdirs = strjoin(vars(indMosaicSubdirs(1:levelMosaicSubdirs)), ',');
eval(['cnv.forwardRenderDirFull = fullfile(pr.aoReconDir, pr.versEditor,' ...
    'xRenderString, cnv.forwardRenderDirGeneral,' nameMosaicSubdirs ');']);
if (~exist(cnv.forwardRenderDirFull,'dir'))
    mkdir(cnv.forwardRenderDirFull);
end

eval(['cnv.reconRenderDirFull = fullfile(pr.aoReconDir, pr.versEditor,' ...
    'xRenderString, cnv.reconRenderDirGeneral,' nameMosaicSubdirs ');']);
if (~exist(cnv.reconRenderDirFull,'dir'))
    mkdir(cnv.reconRenderDirFull);
end

% %% Build mosaic montage directories
% % Build the nested render directories for forward and recon conditions
% cnv.forwardMontageDirFull = fullfile(pr.aoReconDir,pr.versEditor,'xRenderStructures', ...
%     cnv.forwardRenderDirFirst,'mosaicMontages');
% if (~exist(cnv.forwardMontageDirFull,'dir'))
%     mkdir(cnv.forwardMontageDirFull);
% end
% 
% cnv.reconMontageDirFull = fullfile(pr.aoReconDir,pr.versEditor,'xRenderStructures', ...
%     cnv.reconRenderDirFirst,'mosaicMontages');
% if (~exist(cnv.reconMontageDirFull,'dir'))
%     mkdir(cnv.reconMontageDirFull);
% end




