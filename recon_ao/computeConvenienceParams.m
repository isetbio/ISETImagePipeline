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
    sum(pr.displayScaleFactor),pr.sparsePriorStr,pr.eccXDegs,pr.eccYDegs, ...
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
summaryFigsString = 'summaryFigsX';
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

%% Build summary figs routine output information
%
% Need to tell it where it should write the summary figure output.
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

%% Build grabImageInfo routine information
%
% Need to tell it where to loop over conditions that it will summarize
%
% This directory is in same output cascade and ends one level up from
% where the recons get written for each mosaic/stim combination that we
% want to summarize. This takes us to the level of the propL values which
% will be cycled through, and then within that the grabImageInfo.m routine
% will cycle through each of the stimulus colors.
levelMosaicSubdirs = length(indMosaicSubdirs);

% Changing this to 2 for now while renaming files, this should be 1 though
levelStimSubdirs = length(indStimSubdirs)-2;
nameMosaicSubdirs = strjoin(vars(indMosaicSubdirs(1:levelMosaicSubdirs)), ',');
nameStimSubdirs = strjoin(vars(indStimSubdirs(1:levelStimSubdirs)), ',');
eval(['cnv.outputSubdirImageInfo = fullfile(pr.aoReconDir, pr.versEditor,' ...
    'cnv.outputDirGeneral,' nameMosaicSubdirs ',' nameStimSubdirs ...
    ');']);

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

%% Patch to put old naming scheme and current naming scheme in agreement.
%
% Desired updates are as follows:
% Last level change to: stimColor_r_g_b
% Second to last level change to: centerProps_#L
% All of the names are the same exact length, so can afford to do something
% hacky and pull out exact index values to piecemeal the right name back
% together.

outputSubdirRegionProps = dir(cnv.outputSubdirImageInfo);
subDirNames = extractfield(outputSubdirRegionProps, 'name');
subDirNamesFullString = cell2mat(subDirNames);
if contains(subDirNamesFullString, '.DS_Store')
    startDirIndex = 4;
else
    startDirIndex = 3;
end


for i = startDirIndex:length(outputSubdirRegionProps)
    if outputSubdirRegionProps(i).isdir
        if contains(outputSubdirRegionProps(i).name, 'regionProp')

            % Go one step further first and rename all the values
            % within this folder


            outputSubdirStimColor = dir(fullfile(outputSubdirRegionProps(i).folder, ...
                outputSubdirRegionProps(i).name));

            subDirNames2 = extractfield(outputSubdirStimColor, 'name');
            subDirNamesFullString2 = cell2mat(subDirNames2);
            if contains(subDirNamesFullString2, '.DS_Store')
                startDirIndex2 = 4;
            else
                startDirIndex2 = 3;
            end

            for j = startDirIndex2:length(outputSubdirStimColor)
                if outputSubdirStimColor(j).isdir
                    updatedName2 = [outputSubdirStimColor(j).name(1:10) ...
                        outputSubdirStimColor(j).name(18:37)];
                    oldPath2 = fullfile(outputSubdirStimColor(j).folder, ...
                        outputSubdirStimColor(j).name);
                    updatedPath2 = fullfile(outputSubdirStimColor(j).folder, ...
                        updatedName2);
                    movefile(oldPath2,updatedPath2)
                end
            end



            % Fix the name of the prop directory holding the stim info
            updatedName = ['centerProps' outputSubdirRegionProps(i).name(12:17)];
            oldPath = fullfile(outputSubdirRegionProps(i).folder, ...
                outputSubdirRegionProps(i).name);
            updatedPath = fullfile(outputSubdirRegionProps(i).folder, ...
                updatedName);
            movefile(oldPath,updatedPath)

        end
    end
end






% cnv.renderDirUpdate = fullfile(pr.aoReconDir, pr.versEditor, ...
%     xRenderString);
% 
% outputSubdirGeneral = dir(cnv.renderDirUpdate);
% subDirNames = extractfield(outputSubdirGeneral, 'name');
% subDirNamesFullString = cell2mat(subDirNames);
% if contains(subDirNamesFullString, '.DS_Store')
%     startDirIndex = 4;
% else
%     startDirIndex = 3;
% end
% 
% 
% 
% for i = startDirIndex:length(outputSubdirGeneral)
%     if outputSubdirGeneral(i).isdir
% 
%         outputSubdirSize = dir(fullfile(outputSubdirGeneral(i).folder, ...
%             outputSubdirGeneral(i).name));
% 
%         subDirNames2 = extractfield(outputSubdirSize, 'name');
%         subDirNamesFullString2 = cell2mat(subDirNames2);
%         if contains(subDirNamesFullString2, '.DS_Store')
%             startDirIndex2 = 4;
%         else
%             startDirIndex2 = 3;
%         end
% 
% 
%         for j = startDirIndex2:length(outputSubdirSize)
%             if outputSubdirSize(j).isdir
% 
%                 outputSubdirRegionProps = dir(fullfile(outputSubdirSize(j).folder, ...
%                     outputSubdirSize(j).name));
% 
%                 subDirNames3 = extractfield(outputSubdirRegionProps, 'name');
%                 subDirNamesFullString3 = cell2mat(subDirNames3);
%                 if contains(subDirNamesFullString3, '.DS_Store')
%                     startDirIndex3 = 4;
%                 else
%                     startDirIndex3 = 3;
%                 end
% 
% 
%                 for k = startDirIndex3:length(outputSubdirRegionProps)
%                     if outputSubdirRegionProps(k).isdir
% 
%                         outputSubdirVariant = dir(fullfile(outputSubdirRegionProps(k).folder, ...
%                             outputSubdirRegionProps(k).name));
% 
%                         subDirNames4 = extractfield(outputSubdirVariant, 'name');
%                         subDirNamesFullString4 = cell2mat(subDirNames4);
%                         if contains(subDirNamesFullString4, '.DS_Store')
%                             startDirIndex4 = 4;
%                         else
%                             startDirIndex4 = 3;
%                         end
% 
%                         for m = startDirIndex4:length(outputSubdirVariant)
% 
%                             outputSubdirCenterProps = dir(fullfile(outputSubdirVariant(m).folder, ...
%                                 outputSubdirVariant(m).name));
% 
%                             % Fix the name of the prop directory holding the stim info
%                             updatedName = ['centerProps' outputSubdirCenterProps(m).name(12:17) '.mat'];
%                             oldPath = fullfile(outputSubdirCenterProps(m).folder, ...
%                                 outputSubdirCenterProps(m).name);
%                             updatedPath = fullfile(outputSubdirCenterProps(m).folder, ...
%                                 updatedName);
%                             movefile(oldPath,updatedPath)
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% 
% 
% 
% end