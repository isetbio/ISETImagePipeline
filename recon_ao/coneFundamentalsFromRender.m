%% Initialize
close all; clear all ieInit

%% Setup / Simulation parameters
nPixels = 58;
fieldSizeMinutes = 30;
fieldSizeDegs = fieldSizeMinutes/60;
eccXDegs = 2.0;
eccYDegs = 0.0;

% This will allow us to load in project specific precomputed information.
% Also records initials of version editors, otherwise set to 'main'
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
helpDir = '/helperFiles';
versEditor = '_coneFunds';

%% Create scenes with equal quantal monochromatic light
% 
% Input wavelength range and initialize stim input and fundamental output
% vectors
wls = 400:700;
fullStim = ones(3, length(wls));
coneQuantalFundamentals = fullStim; 

% Create the scenes and convert results to RGB values
for lambda = wls
    scene = sceneCreate('uniform monochromatic', lambda); 
    scene = sceneSet(scene, 'name', ['narrow band ' int2str(lambda)]);
%     sceneWindow(scene)
    stimulusImageRGB = sceneGet(scene, 'rgb image');
    stimulusImageRGB = stimulusImageRGB(1:nPixels, 1:nPixels, :);
    stimRGB = squeeze(stimulusImageRGB(1,1,:)); 
    fullStim(:, (lambda - min(wls) + 1)) = stimRGB;
end

%% Copy relevant portions of aoStimReconRunMany and aoStimRecon 
% 
% Line 58 - Loop begins, Line 254 - Create stimulus scene using converted 
% RGB values, Line 272 - Capture excitation quanta for cone classes

% Begin the copy of aoStimReconRunMany, only manipulation is of stimulus
% RGb values

fullStim = [0.80;0.65;0.10]; % Manual overwrite for now, remove later
displayName = 'mono';
stimSizeDegs = 24/60; 
stimBgVal = 0.1; stimRValList = fullStim(1,:); 
stimGValList = fullStim(2,:); stimBValList = fullStim(3,:);
sparsePriorStr = 'conventional';
regPara = 0.001;
stride = 2;

forwardAORender = true; reconAORender = true; 
forwardDefocusDiopters = 0; reconDefocusDiopters = 0;
forwardChrom = "chromNorm"; reconChrom = "chromNorm";

% Create the loop to go over RGB values at each wavelength
for cc = 1:length(stimRValList)
    stimRVal = stimRValList(cc);
    stimGVal = stimGValList(cc);
    stimBVal = stimBValList(cc);

    % Begin the copy from aoRecon
    aoReconDir = getpref('ISETImagePipeline','aoReconDir');
    helpDir = '/helperFiles';
    versEditor = '_coneFunds';

    % Gamma table parameters
    displayGammaBits = 12;
    displayGammaGamma = 2;
    switch (displayName)
        case 'conventional'
            displayFieldName = 'CRT12BitDisplay';
            overwriteDisplayGamma = false;
        case 'mono'
            displayFieldName = 'monoDisplay';
            overwriteDisplayGamma = true;
        otherwise
            error('Unknown display specified');
    end

    % Sparse prior name
    sparsePriorName = [sparsePriorStr 'SparsePrior.mat'];

    % Determine forward pupil diameter, allowing it to differ in AO case
    if (forwardAORender)
        forwardPupilDiamMM = 7;
        forwardAOStr = ['AO' num2str(forwardPupilDiamMM)];
    else
        forwardPupilDiamMM = 3;
        forwardAOStr = ['NOAO' num2str(forwardPupilDiamMM)];
    end
    
    if (reconAORender)
        reconPupilDiamMM = 7;
        reconAOStr = ['AO' num2str(reconPupilDiamMM)];
    else
        reconPupilDiamMM = 3;
        reconAOStr = ['NOAO' num2str(reconPupilDiamMM)];
    end

    % Force build and save
    buildNewForward = false;
    buildNewRecon = false;

    % Determine which method will be used for the reconstruction: ISETBIO or
    % Render Matrix
    reconstructfromRenderMatrix = true;
    if (reconstructfromRenderMatrix)
        exciteSource = 'renderMatrix';
    end

    % Establish if a random seed will be used when building the cone mosaic for
    % the render matrices
    forwardRandSeed = false;
    reconRandSeed = false;
    if (forwardRandSeed)
        forwardSeedStr = 'rand';
    else
        forwardSeedStr = 'noRand';
    end
    if (reconRandSeed)
        reconSeedStr = 'rand';
    else
        reconSeedStr = 'noRand';
    end

    % Optional chromaticities to be used in forward and recon cone mosaics,
    % with normal being trichromatic
    if strcmp(forwardChrom, 'chromProt')
        replaceForwardCones = true;
        forwardStartCones = 1;
        forwardNewCones = cMosaic.MCONE_ID;
    elseif strcmp(forwardChrom, 'chromDeut')
        replaceForwardCones = true;
        forwardStartCones = 2;
        forwardNewCones = cMosaic.LCONE_ID;
    elseif strcmp(forwardChrom, 'chromTrit')
        replaceForwardCones = true;
        forwardStartCones = 3;
        forwardNewCones = cMosaic.LCONE_ID;
    elseif strcmp(forwardChrom, 'chromAllL')
        replaceForwardCones = true;
        forwardStartCones = [2 3];
        forwardNewCones = cMosaic.LCONE_ID;
    elseif strcmp(forwardChrom, 'chromAllM')
        replaceForwardCones = true;
        forwardStartCones = [1 3];
        forwardNewCones = cMosaic.MCONE_ID;
    elseif strcmp(forwardChrom, 'chromNorm')
        replaceForwardCones = false; 
        forwardStartCones = [];
        forwardNewCones = [];
    else
        error('Unrecognized chromaticity input')
    end
    if strcmp(reconChrom, 'chromProt')
        replaceReconCones = true;
        reconStartCones = 1;
        reconNewCones = cMosaic.MCONE_ID;
    elseif strcmp(reconChrom, 'chromDeut')
        replaceReconCones = true;
        reconStartCones = 2;
        reconNewCones = cMosaic.LCONE_ID;
    elseif strcmp(reconChrom, 'chromTrit')
        replaceReconCones = true;
        reconStartCones = 3;
        reconNewCones = cMosaic.LCONE_ID;
    elseif strcmp(reconChrom, 'chromAllL')
        replaceReconCones = true;
        reconStartCones = [2 3];
        reconNewCones = cMosaic.LCONE_ID;
    elseif strcmp(reconChrom, 'chromAllM')
        replaceReconCones = true;
        reconStartCones = [1 3];
        reconNewCones = cMosaic.MCONE_ID;
    elseif strcmp(reconChrom, 'chromNorm')
        replaceReconCones = false;
        reconStartCones = []; 
        reconNewCones = [];
    else
        error('Unrecognized chromaticity input')
    end

    % Set render filennames
    if (forwardAORender)
        forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters, forwardSeedStr, forwardChrom);
    else
        forwardRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,forwardPupilDiamMM,forwardDefocusDiopters, forwardSeedStr, forwardChrom);
    end
    if (reconAORender)
        reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_AO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,reconPupilDiamMM,reconDefocusDiopters, reconSeedStr, reconChrom);
    else
        reconRenderStructureName = sprintf('%sDisplayRender_%d_%0.2f_%0.2f_%d_%d_NOAO_%0.2f_%s_%s.mat',displayName,fieldSizeMinutes,eccXDegs,eccYDegs,nPixels,reconPupilDiamMM,reconDefocusDiopters, reconSeedStr, reconChrom);
    end

    % Build render matrices/files or load from existing cache
    if (buildNewForward || ~exist(fullfile(aoReconDir, helpDir, forwardRenderStructureName),'file'))
        renderStructure = buildRenderStruct(aoReconDir, helpDir, eccXDegs, eccYDegs, ...
            fieldSizeDegs, nPixels, forwardPupilDiamMM, forwardAORender, forwardDefocusDiopters, ...
            overwriteDisplayGamma, displayName, displayFieldName, displayGammaBits, ...
            displayGammaGamma, forwardRandSeed, replaceForwardCones, forwardStartCones, forwardNewCones);
        save(fullfile(aoReconDir, helpDir,forwardRenderStructureName),'renderStructure');
        forwardRenderStructure = renderStructure; clear renderStructure;
    else
        clear forwardRenderStructure;
        load(fullfile(aoReconDir, helpDir, forwardRenderStructureName),'renderStructure');
        forwardRenderStructure = renderStructure; clear renderStructure; 
        grabRenderStruct(forwardRenderStructure, eccXDegs, eccYDegs, fieldSizeDegs, ...
            nPixels, forwardPupilDiamMM, forwardAORender, forwardDefocusDiopters)
    end

    % Set forward variables from loaded/built structure
    forwardRenderMatrix = forwardRenderStructure.renderMatrix;
    forwardConeMosaic = forwardRenderStructure.theConeMosaic;
    forwardOI = forwardConeMosaic.PSF;
    
    % Set display variable
    theForwardDisplay = forwardRenderStructure.theDisplay;
    
    % Clear forward render structure
    clear forwardRenderStructure;

    % Setup output directories
    outputMainName = sprintf('%s_%s_%0.2f_%0.2f_%d_%d_%s_%s%s', ...
        forwardAOStr,reconAOStr,forwardDefocusDiopters,reconDefocusDiopters,nPixels,fieldSizeMinutes,displayName,sparsePriorStr, versEditor);
    outputSubName = sprintf('%0.1f_%0.4f_%d_%0.2f_%0.2f_%0.2f_%0.2f_%s_%s_%s',60*stimSizeDegs, regPara,stride,stimBgVal,stimRVal,stimGVal,stimBVal, exciteSource, forwardChrom, reconChrom);
    outputDir = fullfile(aoReconDir,outputMainName,outputSubName);
    if (~exist(outputDir,'dir'))
        mkdir(outputDir);
    end
    
    % Show cone mosaic
    forwardConeMosaic.visualizeMosaic();
    
    % Generate an image stimulus
    %
    % Stimulus size in retinal degrees should not exceed 'fieldSizeDegs'
    stimSizeFraction = stimSizeDegs / fieldSizeDegs;
    if (stimSizeFraction > 1)
        error('Stimulus size too big given field size');
    end
    idxLB = round(nPixels * (0.5 - stimSizeFraction / 2));
    idxUB = round(nPixels * (0.5 + stimSizeFraction / 2));
    idxRange = idxLB:idxUB;
    
    % Set image pixels
    stimulusImageRGB = ones(nPixels, nPixels, 3) * stimBgVal;
    stimulusImageRGB(idxRange, idxRange, 1) = stimRVal;
    stimulusImageRGB(idxRange, idxRange, 2) = stimGVal;
    stimulusImageRGB(idxRange, idxRange, 3) = stimBVal;
    
    % Show the stimulus by creating an ISETBio scene
    meanLuminanceCdPerM2 = [];
    [stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
        meanLuminanceCdPerM2, forwardConeMosaic.Display);
    stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
    visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);

    % Compute excitations using forward render matrix.
    %
    % We would like these to agree with direct rendering, but
    % for reasons we are slowly coming to understand (quantization
    % is a non-linear process), they don't always.
    forwardExcitationsToStimulusRenderMatrix = forwardRenderMatrix*stimulusImageLinear(:);
    

    % Calculate the fundamentals
    % 
    % Grab the excitation values for each cone and calculate their mean as
    % quantal sensitivities
    sConeInd = [forwardConeMosaic.Mosaic.sConeIndices];
    sConeExcite = forwardExcitationsToStimulusRenderMatrix(sConeInd);
    sQuantalSensitivity = mean(sConeExcite);
    
    mConeInd = [forwardConeMosaic.Mosaic.mConeIndices];
    mConeExcite = forwardExcitationsToStimulusRenderMatrix(mConeInd);
    mQuantalSensitivty = mean(mConeExcite);
    
    lConeInd = [forwardConeMosaic.Mosaic.lConeIndices];
    lConeExcite = forwardExcitationsToStimulusRenderMatrix(lConeInd);
    lQuantalSensitivity = mean(lConeExcite);

    % Combine sensitivities into fundamental matrices and convert to energy
    coneQuantalFundamentals(:, cc) = [sQuantalSensitivity; ...
        mQuantalSensitivty; lQuantalSensitivity];
end


% Placeholder below for conversion of quanta to energy
% coneEnergyFundamentals = Quanta2Energy(coneQuantalFundamentals, 0);




