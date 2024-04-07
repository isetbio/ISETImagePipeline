%% aoStimReconMany
%
% Description:
%    Call aoStimRecon with many combinations of parameters.
%
% See also: aoStimRecon

% History:
%   08/15/22  dhb  Wrote after converting aoStimRecon to a function
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options
%   09/22/22  chr  Convert to its own dichrom file
%   09/27/22  chr  Incorporate inputs for stimulus centering position
%   10/05/22  dhb  Lots of changes for parallel
%   02/20/23  chr  Updates based on Tuten Meeting

%% Clear
clear; close all;

%% Control size of parpool, otherwise may crush memory
% thePool = gcp('nocreate');
% if (isempty(thePool))
%     parpool(5);
% end

%% Portion for ReRunning Figures from previous simulations after updates
rrf = struct;
rrf.rerunImages = false;
rrf.startDisplayName = 'mono';
rrf.viewingDisplayName = 'conventional';
rrf.stimDispScale = 1;
rrf.reconDispScale = 1;

%% Set defaults in prBase
prBase = prBaseDefaults;

%% Version editor string
%
% Helps us keep different calcs separate
prBase.versEditor = 'small_quads_chr';

%% Parameters
%
% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
prBase.displayName = 'mono';
prBase.displayGammaBits = 12;
prBase.displayGammaGamma = 2;
displayScaleFactorList = [1];

%% Spatial parameters
%
% Common to forward and recon models
prBase.nPixels = 50;
prBase.trueCenter = round(prBase.nPixels/2);

%% Mosaic general parameters
prBase.fieldSizeMinutes = 30;
prBase.eccXDegs = 2.0;
prBase.eccYDegs = 0.0;
prBase.forwardRandSeed = false;
prBase.reconRandSeed = false;
prBase.forwardEccVars = false;
prBase.reconEccVars = false;
prBase.reconstructfromRenderMatrix = true;
prBase.addPoissonNoise = false;

%% Choose your journey
%
% Select what you would like to do, for efficiency's sake only recommend
% having one set to true at a time (reconstruct, renderMatrices, or mosaic
% montages)
runReconstructions = true;
buildRenderMatrix = true;
buildMosaicMontages = false;

% The two buildNew flags here force a build of existing matrices, while
% if they are false and we are building, only ones that don't yet exist
% are built.
if buildRenderMatrix
    buildNewForward = false;
    buildNewRecon = false;
end

%% Mosaic cone parameters
%
% Determine if we're going to make a fancy mosaic with cones specified and
% if we want to visualize region bounds. If equals false, default is base
% mosaic. Also put wls calc at top level.
prBase.useCustomMosaic = true;
prBase.viewBounds = false;
prBase.wls = (400:1:700)';

% These are the specific values taken in by the AO script, for this project
% want it to be relatively limited for the sake of speed.
prBase.stimSizeDegsList = [10] / 60;
prBase.focalRegionList = ["center"];
prBase.focalPropLList = [0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];
prBase.focalVariantList = 1;

% Set default variant and proportion L and S cones. Note that throughout
% the simulations, these values will hold and only one per group will be
% switched to the focal value (i.e., if focal variant is 2
prBase.regionVariant = [1 1 1];
prBase.propL = [0.5 0.5 0.5];
prBase.propS = [0.10 0.10 0.10];
%% Stimulus color
%
% We can either specify an explicit list of RGB values, or generate an
% isoluminant series that varies between full green and full red.
% Which we do is controlled by the isoLumRG flag.  Specific parameters
% for each possibility are defined within the corresponding branch of the
% conditional just below.
isoLumRG = true;
if (~isoLumRG)
    % These are rgb values (linear, before gamma correction)
    prBase.stimBgVal = 0.3;
    stimrValList = 0.80;
    stimgValList = 0.65;
    stimbValList = 0.10;
else
    nEquiLumStimuli = 11;
    monoBgScale = 0.5;
    switch (prBase.displayName)
        case 'conventional'
            displayFieldName = 'CRT12BitDisplay';
            overwriteDisplayGamma = false;
            prBase.stimBgVal = [0.1 0.1 0.1];
        case 'mono'
            displayFieldName = 'monoDisplay';
            overwriteDisplayGamma = true;

            % The following values are the approximations for a 0.5 uniform
            % field on a conventional display using the tutorial. These can
            % be scaled accordingly for a lighter or darker background on
            % the interval 0-1/max(monoGray). Temporary patch, once
            % everything else is cleaned may consider replacing with a call
            % to RGBRenderAcrossDisplays
            monoGray = [0.444485445018075; ...
                0.525384925401570; ...
                0.554173733733909];

            % DON'T FORGET TO PUT IN AN ACTUAL VALUE FOR THE "1" HERE
            prBase.stimBgVal = monoGray * monoBgScale;
        otherwise
            error('Unknown display specified');
    end

    % Load the appropriate display
    theDisplayLoad = load(fullfile(prBase.aoReconDir, 'displays', [prBase.displayName 'Display.mat']));
    eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);

    if (overwriteDisplayGamma)
        gammaInput = linspace(0,1,2^prBase.displayGammaBits)';
        gammaOutput = gammaInput.^prBase.displayGammaGamma;
        theDisplay.gamma = gammaOutput(:,[1 1 1]);
    end

    % Using the 1 nm sampling to agree w/ tutorial
    theDisplay = displaySet(theDisplay, 'wave', prBase.wls);

    % Get information we need to render scenes from their spectra through
    % the recon display.
    theXYZStruct = load('T_xyz1931');
    T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,prBase.wls);
    M_rgbToXYZ = T_XYZ*displayGet(theDisplay,'spd primaries')*(prBase.wls(2)-prBase.wls(1));
    M_XYZTorgb = inv(M_rgbToXYZ);
    rPrimaryLuminance = M_rgbToXYZ(2,1);
    gPrimaryLuminance = M_rgbToXYZ(2,2);

    % Compute as set of equally spaced r/(r+g) values that lead
    % to equal luminance stimuli.
    thePrimaries = displayGet(theDisplay,'spd primaries');
    rOverRPlusG = linspace(1,0,nEquiLumStimuli);
    gOverRPlusG = 1-rOverRPlusG;
    gPrimaryAdjust = rPrimaryLuminance/gPrimaryLuminance;
    rRaw = rOverRPlusG;
    gRaw = gOverRPlusG;
    gAdjust = gRaw*rPrimaryLuminance/gPrimaryLuminance;
    b = 0.001;
    equiLumrgbValues = [rRaw ; gAdjust; b*ones(size(rRaw))];
    inputLinearrgbValues = 0.5*equiLumrgbValues/max(equiLumrgbValues(:));

    stimrValList = inputLinearrgbValues(1,:);
    stimgValList = inputLinearrgbValues(2,:);
    stimbValList = inputLinearrgbValues(3,:);
end

% Check that all channels receive same number of inputs
if (length(stimgValList) ~= length(stimrValList) || length(stimbValList) ~= length(stimrValList))
    error('Stimulus value lists must have same length');
end

%% Positions
%
% Input desired x and y position for stimulus to be centered over. Function
% will end if values exceed pixel limits.
%
% Position specified in pixels, could consider specifying in minutes.
pixelsPerMinute = prBase.nPixels/prBase.fieldSizeMinutes;
shiftInMinutesListX = [0];
shiftInMinutesListY = [0];
fullSquareShift = false;

% Convert the shifts to pixel positions
shiftInPixelsListX = round(pixelsPerMinute*shiftInMinutesListX);
shiftInPixelsListY = round(pixelsPerMinute*shiftInMinutesListY);

% Consider including the term "-(prBase.nPixels * 0.1)" to better center
% the quad stim at larger values so not pushing against the outer edge.
centerXPosition = prBase.trueCenter + shiftInPixelsListX;
centerYPosition = prBase.trueCenter + shiftInPixelsListY;
prBase.stimCenter = [centerXPosition ; centerYPosition];

% Loop through created pixel positions if want to create a square grid of
% movement instead of default horizontal shift
if (fullSquareShift)
    centerXPosition = repelem(prBase.stimCenter(1,:), length(shiftInMinutesList));
    centerYPosition = repmat(prBase.stimCenter(1,:), [1,length(shiftInMinutesList)]);
    prBase.stimCenter = [centerXPosition; centerYPosition];
end
deltaCenterList = [prBase.stimCenter - prBase.trueCenter];

%% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional display.
prBase.sparsePriorStr = 'conventional';

%% Reconstruction parameters
%
% Should cycle through a few of these regs to optimize for 58x58 pixels
% Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
regParaList = 0.1;
prBase.stride = 2;
prBase.maxReconIterations = 2000;
prBase.whiteNoiseStarts = 0;
prBase.pinkNoiseStarts = 1;
prBase.sparsePriorPatchStarts = 0;
prBase.stimulusStart = false;
% Should note, the line below only works if have >= 3 starting fields
prBase.uniformStartVals = []; % [[0.5 0.5 0.5]' [0 0.5 0]' [0 0 0.5]'];
prBase.boundedSearch = false;

% Use AO in forward rendering? And determine optics pupil size
prBase.forwardAORender = true;
prBase.forwardNoLCA = true;
prBase.reconAORender = false;
prBase.reconNoLCA = false;
reconPupilDiamListMM =  [2];
forwardPupilDiamListMM = [7];

% Define optics.  Subject only matters if we use a database.  Ignored for
% Marimont and Wandell.  For database, subjectID of 0 means diffraction
% limited.
%
% Databases are 'Polans2015' and 'Artal2012'
prBase.forwardSubjectID = 6;
prBase.forwardZernikeDataBase = 'Polans2015';
prBase.reconSubjectID = 6;
prBase.reconZernikeDataBase = 'Polans2015';

% Residual defocus for forward and recon rendering, of equal sizes
forwardDefocusDioptersList = [0.05];
reconDefocusDioptersList = [0.0];

%% Set up list conditions
runIndex = 1;
for ss = 1:length(prBase.stimSizeDegsList)
    for gg = 1:length(prBase.focalRegionList)
        for oo = 1:length(prBase.focalPropLList)
            for vv = 1:length(prBase.focalVariantList)
                for cc = 1:length(stimrValList)
                    for yy = 1:size(deltaCenterList,2)
                        for ff = 1:length(forwardDefocusDioptersList)
                            for rr = 1:length(regParaList)
                                for pp = 1:length(forwardPupilDiamListMM)
                                    for dsf = 1:length(displayScaleFactorList)
                                        % These parameters do not affect mosaics or
                                        % render matrices.
                                        stimrVal(runIndex) = stimrValList(cc);
                                        stimgVal(runIndex) = stimgValList(cc);
                                        stimbVal(runIndex) = stimbValList(cc);
                                        stimCenter(:,runIndex) = deltaCenterList(:,yy);
                                        regPara(runIndex) = regParaList(rr);
                                        displayScaleFactor(runIndex) = displayScaleFactorList(dsf);

                                        % These do affect mosaics because we
                                        % design mosaics to have desired properties
                                        % within the stimulus region and regions
                                        % adjacent to it. This is taken into account
                                        % when we build montages of mosaics.
                                        stimSizeDegs(runIndex) = prBase.stimSizeDegsList(ss);
                                        focalRegion(runIndex) = prBase.focalRegionList(gg);
                                        focalPropL(runIndex) = prBase.focalPropLList(oo);
                                        focalVariant(runIndex) = prBase.focalVariantList(vv);

                                        % These parameters do by their nature directly affect either the mosaic
                                        % or the render matrices beyond the scope of our current project
                                        forwardDefocusDiopters(runIndex) = forwardDefocusDioptersList(ff);
                                        reconDefocusDiopters(runIndex) = reconDefocusDioptersList(ff);
                                        forwardPupilDiamMM(runIndex) = forwardPupilDiamListMM(pp);
                                        reconPupilDiamMM(runIndex) = reconPupilDiamListMM(pp);

                                        % Bump condition index
                                        runIndex = runIndex + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Set the multeRenderMatrixParams flag properly based on condition lists
%
% Most of the time, we only need to build one set of reconstruction
% mosaics, because many of the things we could vary are not in fact varied
% across different reconstruction conditions.  For example, we typically
% use the same regularization parameter across trials. Setting this flag
% to false prevents us from making the same mosaics and render matrices
% over and over and over again.
%
% Note however, that our build routine builds entire sets of mosaics and
% render matrices that take the stimulus size into account, so we can still
% keep this flag as false even when that list has length greater than 1.
multRenderMatrixParams = false;
if (    length(forwardDefocusDioptersList) == 1 & ...
        length(reconDefocusDioptersList) == 1 & ...
        length(forwardPupilDiamListMM) == 1 & ...
        length(reconPupilDiamListMM) == 1 )
    multRenderMatrixParams = false;
else
    multRenderMatrixParams = true;
end

%% Render Structures
%
% Build the render structures and note that being done in bulk so only run
% over length 1 since builders use the full value lists. This is not the
% case if the core render matrix params above are altered, then must cycle
% through. (Suspect this might still have overcalculation/redundancy but 
% explore more later)
if buildRenderMatrix
    if multRenderMatrixParams
        for pp = 1:length(regPara)
            % Set up paramters structure for this loop, filling in fields that come
            % out of lists precreated above.
            pr = prFromBase(prBase,pp,stimSizeDegs,stimrVal,stimgVal,stimbVal, ...
                stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
                forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
                focalRegion, focalPropL, focalVariant);
            cnv = computeConvenienceParams(pr);

            % Build foward cone mosaic and render matrix if needed
            if (buildNewForward || ~exist(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'file'))
                renderStructure = buildRenderStruct(pr, cnv, "forward");
            end

            % Build recon cone mosaic and render structure if needed
            if (buildNewRecon || ~exist(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'file'))
                renderStructure = buildRenderStruct(pr, cnv, "recon");
            end
        end
    else
        for pp = 1
            % Set up paramters structure for this loop, filling in fields that come
            % out of lists precreated above.
            pr = prFromBase(prBase,pp,stimSizeDegs,stimrVal,stimgVal,stimbVal, ...
                stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
                forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
                focalRegion, focalPropL, focalVariant);
            cnv = computeConvenienceParams(pr);

            % Build foward cone mosaic and render matrix if needed
            if (buildNewForward || ~exist(fullfile(cnv.reconRenderDirFull, cnv.renderName),'file'))
                renderStructure = buildRenderStruct(pr, cnv, "forward");
            end

            % Build recon cone mosaic and render structure if needed
            if (buildNewRecon || ~exist(fullfile(cnv.reconRenderDirFull, cnv.renderName),'file'))
                renderStructure = buildRenderStruct(pr, cnv, "recon");
            end
        end
    end
end

%% Visualize mosaics
%
% Building mosaics themselves is fast, and we control the random number
% generator. Follows the same format as for the render matrices. 
if buildMosaicMontages
    if multRenderMatrixParams
        for pp = 1:length(regPara)
            % Set up paramters structure for this loop, filling in fields that come
            % out of lists precreated above.
            pr = prFromBase(prBase,pp,stimSizeDegs,stimrVal,stimgVal,stimbVal, ...
                stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
                forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
                focalRegion, focalPropL, focalVariant);
            %         pr.quadSelect = quadSelect(:,pp);
            cnv = computeConvenienceParams(pr);
            buildMosaicMontage(pr, cnv, "forward");
            buildMosaicMontage(pr, cnv, "recon");
        end
    else
        for pp = 1
            % Set up paramters structure for this loop, filling in fields that come
            % out of lists precreated above.
            pr = prFromBase(prBase,pp,stimSizeDegs,stimrVal,stimgVal,stimbVal, ...
                stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
                forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
                focalRegion, focalPropL, focalVariant);
            %         pr.quadSelect = quadSelect(:,pp);
            cnv = computeConvenienceParams(pr);
            buildMosaicMontage(pr, cnv, "forward");
            buildMosaicMontage(pr, cnv, "recon");
        end
    end
end

%% Run aoStimRecon.m
% THIS SHOULD BE A PARFOR AFTERWARDS DON'T FORGET
if runReconstructions
    parfor pp = 1:length(regPara)

        % Set up paramters structure for this loop, filling in fields that come
        % out of lists above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimrVal,stimgVal,stimbVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM, ...
            displayScaleFactor, focalRegion, focalPropL, focalVariant);
        pr.quadSelect = quadSelect(:,pp);

        % Compute convenience parameters
        cnv = computeConvenienceParams(pr);

        % Call the driving function
        aoStimRecon(pr,cnv, rrf);
    end
end