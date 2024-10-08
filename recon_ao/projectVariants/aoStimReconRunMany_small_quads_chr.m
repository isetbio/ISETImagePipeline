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
%     parpool(2);
% end

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
prBase.viewingDisplayName = 'conventional';
prBase.displayGammaBits = 12;
prBase.displayGammaGamma = 2;
displayScaleFactorList = [1];

%% Stimulus size
prBase.stimSizeDegsList = [1]/60; %[10 7.5 5.5 4.5 3.5 2 1] / 60;

% When we construct mosaics, add this much to the size of the stimulus
% area that we control, to account for effect of forward optical blur
% on the area the stimulus covers.  Just a little.
%
% This is handled by adding this number to the stimulus size in the call
% made by buildRenderStruct to setConeProportions, thus leaving the target
% size alone everywhere else. That is, this only affects the region of the
% central portion of the mosaic where we set the cone proportions.
prBase.forwardOpticalBlurStimSizeExpansionDegs = 1/60;

%% Mosaic information
%
% Mosaics are built to be matched to the stimulus size, generally.  Because
% building the render matrix is very slow, we may want to dissociate this a
% little at some point.  The first step in that direction is simply to
% define a separate parameter for the mosaic size specification and set
% it to be the same as the stimulus size specfication.  We will need to
% spend a little time on directory naming and parameter passing to actually
% dissociate these two things, so this next line should not be changed
% until that work is done.
prBase.mosaicStimSizeDegsList = prBase.stimSizeDegsList;

% At some point we should expose the parameter that controls the size of
% the immediate surround, because we may be interested in how varying that
% affects things. 

% The parameters here overide the default parameters
% for different regions, with those parametrs specified below.
% 
% Here, a list of regions are defined as the focal region, and we can cycle
% through L cone proportions and variants for that region.  Usually this is
% the center, but other choices are possible.
%
% For each focal region specified, the code cycles through the variants in
% the focalVariantList, replacing the corresponding value in
% regionVariantList on each iteration.
prBase.focalRegionList = ["center"];
prBase.focalVariantList = [2];
prBase.focalPropLList = [0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];

% Additional region variant params
%
% Set default variant and proportion L and S cones. Note that throughout
% the simulations, these values will hold and only one per group will be
% switched to the focal value 

% L proportions switch.  Sets cone proportions matched to region
% variant number.  This wasn't our original intent with how we
% would use the regionVariant flexibility, but for the near and
% far surrounds it is turning out to be useful to have a specific
% cone propL associated with each near and far surround region variant.
% That is because the bookkeeping doesn't separately allow us to easily
% track what happens for multiple instances of the same region variant
% number with different cone proportions.  I think this is probably
% OK.
prBase.regionVariant = [1 1 1];
prBase.propL = [0.0 0.0 0.0];
for rr = 2:3
    switch (prBase.regionVariant(rr))
        case {1, 6}
            prBase.propL(rr) = 0.67;
        case {2, 7}
            prBase.propL(rr) = 0;
        case {3, 8}
            prBase.propL(rr) = 1;
        case {4, 9}
            prBase.propL(rr) = 0.1;
        case {5, 10}
            prBase.propL(rr) = 0.9;
        otherwise
            error('Need to specify propL for this regionVariant case');
    end
end

% Set cone proportions for S for all regions.
% 
% Currently we are adjusting this by hand for different
% stim sizes, because we want to make sure we get some
% S cones in the central stimulated area even for the small
% sizes.  To do this, we need to boost the underlying
% proportion.  The conditional is handled in routine
% prFromBase.
prBase.propSLargeTarget = [0.1 0.07 0.07];  
prBase.propSSmallTarget = [0.15 0.07 0.07];
prBase.targetSizeSPropThresholdDegs = 6/60;

% Establish the size of the nearSurround annulus, input as a width value in
% arcmin. Nominal stimulus, pixel adjustents, and padding all encompass the
% center region, whereas this nearSurround begins just beyond that.
prBase.annulusWidthArc = 2; 

% Stimulus variant number
prBase.stimSeriesVariant = 1; 

%% Choose your journey
%
% Select what you would like to do, for efficiency's sake only recommend
% having one set to true at a time (reconstruct, renderMatrices, or mosaic
% montages)
buildRenderMatrix = true;
buildMosaicMontages = false;
runReconstructions = true;
summaryFigs = true;

%% Spatial parameters
%
% Common to forward and recon models
prBase.nPixels = 61;
prBase.trueCenter = round(prBase.nPixels/2);
if (rem(prBase.nPixels,2) ~= 1)
    error('Stimulus size logic requires odd number of pixels');
end

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

%% Force rewrite of existing render matrices
%
% The two buildNew flags here force a build of existing matrices, while
% if they are false and we are building, only ones that don't yet exist
% are built.  Turn these to true with trepidation because the same render
% matrix may get build over and over as one loops through conditions.  
% Probably what you actually want to do is go delete render matrices that
% you would like rebuilt.
if buildRenderMatrix
    buildNewForward = false;
    buildNewRecon = false;
end

%% Adjust desired visualization aspects of summary montages
figReconRows = false;
scaleToMax = true;
zoomToStim = true;
wavelengthUY = [580 570];

%% Mosaic cone parameters
%
% Determine if we're going to make a mosaic with cones specified and
% if we want to visualize region bounds. If equals false, default is base
% mosaic.
prBase.useCustomMosaic = true;
prBase.viewBounds = false;

% Wavelength specificaiton for calculations
prBase.wls = (400:1:700)';

%% Calculate the actual stimulus size given pixel quantization.
%
% It is these we want to use in constructing the mosaics and scenes.
prBase.pixelsPerMinute = prBase.nPixels/prBase.fieldSizeMinutes;
prBase.minutesPerPixel = 1/prBase.pixelsPerMinute;
prBase.availStimSizesPixels = (1:2:prBase.nPixels);
prBase.availStimSizesMinutes = prBase.minutesPerPixel*prBase.availStimSizesPixels;
prBase.availStimSizesDegs = prBase.availStimSizesMinutes/60;
for ss = 1:length(prBase.stimSizeDegsList)
    stimSizeDegsNominal = prBase.stimSizeDegsList(ss);
    [~,index]= min(abs(prBase.availStimSizesDegs - stimSizeDegsNominal));
    prBase.stimSizePixels(ss) = prBase.availStimSizesPixels(index);
    prBase.stimSizeDegs(ss) = prBase.availStimSizesDegs(index);
    stimSizeMinutes(ss) = prBase.stimSizeDegs(ss)*60;
    fprintf('Nominal stimulus %d: %0.3f minutes, actual used %0.3f minutes, %d pixels, minutes per pixel is %0.4f\n', ...
        ss,prBase.stimSizeDegsList(ss)*60,stimSizeMinutes(ss),prBase.stimSizePixels(ss),prBase.minutesPerPixel);
    if (rem(prBase.stimSizePixels(ss),2) ~= 1)
        error('We are assuming odd numer of pixels and stim pixel size');
    end
end

%% Stimulus color
%
% We can either specify an explicit list of RGB values, or generate an
% isoluminant series that varies between full green and full red.
% Which we do is controlled by the isoLumRG flag.  Specific parameters
% for each possibility are defined within the corresponding branch of the
% conditional just below.
isoLumRGAuto = false;
monoBgScale = 0.2;
monoGray = [0.4445; 0.5254; 0.5542];

if (~isoLumRGAuto)
    % These are rgb values (linear, before gamma correction)
    prBase.stimBgVal = monoGray * monoBgScale;
    stimrValList = [0.4615 0.3846 0.3077 0.2308 0.1538 0.1154 0.0769 0.0385 0.0000];
    stimgValList = [0.0081 0.0242 0.0403 0.0565 0.0726 0.0807 0.0888 0.0968 0.1049];
    stimbValList = [0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005];
else
    nEquiLumStimuli = 11;
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

% %% Position shifts
% %
% % Input desired x and y position for stimulus to be centered over. Function
% % will end if values exceed pixel limits.
% %
% % DHB: Currently commented out all the shifting code.  We could put this
% % back with a little thought, but I am inclined to remove completely to
% % simplify reading this whole routine which is already pretty complicated.
% %
% % Position specified in pixels, could consider specifying in minutes.
% % shiftInMinutesListX = [0];
% % shiftInMinutesListY = [0];
% 
% % Convert the shifts to pixel positions.
% %
% % shiftInPixelsListX = round(prBase.pixelsPerMinute*shiftInMinutesListX);
% % shiftInPixelsListY = round(prBase.pixelsPerMinute*shiftInMinutesListY);
% % if (shiftInPixelsListX ~= 0 | shiftInPixelsListY ~= 0)
% %     error('For right now assuming no shift in pixels specified');
% % end
% 
% % Consider including the term "-(prBase.nPixels * 0.1)" to better center
% % the quad stim at larger values so not pushing against the outer edge.
% %
% % DHB: I think prBase.stimCenter is now ignored in favor of the new
% % way of positioning the sitmulus implemented in aoStimRecon.  I don't
% % think it would be hard to put shifts back if we decide we need that.
% % centerXPosition = prBase.trueCenter + shiftInPixelsListX;
% % centerYPosition = prBase.trueCenter + shiftInPixelsListY;

%% Save center parameters as a vector.
prBase.stimCenter = [prBase.trueCenter ; prBase.trueCenter];

%% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional display.
prBase.sparsePriorStr = 'conventional';

%% Reconstruction parameters
%
% Should cycle through a few of these regs to optimize for 58x58 pixels
% Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
regParaList = [0.005]; %[0.1];
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
        for vv = 1:length(prBase.focalVariantList)
            for oo = 1:length(prBase.focalPropLList)
                for cc = 1:length(stimrValList)
                    for ff = 1:length(forwardDefocusDioptersList)
                        for rr = 1:length(regParaList)
                            for pp = 1:length(forwardPupilDiamListMM)
                                for dsf = 1:length(displayScaleFactorList)
                                    % These parameters do not affect mosaics or
                                    % render matrices.
                                    stimrVal(runIndex) = stimrValList(cc);
                                    stimgVal(runIndex) = stimgValList(cc);
                                    stimbVal(runIndex) = stimbValList(cc);
                                    stimCenter(:,runIndex) = prBase.stimCenter; % deltaCenterList(:,yy);
                                    regPara(runIndex) = regParaList(rr);
                                    displayScaleFactor(runIndex) = displayScaleFactorList(dsf);

                                    % These do affect mosaics because we
                                    % design mosaics to have desired properties
                                    % within the stimulus region and regions
                                    % adjacent to it. This is taken into account
                                    % when we build montages of mosaics.
                                    stimSizeDegs(runIndex) = prBase.stimSizeDegs(ss);
                                    stimSizePixels(runIndex) = prBase.stimSizePixels(ss);
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

% Remvoew the final run index bump to match lengths
runIndex = runIndex - 1;

%% Render Structures
%
% Build the render structures and note that being done in bulk so only run
% over length 1 since builders use the full value lists. This is not the
% case if the core render matrix params above are altered, then must cycle
% through. (Suspect this might still have overcalculation/redundancy but
% explore more later)
if buildRenderMatrix
    for pp = 1:runIndex
        % Set up paramters structure for this loop, filling in fields that come
        % out of lists precreated above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimSizePixels,stimrVal,stimgVal,stimbVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
            focalRegion, focalPropL, focalVariant);
        cnv = computeConvenienceParams(pr);
        
        % Build foward cone mosaic and render matrix if needed
        if (buildNewForward || ~exist(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'file'))
            renderStructure = buildRenderStruct(pr, cnv, "forward");
            clear renderStructure;
        end
        
        % Build recon cone mosaic and render structure if needed
        if (buildNewRecon || ~exist(fullfile(cnv.reconRenderDirFull, cnv.renderName),'file'))
            renderStructure = buildRenderStruct(pr, cnv, "recon");
            clear renderStructure;
        end
    end
end

%% Visualize mosaics
%
% Building mosaics themselves is fast, and we control the random number
% generator. Follows the same format as for the render matrices.
if buildMosaicMontages
    for pp = 1
        % Set up paramters structure for this loop, filling in fields that come
        % out of lists precreated above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimSizePixels,stimrVal,stimgVal,stimbVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
            focalRegion, focalPropL, focalVariant);
        cnv = computeConvenienceParams(pr);
        buildMosaicMontage(pr, cnv, "forward");
        buildMosaicMontage(pr, cnv, "recon");
    end
end

%% Run aoStimRecon.m
%
% Helpful tip: If you get weird error crashes that don't seem to make
% sense, try turning the parfor loop below into a normal loop because the
% issue is likley a masked issue within the aoStimRecon file. Just remember
% to turn it back to a parfor loop when finished.
if runReconstructions
    for pp = 1:runIndex        
        % Set up parameters structure for this loop, filling in fields that come
        % out of lists above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimSizePixels,stimrVal,stimgVal,stimbVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
            focalRegion, focalPropL, focalVariant);
        
        % Compute convenience parameters
        cnv = computeConvenienceParams(pr);
        
        % Call the driving function
        aoStimRecon(pr,cnv);
    end
end

%% Build Summary Figs
%
% Integrate the aoStimReconRerunFigs script into this one for a centralized
% region of post processing. Set it up as another option.
%
% Some of this is done rather precariously so care should be taken as the
% project progresses to ensure the things being cycled over are actually
% what we want.
if actuasummaryFigs
    
    % Bookkeeping variables for number of stimuli and propL as dimensions
    % of future plots
    numStim = length(stimrValList);
    numProp = length(prBase.focalPropLList);
    
    % Identify the main variables we're concerned about for these
    % simulations for ease of computation, except for the final output
    % level (i.e. dont include stimrVal despite the fact that we also
    % change color since those directories contain the xRunOutput.m file).
    mainVars = [stimSizeDegs; focalRegion; focalVariant; focalPropL];
    [~, mainVarsInd] = unique(mainVars.', 'rows', 'stable');
    
    % Cycle only over the instances where the main variables change.
    [~, varInd1] = unique(mainVars(1,:));
    varInd1 = [varInd1' (length(mainVars)+1)];

    % Variable One: Stim Size
    for vo = 1:length(varInd1)-1
        holderVars1 = mainVars(: , varInd1(vo):varInd1(vo+1)-1);
        [~, varInd2] = unique(mainVars(2,:));
        varInd2 = [varInd2' (length(holderVars1)+1)];
        varInd2 = [varInd1(vo) - 1 + varInd2];
        
        % Variable Two: Focal Region
        for vt = 1:length(varInd2)-1
            holderVars2 = mainVars(: , varInd2(vt):varInd2(vt+1)-1);
            [~, varInd3] = unique(mainVars(3,:));
            varInd3 = [varInd3' (length(holderVars2)+1)];
            varInd3 = [varInd2(vt) - 1 + varInd3];
            
            % Variable Three: Focal Variant
            for ve = 1:length(varInd3)-1
                holderVars3 = mainVars(: , varInd3(ve):varInd3(ve+1)-1);
                [~, varInd4] = unique(mainVars(4,:));
                varInd4 = [varInd3(ve) - 1 + varInd4'];
                
                fullReconSummary = [];
                
                % Variable Four: Focal Prop L
                for vf = 1:length(varInd4)
                    % Readjust the index value according to the levels that are
                    % actually pertinent.
                    newInd = varInd4(vf);
                    
                    % Set up paramters structure for this loop, filling in fields that come
                    % out of lists above.
                    pr = prFromBase(prBase,newInd,stimSizeDegs,stimSizePixels,stimrVal,stimgVal,stimbVal, ...
                        stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
                        forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor, ...
                        focalRegion, focalPropL, focalVariant);
                    
                    % Compute convenience parameters
                    cnv = computeConvenienceParams(pr);

                    % Patch propInfoFile.m to create the excel sheets with
                    % mosaic information. This function will slow down code
                    % noticeably since loading in each render structure to
                    % get pertinent information (might want to rethink this
                    % so can save the info at render creation). Should only
                    % need to be run once for each render structure though,
                    % after which it will see the file cached and ignore. 
                    propInfoFile(pr,cnv);
        
                    % Call the function to build the summary plots.
                    [stimSummary, reconSummary] = grabImageInfo(pr, cnv, numStim, ...
                        "figReconRows", figReconRows, "scaleToMax", scaleToMax, ...
                        "wls", pr.wls, "zoomToStim", zoomToStim);
                    
                    % Store the collected info in a running cell and utilize when actually
                    % building the full summary figures.
                    fullReconSummary = [fullReconSummary; reconSummary];
                    
                end
                
                % Make and save summary figures
                buildSummaryFigs(pr, cnv, numStim, numProp, ...
                    fullReconSummary, stimSummary, 'wavelengthUY', wavelengthUY);

                % Close up summary figures so they don't pile up
                close all
            end
        end
    end
end

% DHB: This code was crashing so I commented out the whole section.
%
% %% Appendix 1
% %
% % This spot allows for quick visualization of a single mosaic without
% % having to produce an entire mosaic montage. Manually load in the desired
% % render structure from the proper directory. Then uncomment the needed
% % portion below and run in the command window. Be sure to comment it out
% % when finished.
% 
% % % LOAD THE RENDER STRUCTURE FIRST
% % renderStructure.theConeMosaic.visualizeMosaic()
% 
% tempFig = figure; clf;
% forwardConeMosaic.Mosaic.visualize(...
%     'figureHandle', tempFig, ...
%     'axesHandle', [], ...
%     'withSuperimposedOpticalImage', forwardOI, ...
%     'outlinedConesWithIndices', pr.kConeIndices, ...
%     'plotTitle','Forward OI on Forward Mosaic','superimposedOIAlpha',0.7);
% % 
% % Establish the new center point based on comparison with where the OI
% % actually lands. Also include scaling based on the OI
% % pr.eccXDegs = 1.9945; 
% % pr.eccYDegs = 0.0055;
% scaleForOI = 1; %1.23;
% % scales for [2 3.5 10] = [1.55 1.23 1.05]
% 
% hold on; 
% plot(pr.eccXDegs, pr.eccYDegs, '.', 'MarkerSize', 40);
% 
% xBounds = [];
% yBounds = [];
% xBounds(1,:) = [pr.eccXDegs pr.eccXDegs];
% yBounds(1,:) = [pr.eccYDegs pr.eccYDegs];
% centerWidth = scaleForOI * (pr.stimSizeDegs/2);
% xBounds(2,:) = [xBounds(1)-centerWidth, xBounds(2)+centerWidth];
% yBounds(2,:) = [yBounds(1)-centerWidth, yBounds(2)+centerWidth];
% 
% % Outline the stimulus region
% hold on
% rectangle('Position', [xBounds(2,1) yBounds(2,1) ...
%     (xBounds(2,2) - xBounds(2,1)) ...
%     (yBounds(2,2) - yBounds(2,1))], 'LineWidth', 3)
% 
% % 
% % % Scale Factor
% % a = 1.639;
% % b = -1.598;
% % c = 1.009;
% % 
% % % Arcmin Anulus
% % % Data points represent number of arcmin added to the diameter of the
% % % stimulus
% % % a = 1.44;
% % % b = -0.7641;
% % % c = 0.2521;
% % 
% % % Power Fits (2-term)
% % val = a*x^b + c;
% % 
% % 
% % 
% % % Visaulize index numbers over each cone in the mosaic
% % for i=1:length(renderStructure.theConeMosaic.Mosaic.coneTypes)
% %     txt = int2str(i);
% %     t = text((renderStructure.theConeMosaic.Mosaic.coneRFpositionsDegs(i,1)-0.005), ...
% %         renderStructure.theConeMosaic.Mosaic.coneRFpositionsDegs(i,2),txt);
% %     t.FontSize=11;
% %     hold on
% % end
% % 
% % % Call function to manually select and swap cones (see mouseClick)
% % g = @(x, y) mouseClick(x, y, renderStructure.theConeMosaic);
% % set(gcf, 'WindowButtonDownFcn', g)
% % keyboard
% 
% 
% %% Appendix 2
% % 
% % Visualize the OI on the retinal mosaic, to be conducted after running the
% % aoStimRecon portion outside of the parfor loop and inserting a line break
% % somewhere near line 520. Use some easily identifiable and removable stim
% % color like 0.5 0.5 0.5 for a maxiteration of 2 to quickly get through the
% % initial phases. 
% 
% 
% tempFig = figure; clf;
% forwardConeMosaic.Mosaic.visualize(...
%     'figureHandle', tempFig, ...
%     'axesHandle', [], ...
%     'withSuperimposedOpticalImage', forwardOI, ...
%     'outlinedConesWithIndices', pr.kConeIndices, ...
%     'plotTitle','Forward OI on Forward Mosaic','superimposedOIAlpha',0.7);


