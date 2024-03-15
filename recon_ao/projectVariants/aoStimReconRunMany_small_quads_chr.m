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

% To keep account for the different computers simulations might be run on
% simulatneously. macPro = trashcan, iMac = Semin's old computer
prBase.system = 'macPro';

%% Calculate Cone Proportions


% Provide input matrix with desired surrounding annulus widths. Each column
% corresponds to another annulus moving outward. Each row corresponds to a
% different series of annuli layers (labeled as "conditions")
% 
% 
% For reference, 2-3 arcmin seems to agree with OI spread. 
prBase.annWidthArc = [1; 2];

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
runReconstructions = false;
buildRenderMatrix = false;
buildMosaicMontages = true; 

if buildRenderMatrix
buildNewForward = false;
buildNewRecon = false;
end
%% Mosaic cone domain
% Top level domain values of all possible combinations we'll want to
% run. Useful for rapidly building render matrices or viewing mosaic
% montages, but is not sent into the aoScript to avoid overrunning. 
prBase.viewMosaicMontage = false; 
prBase.setProps = true; 
prBase.viewBounds = false; 

prBase.focalVariantDomain = 1; %1:5;
prBase.stimSizeDegsDomain = 10/60;%[2 3.5 10] / 60;
prBase.focalRegionDomain = "center"; %["center" "nearSurround" "distantSurround"];
prBase.focalPropLListDomain = 0.9;%[0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];

%% Mosaic cone parameters

% These are the specific values taken in by the AO script, for this project
% want it to be relatively limited for the sake of speed. 
prBase.focalVariant = 1;
stimSizeDegsList = [3.5] / 60;
prBase.focalRegion = "center"; % Options: center, nearSurround, distantSurround, multiple, global
prBase.focalPropLList = [0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];

%% Stimulus color
% rgb values (before gamma correction)
prBase.stimBgVal = 0.3;
stimRValList = 0.80;
stimGValList = 0.65;
stimBValList = 0.10;

isoLumRG = true;
nEquiLumStimuli = 11;

if (isoLumRG)
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
            prBase.stimBgVal = monoGray * 1;

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
    wls = (400:1:700)';
    theDisplay = displaySet(theDisplay, 'wave', wls);

    % Get information we need to render scenes from their spectra through
    % the recon display.
    theXYZStruct = load('T_xyz1931');
    T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
    M_rgbToXYZ = T_XYZ*displayGet(theDisplay,'spd primaries')*(wls(2)-wls(1));
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

    stimRValList = inputLinearrgbValues(1,:);
    stimGValList = inputLinearrgbValues(2,:);
    stimBValList = inputLinearrgbValues(3,:);
end

% Check that all channels receive same number of inputs
if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
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
quadCenters = round((prBase.nPixels) / 4);
centerXPosition = prBase.trueCenter + quadCenters + shiftInPixelsListX;
centerYPosition = prBase.trueCenter + quadCenters + shiftInPixelsListY;
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
for ss = 1:length(stimSizeDegsList)
    for cc = 1:length(stimRValList)
        for yy = 1:size(deltaCenterList,2)
            for ff = 1:length(forwardDefocusDioptersList)
                for rr = 1:length(regParaList)
                    for dd = 1:length(forwardChromList)%%%%%%%%%%%%%
                        for pp = 1:length(forwardPupilDiamListMM)
                            for qq = 1:length(quadSelectList(1,:))
                                for dsf = 1:length(displayScaleFactorList)

                                    stimSizeDegs(runIndex) = stimSizeDegsList(ss);

                                    stimRVal(runIndex) = stimRValList(cc);
                                    stimGVal(runIndex) = stimGValList(cc);
                                    stimBVal(runIndex) = stimBValList(cc);

                                    stimCenter(:,runIndex) = deltaCenterList(:,yy);

                                    forwardDefocusDiopters(runIndex) = forwardDefocusDioptersList(ff);
                                    reconDefocusDiopters(runIndex) = reconDefocusDioptersList(ff);

                                    regPara(runIndex) = regParaList(rr);

                                    forwardChrom(runIndex) = forwardChromList(dd);
                                    reconChrom(runIndex) = reconChromList(dd);

                                    forwardPupilDiamMM(runIndex) = forwardPupilDiamListMM(pp);
                                    reconPupilDiamMM(runIndex) = reconPupilDiamListMM (pp);

                                    quadSelect(:,runIndex) = quadSelectList(:,qq);

                                    displayScaleFactor(runIndex) = displayScaleFactorList(dsf);

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

%% Render Structures
%
% Either build the render structure OR visualize a montage of possible
% mosaics OR load an existing render struct. NOTE: All paths are mutually
% exclusive, can only do one of the above at a time. Loading is done within
% the aoStimRecon script in the chunk below.
for pp = 1:length(regPara)
    % Set up paramters structure for this loop, filling in fields that come
    % out of lists precreated above.
    pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
        stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
        forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor);
    pr.quadSelect = quadSelect(:,pp);
    cnv = computeConvenienceParams(pr);

    % Build foward cone mosaic and render matrix if needed
    if (buildNewForward || prBase.viewMosaicMontage || ...
            ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'file'))
        %             renderStructure = buildrenderStruct(pr, cnv);
        renderStructure = buildRenderStruct(pr, cnv, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardNoLCA, ...
            pr.forwardDefocusDiopters, pr.forwardRandSeed, cnv.replaceForwardCones, cnv.forwardStartCones, ...
            cnv.forwardNewCones, pr.forwardEccVars, pr.forwardSubjectID, pr.forwardZernikeDataBase, pr.forwardChrom);
        %             save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'renderStructure','-v7.3');
        %             forwardRenderStructure = renderStructure; clear renderStructure;
    end

    % Build recon cone mosaic and render structure if needed
    if (buildNewRecon || prBase.viewMosaicMontage || ...
            ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr, cnv, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconNoLCA, ...
            pr.reconDefocusDiopters, pr.reconRandSeed, cnv.replaceReconCones, cnv.reconStartCones, ...
            cnv.reconNewCones, pr.reconEccVars, pr.reconSubjectID, pr.reconZernikeDataBase, pr.reconChrom);
        %             save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'renderStructure','-v7.3');
        %             reconRenderStructure = renderStructure; clear renderStructure;
    end
end



keyboard();

%% Run aoStimRecon.m
% THIS SHOULD BE A PARFOR AFTERWARDS DON'T FORGET
if reconstruct
    parfor pp = 1:length(regPara)

        % Set up paramters structure for this loop, filling in fields that come
        % out of lists above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor);
        pr.quadSelect = quadSelect(:,pp);

        % Compute convenience parameters
        cnv = computeConvenienceParams(pr);

        % Call the driving function
        aoStimRecon(pr,cnv, rrf);
    end
end

