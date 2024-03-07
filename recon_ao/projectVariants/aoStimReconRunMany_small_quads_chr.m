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
prBase.system = 'chrTrashcan';

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

%% Mosaic parameters
prBase.fieldSizeMinutes = 30;
prBase.eccXDegs = 2.0;
prBase.eccYDegs = 0.0;
prBase.forwardRandSeed = false;
prBase.reconRandSeed = false;
prBase.forwardEccVars = false;
prBase.reconEccVars = false;
prBase.reconstructfromRenderMatrix = true;
prBase.addPoissonNoise = false;

% Mosaic chromatic type, options are:
%    "chromNorm", "chromProt", "chromDeut", "chromTrit",
%    "chromAllL", "chromAllM", "chromAllS", "quadSeq" and number -142
% ["quadSeq132" "quadSeq137" "quadSeq138" "quadSeq139" "quadSeq140" "quadSeq141" "quadSeq142"];
forwardChromList = "quadSeqNew"; %["quadSeq128" "quadSeq129" "quadSeq130" "quadSeq131" "quadSeq132" "quadSeq133" "quadSeq134" "quadSeq135" "quadSeq136"]; 
reconChromList   = "quadSeqNew"; %["quadSeq128" "quadSeq129" "quadSeq130" "quadSeq131" "quadSeq132" "quadSeq133" "quadSeq134" "quadSeq135" "quadSeq136"];  

% Build new sequence by
prBase.quads(1).name  = 'useQuadSeq';
prBase.quads(1).value = true;

% If want to apply percentages to the full mosaic instead of a
% quadrant, set to true and use percentages corresponding to Quadrant 4
% (prBase.quads(5))
fullMosaicPercent = true;

if(prBase.quads(1).value)
    % Initialize storage structure with information on each quadrant
    prBase.quads(2).name = 'Quad1';
    prBase.quads(3).name = 'Quad2';
    prBase.quads(4).name = 'Quad3';
    prBase.quads(5).name = 'Quad4';

    % Enter desired percent as decimal of L cones per region across
    % quadrants. The remaining percent will be made of M cones. Entries
    % should start with outermost regions first and progress inward
    prBase.quads(2).percentL = [0.27];
    prBase.quads(3).percentL = [0.53];
    prBase.quads(4).percentL = [0.71];
    prBase.quads(5).percentL = [0.50];

    % Enter desired percent as decimal of S cones per region across
    % quadrants. Follows same form as above
    prBase.quads(2).percentS = [0.05];
    prBase.quads(3).percentS = [0.05];
    prBase.quads(4).percentS = [0.05];
    prBase.quads(5).percentS = [0.1];

    % Establish initial region boundaries in the x and y direction for all
    % four quadrants based on FOV
    prBase.quads(2).xbounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccXDegs;
    prBase.quads(3).xbounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccXDegs;
    prBase.quads(4).xbounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccXDegs;
    prBase.quads(5).xbounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccXDegs;
    prBase.quads(2).ybounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccYDegs;
    prBase.quads(3).ybounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccYDegs;
    prBase.quads(4).ybounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccYDegs;
    prBase.quads(5).ybounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccYDegs;

    if (fullMosaicPercent)
        prBase.quads(5).xbounds = [-prBase.fieldSizeMinutes/60/2 prBase.fieldSizeMinutes/60/2] + prBase.eccXDegs;
        prBase.quads(5).ybounds = [-prBase.fieldSizeMinutes/60/2 prBase.fieldSizeMinutes/60/2] + prBase.eccYDegs;
    end
end

prBase.quads(6).name  = 'overrideQuadSeq'; % Should change this to customQuadSeq
prBase.quads(6).value = true;

% Add indices of cones to be silenced.
prBase.kConeIndices = [];

% Select which quadrants from the above to activate, of the order:
% Q1 Q2 Q3 Q4 Origin
quadSelectList = [[false false false false true]]';%...

% Force build and save of render structures.  This
% only affects this script, and will typically be false.
buildNewForward = false;
buildNewRecon = false;

%% Stimulus parameters.
%
% Size list parameter in degs, expressed as min/60 (because 60 min/deg)
stimSizeDegsList = [3.5] / 60;

% RGB values (before gamma correction)
prBase.stimBgVal = 0.3;% [0.1054 0.1832 0.1189]
stimRValList = 0.80;%0.1054 ./ [2 4 6 8 10];  %[1];% 1.0 0.0];
stimGValList = 0.65;%0.1832 ./ [2 4 6 8 10];  %[1];% 0.0 1.0];
stimBValList = 0.10;%0.1189 ./ [2 4 6 8 10];  %[0];% 0.0 0.0];

% Overwrite stim values to make isoluminant colors on the RG channel.
% Allow for offset from the true isoLum values based on variability in
% staircase procedure where (-) is more red and (+) is more green,
% (Intervals of 50s? 100s? 1000s?)
isoLumRG = true;
% colorStepRG = [-1729 -1280 -880 -480 0 240 480 720 960 1200 1440 1680 1729];
nEquiLumStimuli = 11;


if (isoLumRG)
    switch (prBase.displayName)
        case 'conventional'
            displayFieldName = 'CRT12BitDisplay';
            overwriteDisplayGamma = false; 
            prBase.stimBgVal = [0.1 0.1 0.1];
%             prBase.stimBgVal = [0.1054 0.1832 0.1189]/10;
        case 'mono'
            displayFieldName = 'monoDisplay';
            overwriteDisplayGamma = true;
            prBase.stimBgVal = [0.1 0.1 0.1]/10;% [0.1054 0.1832 0.1189]/10;
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



%     % Grab the display primary xyz values and pull out luminance column
%     eval(['primariesXYZ = displayGet(theDisplayLoad.' displayFieldName ', ''primaries xyz'');']);
%     lumRGB = primariesXYZ(:,2);
% 
%     % Create offset vectors for R and G-raw values, then multiply raw value
%     % by luminance ratio between channels. Find index where R and G
%     % channels are most similar and add step if desired
%     isoLumR = fliplr(0:0.0001:1); isoLumGRaw = 1 - isoLumR;
%     ratioRG = lumRGB(1) / lumRGB(2); alpha = lumRGB(2) / lumRGB(1);
%     isoLumG = ratioRG .* isoLumGRaw;
% 
%     gMax1 = (1-prBase.stimBgVal(1))/alpha; rMax1 = gMax1 * alpha;
%     gMax2 = (1-prBase.stimBgVal(2)); rMax2 = gMax2 * alpha;
% 
%     if gMax1 > 1 || rMax1 > 1
%         rBound = find(isoLumR < rMax2);
%         isoLumR = isoLumR(min(rBound):end);
%         isoLumG = isoLumG(min(rBound):end);
%     elseif gMax2 > 1 || rMax2 > 1
%         rBound = find(isoLumR < rMax1);
%         isoLumR = isoLumR(min(rBound):end);
%         isoLumG = isoLumG(min(rBound):end);
%     else
%         error('Both calculations with background exceed limits')
%     end
% 
%     difRG = abs(isoLumR - isoLumG)';
%     indRG = find(difRG == min(difRG)) + colorStepRG;

    % Overwrite R and G val list based on isolum conditions. To all three
    % channels add background stim value.
%     stimRValList = isoLumR(indRG) + prBase.stimBgVal(1);
%     stimGValList = isoLumG(indRG) + prBase.stimBgVal(2);
%     stimBValList = zeros(1,length(stimRValList)) + prBase.stimBgVal(3);

%     % Clean workspace
%     clear theDisplayLoad; clear theDisplay; clear displayFieldName; clear primariesXYZ;
%     clear isoLumR; clear isoLumG; clear isoLumGRaw; clear difRG; clear indRG;
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
                    for dd = 1:length(forwardChromList)
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

%% Build render structures we need if they are not cached
for pp = 1:length(regPara)

    % Set up paramters structure for this loop, filling in fields that come
    % out of lists precreated above.
    pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
        stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
        forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor);
    pr.quadSelect = quadSelect(:,pp);
    % Compute convenience parameters
    cnv = computeConvenienceParams(pr);

    % Build foward cone mosaic and render matrix if needed
    if (buildNewForward || ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr, cnv, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardNoLCA, ...
            pr.forwardDefocusDiopters, pr.forwardRandSeed, cnv.replaceForwardCones, cnv.forwardStartCones, ...
            cnv.forwardNewCones, pr.forwardEccVars, pr.forwardSubjectID, pr.forwardZernikeDataBase, pr.forwardChrom);
        save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'renderStructure','-v7.3');
        forwardRenderStructure = renderStructure; clear renderStructure;
    end

    % Build recon cone mosaic and render structure if needed
    if (buildNewRecon || ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr, cnv, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconNoLCA, ...
            pr.reconDefocusDiopters, pr.reconRandSeed, cnv.replaceReconCones, cnv.reconStartCones, ...
            cnv.reconNewCones, pr.reconEccVars, pr.reconSubjectID, pr.reconZernikeDataBase, pr.reconChrom);
        save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'renderStructure','-v7.3');
        reconRenderStructure = renderStructure; clear renderStructure;
    end
end

% THIS SHOULD BE A PARFOR AFTERWARDS DON'T FORGET
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

