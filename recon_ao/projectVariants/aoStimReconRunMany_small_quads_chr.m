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

%% Clear
clear; close all;

%% Control size of parpool, otherwise may crush memory
% thePool = gcp('nocreate');
% if (isempty(thePool))
%     parpool(5);
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
prBase.displayGammaBits = 12;
prBase.displayGammaGamma = 2;
displayScaleFactorList = [1];

%% Spatial parameters
% 
% Common to forward and recon models
prBase.nPixels = 100;
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
%    "chromAllL", "chromAllM", "chromAllS", "quadSeq" and number
%    Currently established quadSeq1 - quadSeq6
forwardChromList = ["quadSeq5"]; 
reconChromList =   ["quadSeq5"];

% Build new sequence by
prBase.quads(1).name = 'useQuadSeq';
prBase.quads(1).value = true;

if(prBase.quads(1).value)
    % Initialize storage structure with information on each quadrant
    prBase.quads(2).name = 'Quad1'; 
    prBase.quads(3).name = 'Quad2'; 
    prBase.quads(4).name = 'Quad3'; 
    prBase.quads(5).name = 'Quad4'; 

    % Enter desired percent as decimal of L cones per region across
    % quadrants. The remaining percent will be made of M cones. Entries 
    % should start with outermost regions first and progress inward
    prBase.quads(2).percentL = [0.45]; 
    prBase.quads(3).percentL = [0.50];
    prBase.quads(4).percentL = [0.55];
    prBase.quads(5).percentL = [0.60]; 

    % Enter desired percent as decimal of S cones per region across
    % quadrants. Follows same form as above
    prBase.quads(2).percentS = [0.10]; 
    prBase.quads(3).percentS = [0.10];
    prBase.quads(4).percentS = [0.10];
    prBase.quads(5).percentS = [0.10]; 

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
end

prBase.quads(6).name = 'overrideQuadSeq';
prBase.quads(6).value = true;

% prBase.kConeIndices = [564 505 493 479 528 565 572 585 570];

% Select which quadrants from the above to activate
quadSelectList = [[true true true true]]';%...
%     [true false false false];...
%     [false true false false];...
%     [false false true false];...
%     [false false false true]]';

% Force build and save of render structures.  This
% only affects this script, and will typically be false.
buildNewForward = false;
buildNewRecon = false;

%% Stimulus parameters.
%
% Size list parameter in degs, expressed as min/60 (because 60 min/deg)
stimSizeDegsList = 1 / 60;

% RGB values (before gamma correction)
prBase.stimBgVal = 0.1;
stimRValList = [1.0]% 1.0 0.0 0.80 0.80 0.60];
stimGValList = [1.0]% 0.0 1.0 0.65 0.45 0.65];
stimBValList = [0.0]% 0.0 0.0 0.10 0.10 0.10];

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
shiftInMinutesList = [-3];
fullSquareShift = false;

% Convert the shifts to pixel positions
shiftInPixelsList = round(pixelsPerMinute*shiftInMinutesList);
quadCenters = round(prBase.nPixels / 4);
centerXPosition = prBase.trueCenter + quadCenters + shiftInPixelsList;
centerYPosition = prBase.trueCenter + quadCenters + (shiftInPixelsList(1)*ones(size(centerXPosition)));
%%%% Change the above back to a plus when switching to multiple runs!!!!
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
regParaList = 0.005; %[0.1 0.005 0.001]; %[0.01 0.005 0.001];   % 0.01 0.1 1];
prBase.stride = 2;
prBase.maxReconIterations = 2000;
prBase.whiteNoiseStarts = 0;
prBase.pinkNoiseStarts = 1;
prBase.sparsePriorPatchStarts = 0;
prBase.stimulusStart = false;
prBase.uniformStartVals = [];
prBase.boundedSearch = false;

% Use AO in forward rendering? And determine optics pupil size
prBase.forwardAORender = true;
prBase.forwardNoLCA = true;
prBase.reconAORender = true;
prBase.reconNoLCA = false;
reconPupilDiamListMM =   6;
forwardPupilDiamListMM = 6;

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
forwardDefocusDioptersList = [0.05 0.05 0.05];% 0.05 0.1]; 
reconDefocusDioptersList = [0.12 0.18 0.24];% 0.05 0.1];

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
        renderStructure = buildRenderStruct(pr.aoReconDir , pr.eccXDegs, pr.eccYDegs, ...
            pr.fieldSizeMinutes/60, pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardNoLCA, pr.forwardDefocusDiopters, ...
            cnv.overwriteDisplayGamma, pr.displayName, cnv.displayFieldName, pr.displayGammaBits,  ...
            pr.displayGammaGamma, pr.forwardRandSeed, cnv.replaceForwardCones, cnv.forwardStartCones, ...
            cnv.forwardNewCones, pr.forwardEccVars, pr.forwardSubjectID, pr.forwardZernikeDataBase, pr.quads);
        save(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'renderStructure','-v7.3');
        forwardRenderStructure = renderStructure; clear renderStructure;
    end

    % Build recon cone mosaic and render structure if needed
    if (buildNewRecon || ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr.aoReconDir , pr.eccXDegs, pr.eccYDegs, ...
            pr.fieldSizeMinutes/60, pr.nPixels, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconNoLCA, pr.reconDefocusDiopters, ...
            cnv.overwriteDisplayGamma, pr.displayName, cnv.displayFieldName, pr.displayGammaBits, ...
            pr.displayGammaGamma, pr.reconRandSeed, cnv.replaceReconCones, cnv.reconStartCones, ...
            cnv.reconNewCones, pr.reconEccVars, pr.reconSubjectID, pr.reconZernikeDataBase, pr.quads);
        save(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'renderStructure','-v7.3');
        reconRenderStructure = renderStructure; clear renderStructure;
    end
end

% Run the reconstructions in parallel
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
    aoStimRecon(pr,cnv);
end