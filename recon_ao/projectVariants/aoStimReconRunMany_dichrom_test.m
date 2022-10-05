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

%% Clear
clear; close all;

%% Version editor string
%
% Helps us keep different calcs separate
prBase.versEditor = 'dichrom_test';

%% Point at directory with data files for this subproject
%
% This will allow us to load in project specific precomputed information.
% Also records initials of version editors, otherwise set to 'main'
prBase.aoReconDir = getpref('ISETImagePipeline','aoReconDir');

%% Parameters
%
% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
prBase.displayName = 'conventional';
prBase.displayGammaBits = 12;
prBase.displayGammaGamma = 2;

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

%% Stimulus parameters.
%
% Size list parameter in degs, expressed as min/60 (because 60 min/deg)
stimSizeDegsList = 0.5; %[24/60];

% RGB values (before gamma correction) 
prBase.stimBgVal = 0.1;
stimRValList = [1.1048e-02 8.2258e-01];
stimGValList = [6.1803e-01 3.4322e-01];
stimBValList = [9.6667e-01 9.7158e-01];

% Check that all channels receive same number of inputs
if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
    error('Stimulus value lists must have same length');
end

% Input desired x and y position for stimulus to be centered over. Function
% will end if values exceed pixel limits. 
%
% Position specified in pixels, could consider specifying in degrees.
centerXPosition = [prBase.trueCenter];
centerYPosition = [prBase.trueCenter];
prBase.stimCenter = [centerXPosition ; centerYPosition];
deltaCenterList = [prBase.stimCenter - prBase.trueCenter];

%% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional display.
prBase.sparsePriorStr = 'conventional';

%% Reconstruction parameters
%
% Should cycle through a few of these regs to optimize for 58x58 pixels
% Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
regParaList = 0.002; %[0.01 0.005 0.001];   % 0.01 0.1 1];
prBase.stride = 2;
prBase.maxReconIterations = 5;
prBase.whiteNoiseStarts = 0;
prBase.pinkNoiseStarts = 0;
prBase.sparsePriorPatchStarts = 0;
prBase.uniformStartVals = [0.5];

% Use AO in forward rendering? Should consider mix-and-match 
%
% This determines pupil diameter which typically differs in AO 
prBase.forwardAORender = false;
prBase.reconAORender = false;

% Residual defocus for forward and recon rendering, of equal sizes
forwardDefocusDioptersList = [0.00];% 0.05 0.1]; 
reconDefocusDioptersList = [0.00];% 0.05 0.1];

% Mosaic chromatic type, options are:
%    "chromNorm", "chromProt", "chromDeut", "chromTrit", 
%    "chromAllL", "chromAllM", "chromAllS"
forwardChromList = ["chromDeut"]; 
reconChromList = ["chromDeut"];

% Force build and save of render structures.  This
% only affects this script, and will typically be false.
buildNewForward = false;
buildNewRecon = false;

%% Set up list conditions
runIndex = 1;
for ss = 1:length(stimSizeDegsList)
    for cc = 1:length(stimRValList)
        for yy = 1:length(deltaCenterList)
            for ff = 1:length(forwardDefocusDioptersList)
                for rr = 1:length(regParaList)
                    for dd = 1:length(forwardChromList)

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

                        runIndex = runIndex + 1;
                    end
                end
            end
        end
    end
end

%% Build render structures we need if they are not cached
% Run them all in parallel
for pp = 1:length(regPara)

    % Set up paramters structure for this loop, filling in fields that come
    % out of lists above.
    pr = prBase;
    pr.stimSizeDegs = stimSizeDegs(pp);
    pr.stimRVal = stimRVal(pp);
    pr.stimGVal = stimGVal(pp);
    pr.stimBVal = stimBVal(pp);
    pr.stimCenter = stimCenter(:,pp);
    pr.forwardDefocusDiopters = forwardDefocusDiopters(pp);
    pr.reconDefocusDiopters = reconDefocusDiopters(pp);
    pr.regPara = regPara(pp);
    pr.forwardChrom = forwardChrom(pp);
    pr.reconChrom = reconChrom(pp);

    % Compute convenience parameters
    cnv = computeConvenienceParams(pr);

    % Build foward cone mosaic and render matrix if needed
    if (buildNewForward || ~exist(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr.aoReconDir , pr.eccXDegs, pr.eccYDegs, ...
            pr.fieldSizeMinutes/60, pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters, ...
            cnv.overwriteDisplayGamma, pr.displayName, cnv.displayFieldName, pr.displayGammaBits, ...
            pr.displayGammaGamma, pr.forwardRandSeed, cnv.replaceForwardCones, cnv.forwardStartCones, ...
            cnv.forwardNewCones, pr.forwardEccVars);
        save(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'renderStructure');
        forwardRenderStructure = renderStructure; clear renderStructure;
    end

    % Build recon cone mosaic and render structure if needed
    if (buildNewRecon || ~exist(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr.aoReconDir , pr.eccXDegs, pr.eccYDegs, ...
            pr.fieldSizeMinutes/60, pr.nPixels, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconDefocusDiopters, ...
            cnv.overwriteDisplayGamma, pr.displayName, cnv.displayFieldName, pr.displayGammaBits, ...
            pr.displayGammaGamma, pr.reconRandSeed, cnv.replaceReconCones, cnv.reconStartCones, ...
            cnv.reconNewCones, pr.reconEccVars);
        save(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'renderStructure');
        reconRenderStructure = renderStructure; clear renderStructure;
    end
end

% Run them all in parallel
for pp = 1:length(regPara)

    % Set up paramters structure for this loop, filling in fields that come
    % out of lists above.
    pr = prBase;
    pr.stimSizeDegs = stimSizeDegs(pp);
    pr.stimRVal = stimRVal(pp);
    pr.stimGVal = stimGVal(pp);
    pr.stimBVal = stimBVal(pp);
    pr.stimCenter = stimCenter(:,pp);
    pr.forwardDefocusDiopters = forwardDefocusDiopters(pp);
    pr.reconDefocusDiopters = reconDefocusDiopters(pp);
    pr.regPara = regPara(pp);
    pr.forwardChrom = forwardChrom(pp);
    pr.reconChrom = reconChrom(pp);

    % Compute convenience parameters
    cnv = computeConvenienceParams(pr);

    % Call the driving function
    aoStimRecon(pr,cnv);
end