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
%   10/05/22  dhb  Lots of changes for parall

%% Clear
clear; close all;

%% Set defaults in prBase
prBase = prBaseDefaults;

%% Version editor string
%
% Helps us keep different calcs separate
prBase.versEditor = 'dichrom_image_chr';

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
prBase.nPixels = 128;
prBase.trueCenter = round(prBase.nPixels/2);

%% Mosaic parameters
prBase.fieldSizeMinutes = 60;
prBase.eccXDegs = 0.0;
prBase.eccYDegs = 0.0;
prBase.forwardRandSeed = false;
prBase.reconRandSeed = false;
prBase.forwardEccVars = false;
prBase.reconEccVars = false;
prBase.reconstructfromRenderMatrix = true;
prBase.addPoissonNoise = false;

%% Stimulus parameters.
%
% Size list parameter in degs, expressed as min/60 (because 60 min/deg)
%
% Because we're specifying an image below, this is a dummy parameter
% in this script.
stimSizeDegsList = 0.5; %[24/60];

% RGB values.  Here we reconstruct a specified image.
%
% This shows how to read in an image and set it up to be reconstructed.
% The example is for a Matlab indexed image.  If it were a straight RGB
% image this would be even easier.
%
%prBase.imageName = 'mandrill';
%prBase.imageType = 'matindexed';
% prBase.imageName = 'zebra';
% prBase.imageType = 'tif';
prBase.imageName = 'onion';
prBase.imageType = 'png';
switch (prBase.imageType)
    case 'tif'
        theImageRGB = imread(fullfile(prBase.aoReconDir,'images',[prBase.imageName '.tif']),'tif');
        prBase.stimBgVal = imresize(theImageRGB,'OutputSize',[prBase.nPixels prBase.nPixels]);
    case 'png'
        theImageRGB = imread([prBase.imageName '.' prBase.imageType]);
        prBase.stimBgVal = imresize(theImageRGB,'OutputSize',[prBase.nPixels prBase.nPixels]);
    case 'matindexed'
        rawImage = load([prBase.imageName '.mat']);
        theImageRGB = zeros(size(rawImage.X,1),size(rawImage.X,2),3);
        for ii = 1:size(rawImage.X,1)
            for jj = 1:size(rawImage.X,2)
                for kk = 1:3
                    theImageRGB(ii,jj,kk) = rawImage.map(rawImage.X(ii,jj),kk);
                end
            end
        end
        prBase.stimBgVal = imresize(theImageRGB,'OutputSize',[prBase.nPixels prBase.nPixels]);
end

% We store the image in the stimBgVal field above, which asStimRecon understands.
% We dummy up R, G, and BVal lists to have one entry each.
stimRValList = [1];
stimGValList = [1];
stimBValList = [1];

% Check that all channels receive same number of inputs
if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
    error('Stimulus value lists must have same length');
end

% Input desired x and y position for stimulus to be centered over. Function
% will end if values exceed pixel limits.
%
% Position specified in pixels, could consider specifying in degrees.
%
% Because we're specifying an image, this is a dummy parameter and the
% list should have length 1.
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
regParaList = 0.01; %[0.01 0.005 0.001];   % 0.01 0.1 1];
prBase.stride = 2;
prBase.maxReconIterations = 5000;
prBase.whiteNoiseStarts = 0;
prBase.pinkNoiseStarts = 1;
prBase.sparsePriorPatchStarts = 0;
prBase.stimulusStart = false;
prBase.uniformStartVals = []; %[ [0.5 0.5 0.5]'  [0.5 0 0]' [0 0.5 0]' [0 0 0.5]' [0 0 0]' [1 1 1]' ];
prBase.boundedSearch = true;

% Use AO in forward rendering? And determine optics pupil size
prBase.forwardAORender = false;
prBase.reconAORender = false;
prBase.forwardPupilDiamMM = 3;
prBase.reconPupilDiamMM = 3;
prBase.forwardSubjectID = 6;
prBase.forwardZernikeDataBase = 'Polans2015';
prBase.reconSubjectID = 6;
prBase.reconZernikeDataBase = 'Polans2015';

% Residual defocus for forward and recon rendering, of equal sizes
forwardDefocusDioptersList = [0.00];% 0.05 0.1];
reconDefocusDioptersList = [0.00];% 0.05 0.1];

% Mosaic chromatic type, options are:
%    "chromNorm", "chromProt", "chromDeut", "chromTrit",
%    "chromAllL", "chromAllM", "chromAllS"
forwardChromList = ["chromDeut" "chromNorm" "chromNorm"];
reconChromList =   ["chromDeut" "chromDeut" "chromNorm"];

% Force build and save of render structures.  This
% only affects this script, and will typically be false.
buildNewForward = false;
buildNewRecon = false;

%% Set up list conditions
runIndex = 1;
for ss = 1:length(stimSizeDegsList)
    for cc = 1:length(stimRValList)
        for yy = 1:size(deltaCenterList,2)
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
for pp = 1:length(regPara)

    % Set up paramters structure for this loop, filling in fields that come
    % out of lists precreated above.
    pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
        stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
        forwardChrom,reconChrom);

    % Compute convenience parameters
    cnv = computeConvenienceParams(pr);

    % Build foward cone mosaic and render matrix if needed
    if (buildNewForward || ~exist(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr.aoReconDir , pr.eccXDegs, pr.eccYDegs, ...
            pr.fieldSizeMinutes/60, pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters, ...
            cnv.overwriteDisplayGamma, pr.displayName, cnv.displayFieldName, pr.displayGammaBits, ...
            pr.displayGammaGamma, pr.forwardRandSeed, cnv.replaceForwardCones, cnv.forwardStartCones, ...
            cnv.forwardNewCones, pr.forwardEccVars, pr.forwardSubjectID, pr.forwardZernikeDataBase);
        save(fullfile(cnv.renderDir , cnv.forwardRenderStructureName),'renderStructure', '-v7.3');
        forwardRenderStructure = renderStructure; clear renderStructure;
    end

    % Build recon cone mosaic and render structure if needed
    if (buildNewRecon || ~exist(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'file'))
        renderStructure = buildRenderStruct(pr.aoReconDir , pr.eccXDegs, pr.eccYDegs, ...
            pr.fieldSizeMinutes/60, pr.nPixels, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconDefocusDiopters, ...
            cnv.overwriteDisplayGamma, pr.displayName, cnv.displayFieldName, pr.displayGammaBits, ...
            pr.displayGammaGamma, pr.reconRandSeed, cnv.replaceReconCones, cnv.reconStartCones, ...
            cnv.reconNewCones, pr.reconEccVars, pr.reconSubjectID, pr.reconZernikeDataBase);
        sizeCheck = whos('renderStructure');
        save(fullfile(cnv.renderDir , cnv.reconRenderStructureName),'renderStructure', '-v7.3');
        reconRenderStructure = renderStructure; clear renderStructure;
    end
end

% Run the reconstructions in parallel
parfor pp = 1:length(regPara)

    % Set up paramters structure for this loop, filling in fields that come
    % out of lists above.
    pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
        stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
        forwardChrom,reconChrom);

    % Compute convenience parameters
    cnv = computeConvenienceParams(pr);

    % Call the driving function
    aoStimRecon(pr,cnv);
end