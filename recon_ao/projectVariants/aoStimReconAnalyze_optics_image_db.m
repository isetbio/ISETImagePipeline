%% aoStimAnalyze
%
% Description:
%    Call aoStimRecon with many combinations of parameters.
%
% See also: aoStimRecon

% History:
%   11/27/22 dhb  Wrote analysis version.

%% Clear
clear; close all;

%% Set defaults in prBase
prBase = prBaseDefaults;

%% Version editor string
%
% Helps us keep different calcs separate
prBase.versEditor = 'optics_image_db';

%% These parameters are the ones to vary
whichCase = 5;
switch (whichCase)
    case 1
        displayScaleFactorList = [10];
        forwardPupilDiamListMM = [3 3   3 3   3];
        reconPupilDiamListMM =   [2 2.5 3 3.5 4];
        forwardDefocusDioptersList = zeros(size([-2 -1.5 -1 -0.5 0 0.5 1 1.5 2]));
        reconDefocusDioptersList =              [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
        prBase.addPoissonNoise = false;
    case 2
        displayScaleFactorList = [10];
        forwardPupilDiamListMM = [3 3   3 3   3];
        reconPupilDiamListMM =   [2 2.5 3 3.5 4];
        forwardDefocusDioptersList = 0.5*ones(size([-2 -1.5 -1 -0.5 0 0.5 1 1.5 2]));
        reconDefocusDioptersList =                 [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
        prBase.addPoissonNoise = false;
    case 3
        displayScaleFactorList = [1];
        forwardPupilDiamListMM = [3 3   3 3   3];
        reconPupilDiamListMM =   [2 2.5 3 3.5 4];
        forwardDefocusDioptersList = 0.5*ones(size([-2 -1.5 -1 -0.5 0 0.5 1 1.5 2]));
        reconDefocusDioptersList =                 [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
        prBase.addPoissonNoise = false;
    case 4
        displayScaleFactorList = [10];
        reconPupilDiamListMM =   [2 2.5 3 3.5 4];
        forwardPupilDiamListMM = 3.5*ones(size(reconPupilDiamListMM));
        reconDefocusDioptersList =                 [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
        forwardDefocusDioptersList = -0.5*ones(size(reconDefocusDioptersList));
        prBase.addPoissonNoise = true;
    case 5
        displayScaleFactorList = [1];
        reconPupilDiamListMM =   [2 2.5 3 3.5 4];
        forwardPupilDiamListMM = 3.5*ones(size(reconPupilDiamListMM));
        reconDefocusDioptersList =                 [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
        forwardDefocusDioptersList = -0.5*ones(size(reconDefocusDioptersList));
        prBase.addPoissonNoise = true;

end

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
stimSizeDegsList = 0;

% RGB values.  Here we reconstruct a specified image.
%
% This shows how to read in an image and set it up to be reconstructed.
% The example is for a Matlab indexed image.  If it were a straight RGB
% image this would be even easier.
prBase.imageName = 'puppy';
prBase.imageType = 'jpeg';
switch (prBase.imageType)
    case 'jpeg'
        % JPEG options include 'butterfly', 'dragonfly', 'panda', 'pizza'
        % 'puppy', and 'sunflower'.  Not all of these images are square, so
        % we select out the largest square in the center before resizing.
        theImageRGB = imread(fullfile(prBase.aoReconDir,'images',[prBase.imageName '.jpeg']),'jpeg');
        [m,n,k] = size(theImageRGB);
        minDim = min([m,n]);
        mSpace = minDim/2; nSpace = minDim/2;
        lowM = round(m/2-mSpace)+1; highM = lowM+minDim-1; lowN = round(n/2-nSpace)+1; highN = lowN+minDim-1;
        prBase.stimBgVal = imresize(theImageRGB(lowM:highM,lowN:highN,:),'OutputSize',[prBase.nPixels prBase.nPixels]);
        %figure; imshow(prBase.stimBgVal);
    case 'tif'
        % Options are 'zebra'
        theImageRGB = imread(fullfile(prBase.aoReconDir,'images',[prBase.imageName '.tif']),'tif');
        prBase.stimBgVal = imresize(theImageRGB,'OutputSize',[prBase.nPixels prBase.nPixels]);
    case 'png'
        % Options are 'onion'
        theImageRGB = imread([prBase.imageName '.' prBase.imageType]);
        prBase.stimBgVal = imresize(theImageRGB,'OutputSize',[prBase.nPixels prBase.nPixels]);
    case 'matindexed'
        % Options are 'mandrill'
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
% Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
regParaList = 0.0005; %[0.05 0.01 0.005 0.001 0.0005 0.0001]; %[0.01 0.005 0.001];   % 0.01 0.1 1];
prBase.stride = 4;
prBase.maxReconIterations = 1000;
prBase.whiteNoiseStarts = 0;
prBase.pinkNoiseStarts = 1;
prBase.sparsePriorPatchStarts = 0;
prBase.stimulusStart = false;
prBase.uniformStartVals = []; %[ [0.5 0.5 0.5]'  [0.5 0 0]' [0 0.5 0]' [0 0 0.5]' [0 0 0]' [1 1 1]' ];
prBase.boundedSearch = false;

% Use AO in forward rendering? And determine optics pupil size
prBase.forwardAORender = false;
prBase.reconAORender = false;

% Define optics.  Subject only matters if we use a database.
%
% Databases are 'Polans2015' and 'Artal2012'
prBase.forwardSubjectID = 6;
prBase.forwardZernikeDataBase = 'Artal2012';
prBase.reconSubjectID = 6;
prBase.reconZernikeDataBase = 'Artal2012';
% prBase.forwardSubjectID = 0;
% prBase.forwardZernikeDataBase = 'MarimontWandell';
% prBase.reconSubjectID = 0;
% prBase.reconZernikeDataBase = 'MarimontWandell';



% Mosaic chromatic type, options are:
%    "chromNorm", "chromProt", "chromDeut", "chromTrit",
%    "chromAllL", "chromAllM", "chromAllS"
% forwardChromList = ["chromDeut" "chromNorm" "chromNorm"];
% reconChromList =   ["chromDeut" "chromDeut" "chromNorm"];
forwardChromList = ["chromNorm"];
reconChromList =   ["chromNorm"];

% Turn off quads for these calculations
prBase.quads(1).name = 'useQuadSeq';
prBase.quads(1).value = false;

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
                        for pp = 1:length(forwardPupilDiamListMM)
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

% Analyze the reconstructions.  Gather data from each.
for pp = 1:length(regPara)

    % Set up paramters structure for this loop, filling in fields that come
    % out of lists above.
    pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
        stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
        forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor);

    % Compute convenience parameters
    cnv = computeConvenienceParams(pr);

    % Call the driving function
    fprintf('Loading data file %d of %d\n',pp,length(regPara));
    theData = load(fullfile(cnv.outputDir,'xRunOutput.mat'));
    stimulusImageRGB{pp} = theData.stimulusImageRGB;
    stimulusImageLinear{pp} = theData.stimulusImageLinear;
    % reconImageLinear = theData.reconImageLinearTemp;
    reconLogLikelihoods(pp) = theData.multistartStruct.reconLogLikelihoods(theData.reconIndex);
    reconLogPriors(pp) = theData.multistartStruct.reconLogPriors(theData.reconIndex);
    reconLosses(pp) = theData.multistartStruct.reconLosses(theData.reconIndex);
    reconScaleFactor(pp) = theData.reconScaleFactor(theData.reconIndex);
    forwardExcitations(:,pp) = theData.forwardExcitationsToStimulusUse;
    reconExcitations(:,pp) = theData.multistartStruct.reconPreds(:,theData.reconIndex);
    runForwardPupilDiameterMM(pp) = pr.forwardPupilDiamMM;
    runReconPupilDiameterMM(pp)= pr.reconPupilDiamMM;
    runForwardDefocusDiopters(pp) = pr.forwardDefocusDiopters;
    runReconDefocusDiopters(pp) = pr.reconDefocusDiopters;
    R = corrcoef(forwardExcitations(:,pp),reconExcitations(:,pp));
    corrs(pp) = R(1,2);
end

% Get the forward values used
forwardPupilVal = unique(runForwardPupilDiameterMM);
if (length(forwardPupilVal) ~= 1)
    error('More than one forward pupil value');
end
forwardDefocusVal = unique(runForwardDefocusDiopters);
if (length(forwardDefocusVal) ~= 1)
    error('More than one forward defocus value');
end
index = find(runReconPupilDiameterMM == forwardPupilVal & runReconDefocusDiopters == forwardDefocusVal);
negLogLossAtForwardVals = -reconLosses(index);

% Get the set of recon pupil values used
reconPupilVals = unique(runReconPupilDiameterMM);
reconDefocusVals = unique(runReconDefocusDiopters);
for xx = 1:length(reconPupilVals)
    for yy = 1:length(reconDefocusVals)
        X(xx,yy) = reconPupilVals(xx);
        Y(xx,yy) = reconDefocusVals(yy);
        index = find(runReconPupilDiameterMM == reconPupilVals(xx) & runReconDefocusDiopters == reconDefocusVals(yy));
        Z(xx,yy) = -reconLosses(index);
    end
end

% Plot the negative log loss surface
figure; hold on
surf(X,Y,Z);

% Indicate location of max neg log loss
[~,index] = max(-reconLosses);
plot3(runReconPupilDiameterMM(index),runReconDefocusDiopters(index),-reconLosses(index),'ro','MarkerFaceColor','r','MarkerSize',12);

% Indicate forward vals
plot3(forwardPupilVal,forwardDefocusVal,negLogLossAtForwardVals,'go','MarkerFaceColor','g','MarkerSize',8);

% Tidy up plot
zlim([0.9*max(-reconLosses), 1.1*max(-reconLosses)]);
xlabel('Pupil Diameter MM');
ylabel('Defocus Diopters');
title('Neg Log Loss');    
view(-20,60);

% figure; plot(reconPupilDiamMM,-reconLosses,'ro'); title('Loss');
% figure; plot(reconPupilDiamMM,reconLogPriors,'ro'); title('Prior');
% figure; plot(reconPupilDiamMM,reconLogLikelihoods,'ro'); title('Likelihood');
% figure; plot(reconPupilDiamMM,corrs,'ro'); title('Correlation');