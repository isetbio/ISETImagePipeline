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

%% Parameters
%
% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
displayName = 'conventional';
versEditor = 'dichrom_db';

%% Spatial parameters
% 
% Common to forward and recon models
nPixels = 100;
trueCenter = round(nPixels/2);
forwardEccVars = false;
reconEccVars = false;

%% Stimulus parameters.
%
% Size list parameter in degs, expressed as min/60 (because 60 min/deg)
stimSizeDegsList = 0.5; %[24/60];

% RGB values (before gamma correction) 
stimBgVal = 0.1;
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
centerXPosition = [trueCenter trueCenter];
centerYPosition = [trueCenter trueCenter];
stimCenter = [centerXPosition ; centerYPosition];
deltaCenter = stimCenter - trueCenter;

%% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional display.
sparsePriorStr = 'conventional';

%% Reconstruction parameters
%
% Should cycle through a few of these regs to optimize for 58x58 pixels
% Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
regParaList = [0.01 0.005 0.001];   % 0.01 0.1 1];
stride = 2;
maxReconIterations = 1000;

% Use AO in forward rendering? Should consider mix-and-match 
%
% This determines pupil diameter which typically differs in AO 
forwardAORender = false;
reconAORender = false;

% Residual defocus for forward and recon rendering, of equal sizes
forwardDefocusDioptersList = [0.00];% 0.05 0.1]; 
reconDefocusDioptersList = [0.00];% 0.05 0.1];

% Mosaic chromatic type, options are:
%    "chromNorm", "chromProt", "chromDeut", "chromTrit", 
%    "chromAllL", "chromAllM", "chromAllS"
forwardChromList = ["chromDeut", "chromNorm", "chromNorm"]; 
reconChromList = ["chromDeut", "chromDeut", "chromNorm"];

%% Set up list conditions
runIndex = 1;
for ss = 1:length(stimSizeDegsList)
    for cc = 1:length(stimRValList)
        for yy = 1:length(deltaCenter)
            for ff = 1:length(forwardDefocusDioptersList)
                for rr = 1:length(regParaList)
                    for dd = 1:length(forwardChromList)

                        stimSizeDegs(runIndex) = stimSizeDegsList(ss);

                        stimRVal(runIndex) = stimRValList(cc);
                        stimGVal(runIndex) = stimGValList(cc);
                        stimBVal(runIndex) = stimBValList(cc);

                        stimCenter(:,runIndex) = deltaCenter(:,yy);

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

% Run them all in parallel
parfor pp = 1:length(regPara)
    aoStimRecon(displayName,sparsePriorStr,...
        forwardAORender, reconAORender, ...
        forwardDefocusDiopters(pp), reconDefocusDiopters(pp), ...
        stimSizeDegs(pp),stimBgVal,stimRVal(pp),stimGVal(pp),stimBVal(pp),...
        regPara(pp), stride, forwardChrom(pp), reconChrom(pp), ...
        stimCenter(:,pp), trueCenter, nPixels, versEditor, maxReconIterations, ...);
        forwardEccVars, reconEccVars);
end