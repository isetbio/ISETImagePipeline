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
versEditor = 'dichrom';

%% Spatial parameters
% 
% Common to forward and recon models
nPixels = 58;
trueCenter = round(nPixels/2);
eccVars = true;

%% Stimulus parameters.
%
% Size list parameter in degs, expressed as min/60 (because 60 min/deg)
stimSizeDegsList = [24/60];

% RGB values (before gamma correction) 
stimBgVal = 0.1;
stimRValList = [0.0110 0.9499];
stimGValList = [0.6180 0.1729];
stimBValList = [0.9667 0.9732];

% Check that all channels receive same number of inputs
if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
    error('Stimulus value lists must have same length');
end

% Input desired x and y position for stimulus to be centered over. Function
% will end if values exceed pixel limits. 
%
% Position specified in pixels, could consider specifying in degrees.
centerXPosition = [trueCenter];
centerYPosition = [trueCenter];
stimCenter = [centerXPosition; centerYPosition];
deltaCenter = stimCenter - trueCenter;

%% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional
%                           display.
sparsePriorStr = 'conventional';

%% Reconstruction parameters
%
% Should cycle through a few of these regs to optimize for 58x58 pixels
% Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
regParaList = [0.001];   % 0.01 0.1 1];
stride = 2;

% Use AO in forward rendering? Should consider mix-and-match 
%
% This determines pupil diameter which typically differs in AO 
forwardAORender = [true];
reconAORender = [true];

% Residual defocus for forward and recon rendering, of equal sizes
forwardDefocusDioptersList = [0.00];% 0.05 0.1]; 
reconDefocusDioptersList = [0.00];% 0.05 0.1];

% Mosaic chromatic type, options are:
%    "chromNorm", "chromProt", "chromDeut", "chromTrit", 
%    "chromAllL", "chromAllM", "chromAllS"
forwardChromList = ["chromDeut", "chromNorm", "chromNorm"]; 
reconChromList = ["chromDeut", "chromDeut", "chromNorm"];

%% Run through specified list conditions
for ss = 1:length(stimSizeDegsList)
    stimSizeDegs = stimSizeDegsList(ss);
    for cc = 1:length(stimRValList)
        stimRVal = stimRValList(cc);
        stimGVal = stimGValList(cc);
        stimBVal = stimBValList(cc);
        for yy = 1:length(deltaCenter)
            stimCenter = deltaCenter(:,yy);
            for ff = 1:length(forwardDefocusDioptersList)
                forwardDefocusDiopters = forwardDefocusDioptersList(ff);
                reconDefocusDiopters = reconDefocusDioptersList(ff);
                for rr = 1:length(regParaList)
                    regPara = regParaList(rr);
                    for dd = 1:length(forwardChromList)
                        forwardChrom = forwardChromList(dd);
                        reconChrom = reconChromList(dd);
                        aoStimRecon(displayName,sparsePriorStr,...
                            forwardAORender, reconAORender, ...
                            forwardDefocusDiopters, reconDefocusDiopters, ...
                            stimSizeDegs,stimBgVal,stimRVal,stimGVal,stimBVal,...
                            regPara,stride, forwardChrom, reconChrom, ...
                            stimCenter, nPixels, trueCenter, eccVars, versEditor);
                    end
                end
            end
        end
    end
end