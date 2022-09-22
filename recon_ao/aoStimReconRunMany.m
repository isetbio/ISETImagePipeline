%% aoStimReconMany
%
% Description:
%    Call aoStimRecon with many combinations of parameters.
%
% See also: aoStimRecon

% History:
%   08/15/22  dhb  Wrote after converting aoStimRecon to a function
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options

%% Clear
clear; close all;

%% Parameters
%
% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
displayName = 'mono';

% Stimulus parameters.
%
% Size list parameter in degs, but expressed as min/60 (because 60 min/deg)
stimSizeDegsList = [24/60]

% RGB values (before gamma correction) 
stimBgVal = 0.1;
% Below are new values, first is for uniform field recon in dichrom, second
% is the corresponding metamer in dichrom conditions
stimRValList = 0.80; %[0.1620 0.9281]; % 0.6466];
stimGValList = 0.65; %[0.8461 0.6745]; % 0.7015];
stimBValList = 0.10; %[0.9490 0.9492]; % 0.0994];
if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
    error('Stimulus value lists must have same length');
end

% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional
%                           display.
sparsePriorStr = 'conventional';

% Reconstruction parameters
regParaList = [0.001];% 0.01 0.1 1];
stride = 2;

% Use AO in forward rendering?
%
% This determines pupil diameter which typically differs in AO 
forwardAORender = [true];
reconAORender = [true];

% Residual defocus for forward and recon rendering, of equal sizes
forwardDefocusDioptersList = [0.00];% 0.05 0.1]; 
reconDefocusDioptersList = [0.00];% 0.05 0.1];

% Establish chromaticity for forward and recon mosaic, with string options:
% "chromNorm", "chromProt", "chromDeut", "chromTrit", "chromAllL", "chromAllM"
forwardChromList = ["chromNorm"];% "chromNorm" "chromNorm"]; 
reconChromList = ["chromNorm"];% "chromDeut" "chromNorm"];

% Establish a slide list which adjusts the stimulus position across the
% mosaic diagonally in non-overlapping portions
slideList = -3:3; 



%% Run through specified list conditions
for ss = 1:length(stimSizeDegsList)
    stimSizeDegs = stimSizeDegsList(ss);
    for cc = 1:length(stimRValList)
        stimRVal = stimRValList(cc);
        stimGVal = stimGValList(cc);
        stimBVal = stimBValList(cc);
        for ff = 1:length(forwardDefocusDioptersList)
            forwardDefocusDiopters = forwardDefocusDioptersList(ff);
            reconDefocusDiopters = reconDefocusDioptersList(ff);
            for rr = 1:length(regParaList)
                regPara = regParaList(rr);
                for dd = 1:length(forwardChromList)
                    forwardChrom = forwardChromList(dd);
                    reconChrom = reconChromList(dd);
                    for yy = 1:length(slideList)
                        slide = slideList(yy);
                        aoStimRecon(displayName,sparsePriorStr,...
                            forwardAORender, reconAORender, ...
                            forwardDefocusDiopters, reconDefocusDiopters, ...
                            stimSizeDegs,stimBgVal,stimRVal,stimGVal,stimBVal,...
                            regPara,stride, forwardChrom, reconChrom, slide);
                    end
                end
            end
        end
    end
end