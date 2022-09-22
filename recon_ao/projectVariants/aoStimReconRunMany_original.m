%% aoStimReconMany
%
% Description:
%    Call aoStimRecon with many combinations of parameters.
%
% See also: aoStimRecon

% History:
%   08/15/22  dhb  Wrote after converting aoStimRecon to a function

%% Clear
clear; close all;

%% Parameters
%
% Display, options are:
%    'conventional'    - A conventional display
%    'mono'            - A display with monochromatic primaries
displayName = 'mono';

% Stimulus parameters
stimSizeDegsList = [0.4 10/60];
stimBgVal = 0.1;
stimRValList = [0.80  0.70];
stimGValList = [0.65  0.70];
stimBValList = [0.10  0.10];
if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
    error('Stimulus value lists must have same length');
end

% Prior parameters
%
% conventionalSparsePrior - from the paper, images analyzed on conventional
%                           display.
sparsePriorStr = 'conventional';

% Reconstruction parameters
regParaList = [0.001 0.01 0.1 1];
stride = 2;

% Use AO in forward rendering?
%
% This determines pupil diameter which typically differs in AO
forwardAORender = true;

% Residual defocus for forward rendering
forwardDefocusDioptersList = [0.00 0.05 0.1];

%% Run through specified list conditions
for ss = 1:length(stimSizeDegsList)
    stimSizeDegs = stimSizeDegsList(ss);
    for cc = 1:length(stimRValList)
        stimRVal = stimRValList(cc);
        stimGVal = stimGValList(cc);
        stimBVal = stimBValList(cc);
        for ff = 1:length(forwardDefocusDioptersList)
            forwardDefocusDiopters = forwardDefocusDioptersList(ff);
            for rr = 1:length(regParaList)
                regPara = regParaList(rr);
                aoStimRecon(displayName,sparsePriorStr,...
                    forwardAORender, ...
                    forwardDefocusDiopters, ...
                    stimSizeDegs,stimBgVal,stimRVal,stimGVal,stimBVal,...
                    regPara,stride);
            end
        end
    end
end