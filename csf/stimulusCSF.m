function imageStim = stimulusCSF(stimType, stimCrst, stimFreq) 
% Illustrate how to generate scenes depicting various cone-specific stimuli
%
% Syntax:
%   t_generateConeSpecificStimuli
%
% Description:
%    Simple script that demonstrates how to generate scenes depicting 
%    various cone-specific stimuli on specific background. The script uses 
%    plotlab,for plotting the results. 
%    Plotlab is freely available at: https://github.com/npcottaris/plotlab
%
% Inputs:
%    stimType: direction of chromatic modulation for the stimulus
%    {'S-(L+M)', 'L-M', 'L+M', 'L+M+S', 'all-different'}  
%    stimCrst: stimulus contrast
%    stimFreq: stimulus frequency
%
% Outputs:
%    None.
% History:
%    03/28/20  npc  Wrote it.
%    09/16/20  lq   Changed it to a function

    %% Specify the desired mean (x,y) chromaticity and luminance (background)
    % Horizon light - see https://en.wikipedia.org/wiki/Standard_illuminant
    background.xyChroma = [0.345 0.358];
    % 40 cd/m2
    background.luminance = 40;
    % In vector form
    background.xyY = [background.xyChroma(1) background.xyChroma(2) background.luminance];
    
    %% Specify stimulus type
    test.spatialModulationParams = getModulationParamsForStimType(stimType, stimCrst, stimFreq);
            
    %% Stimulus field of view
    fieldOfViewDegs = 4;
    
    %% Presentation display
    presentationDisplay = generatePresentationDisplay(); 
    
    %% Background XYZ tri-stimulus values
    background.XYZ = (xyYToXYZ(background.xyY(:)))';
    
    %% Background linear RGB primary values for the presentation display
    background.RGB = imageLinearTransform(background.XYZ, inv(displayGet(presentationDisplay, 'rgb2xyz')));
    
    %% Background LMS excitations
    background.LMS = imageLinearTransform(background.RGB, displayGet(presentationDisplay, 'rgb2lms'));
    
    %% Stimulus spatial modulation of the L-, M-, and S-cone contrast
    test.LMScontrastImage = generateSpatialContrastPatterns(fieldOfViewDegs, test.spatialModulationParams);
    
    %% Stimulus LMS excitations image for the given background and spatial modulation
    test.LMSexcitationImage = bsxfun(@times, (1+test.LMScontrastImage), reshape(background.LMS, [1 1 3]));
    
    %% Stimulus linear RGB primaries image
    test.RGBimage = imageLinearTransform(test.LMSexcitationImage, inv(displayGet(presentationDisplay, 'rgb2lms')));
    
    %% Make sure we are in gamut (no subpixels with primary values outside of [0 1]
    assert((numel(find(test.RGBimage>1))==0)&&(numel(find(test.RGBimage<0))==0), ...
        sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
        numel(find(test.RGBimage>1)), numel(find(test.RGBimage<0))));
    
    %% Gamma correction
    % Gamma correct linear RGB values (primaries) through the display's
    % gamma table to get RGB settings values
    test.RGBimageGammaCorrected = ieLUTLinear(test.RGBimage, ...
        displayGet(presentationDisplay, 'inverse gamma'));
    
    test.RGBimageGammaCorrected = test.RGBimageGammaCorrected / max(test.RGBimageGammaCorrected(:));

    %% Generate scene corresponding to the test stimulus on the presentation display
    theScene = sceneFromFile(test.RGBimageGammaCorrected,'rgb',background.luminance, presentationDisplay);
    theScene = sceneSet(theScene, 'h fov', fieldOfViewDegs);

    %% Verify that we have the desired cone contrast profiles
    % Get the emitted radiance image
    emittedRadianceImage = sceneGet(theScene, 'energy');
    
    %% Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayGet(presentationDisplay, 'wave'));
    
    %% Compute the LMS cone contrasts of the emitted radiance image
    test.achievedLMScontrastImage = computeLMScontrastImage(emittedRadianceImage, coneFundamentals);

    % Return RGM image of the stimulus
    imageStim = test.RGBimage;
end

function LMScontrastImage = computeLMScontrastImage(radianceImage, coneFundamentals)
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
        
    coneExcitationsImage = radianceImage * coneFundamentals;
    coneExcitationsBackground = coneExcitationsImage(1,:);

    LMScontrastImage = bsxfun(@times, ...
        bsxfun(@minus, coneExcitationsImage, coneExcitationsBackground), ...
        1./coneExcitationsBackground);
    
    LMScontrastImage = reshape(LMScontrastImage, [rowsNum colsNum 3]);
end

function LMScontrastImage = generateSpatialContrastPatterns(fieldOfViewDegs, spatialModulationParams)
    x = linspace(-fieldOfViewDegs/2,fieldOfViewDegs/2,256);
    [X,Y] = meshgrid(x);
    coneIDs = keys(spatialModulationParams);
    for coneIndex = 1:numel(coneIDs)
        p = spatialModulationParams(coneIDs{coneIndex});
        contrastPattern = cos(2*pi*(p.spatialFrequency(1)*(X-p.gaborPos(1)) + p.spatialFrequency(2)*(Y-p.gaborPos(2))));
        contrastEnvelope = exp(-((X-p.gaborPos(1))/p.gaborSigma(1)).^2) .* exp(-((Y-p.gaborPos(2))/p.gaborSigma(2)).^2);
        LMScontrastImage(:,:,coneIndex) = p.maxContrast * contrastPattern .* contrastEnvelope;
    end
end

function presentationDisplay = generatePresentationDisplay()
    % Generate presentation display
    presentationDisplay = displayCreate('LCD-Apple');
    bitDepth = 16;
    gammaT = 'nonlinear';
    if (strcmp(gammaT, 'linear'))
        N = 2^bitDepth;
        gTable = repmat(linspace(0, 1, N), 3, 1)';
    else
        x = linspace(0,1,2^bitDepth );
        gTable = x(:).^2.1;
        gTable = repmat(gTable, [1,3]);
    end
    presentationDisplay = displaySet(presentationDisplay, 'gTable', gTable);
end

function spatialModulationParams = getModulationParamsForStimType(stimType, stimCrst, stimFreq)
    spatialModulationParams = containers.Map();
    switch (stimType)
        case 'all-different'
            spatialModulationParams('L') = struct('maxContrast',  0.1 * stimCrst, 'gaborPos', [0 0], 'gaborSigma', [Inf 0.6], 'spatialFrequency', [1.5 0]);
            spatialModulationParams('M') = struct('maxContrast', -0.1 * stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [0 1]);
            spatialModulationParams('S') = struct('maxContrast', -0.8 * stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.2 1.2], 'spatialFrequency', [0.6 0]);
        case 'S-(L+M)'
            spatialModulationParams('L') = struct('maxContrast', -0.2 * stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('M') = struct('maxContrast', -0.2 * stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.8 * stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
        case 'L-M'
            spatialModulationParams('L') = struct('maxContrast',  stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('M') = struct('maxContrast', -stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.0, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
        case 'L+M'
            spatialModulationParams('L') = struct('maxContrast',  stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('M') = struct('maxContrast',  stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.0, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
        case 'L+M+S'
            spatialModulationParams('L') = struct('maxContrast',  stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('M') = struct('maxContrast',  stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
            spatialModulationParams('S') = struct('maxContrast',  stimCrst, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [stimFreq 0]);
        otherwise
            error('Not know stimType: ''%s''.', stimType);
    end
end
