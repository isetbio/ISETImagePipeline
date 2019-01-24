%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
fovDeg = 1;
theMosaic = coneMosaicHex(5, ...              % hex lattice sampling factor
   'fovDegs', fovDeg, ...                     % match mosaic width to stimulus size
   'eccBasedConeDensity', true, ...           % cone density varies with eccentricity
   'eccBasedConeQuantalEfficiency', false, ... % cone quantal efficiency varies with eccentricity
   'integrationTime', 50, ...               % 0.1s integration time
   'maxGridAdjustmentIterations', 50);        % terminate iterative lattice adjustment after 50 iterations

theMosaic.visualizeGrid(...
    'backgroundColor', [1 1 1], ...
    'ticksInVisualDegs', true);

%% Create scene 
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.57);
meanLuminanceCdPerM2 = 100;

% Load picture
fileName = fullfile(isetbioDataPath, 'images', 'rgb', 'eagle.jpg');
realizedStimulusScene = sceneFromFile(fileName, 'rgb', ...
    meanLuminanceCdPerM2, presentationDisplay);

% Set the angular scene width
realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', fovDeg);

% Visualize different aspects of the generated scene
visualizeScene(realizedStimulusScene);

%% Human optics
theOI = oiCreate('wvf human');
theOI = oiCompute(theOI, realizedStimulusScene);

% Visualize different aspects of the computed optical image
visualizeOpticalImage(theOI, 'displayRetinalContrastProfiles', true);

%% Cone excitations
nTrialsNum = 10;
emPath = zeros(nTrialsNum, 1, 2);

% Compute mosaic excitation responses
coneExcitations = theMosaic.compute(theOI, 'emPath', emPath);
visualizeConeMosaicResponses(theMosaic, coneExcitations, 'R*/cone/tau');

%% Default reconstruction / demosaic