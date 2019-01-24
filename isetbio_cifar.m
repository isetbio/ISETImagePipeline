%% Create scene 
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.57);
meanLuminanceCdPerM2 = 100;

% Load picture
realizedStimulusScene = sceneFromFile(image, 'rgb', ...
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
