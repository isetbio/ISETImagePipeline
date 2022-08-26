
nPixels = 100;
fieldSizeDegs = 1;

stimBgVal = 0.1;
stimRVal = [0.80] + stimBgVal;  
stimGVal = [0.65] + stimBgVal; 
stimBVal = [0.10] + stimBgVal;
stimRGB = [stimRVal; stimGVal; stimBVal];

stimulusImageRGB = ones(nPixels, nPixels, 3) * stimBgVal;
stimulusImageRGB(:, :, 1) = stimRGB(1);
stimulusImageRGB(:, :, 2) = stimRGB(2);
stimulusImageRGB(:, :, 3) = stimRGB(3);

displayName = 'conventional';
displayFieldName = 'CRT12BitDisplay';
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
theDisplayLoad = load(fullfile(aoReconDir,[displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);


% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);

% Grab the primaries of the display
B_primary = theDisplay.spd;
wls = theDisplay.wave;

% Build a mosaic
theConeMosaic = ConeResponseCmosaic(0, 0, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', true);
theMosaic = theConeMosaic.Mosaic;
coneWls = theMosaic.wave;
coneQE = theMosaic.qe';
coneQESpline = SplineCmf(coneWls,coneQE,wls);
coneFundamentals = EnergyToQuanta(wls,coneQESpline')';
figure; clf; hold on; plot(wls,coneFundamentals','r');

% Compute cone excitations from RGB values
coneExcitations = (coneFundamentals*B_primary)*stimRGB;

% Perturb M cones
coneExcitations1 = coneExcitations + [0 0.05*coneExcitations(2) 0]';
stimRGB1 = inv(coneFundamentals*B_primary)*coneExcitations1;
