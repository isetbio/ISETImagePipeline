close all; clear all
% Establish sizing framework
nPixels = 100;
fieldSizeDegs = 1;

% Create stimulus values
stimBgVal = 0.0;
stimRVal = [0.80] + stimBgVal;  
stimGVal = [0.65] + stimBgVal; 
stimBVal = [0.10] + stimBgVal;
stimRGB = [stimRVal; stimGVal; stimBVal];

% Apply values to full image
stimulusImageRGB = ones(nPixels, nPixels, 3) * stimBgVal;
stimulusImageRGB(:, :, 1) = stimRGB(1);
stimulusImageRGB(:, :, 2) = stimRGB(2);
stimulusImageRGB(:, :, 3) = stimRGB(3);

% Create the display for processing
displayName = 'mono';
displayFieldName = 'monoDisplay';
aoReconDir = getpref('ISETImagePipeline','aoReconDir'); helpDir = '/helperFiles';
theDisplayLoad = load(fullfile(aoReconDir,helpDir,[displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);

% Overwrite Gamma portion that was included in the building of mono cone
% mosaics, must adjust accordingly if using a 'conventional' displayName 
displayGammaBits = 12;
displayGammaGamma = 2;
gammaInput = linspace(0,1,2^displayGammaBits)';
gammaOutput = gammaInput.^displayGammaGamma;
theDisplay.gamma = gammaOutput(:,[1 1 1]);


% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);

% Explicitly do gamma correction and check below.  We'll need to be able to
% do this lower down.
stimulusImageRGB1 = gammaCorrection(stimulusImageLinear, theDisplay);

% And becuase it's here, do inverse gamma correction the ISETBio way.
% Note that currently there is an invGammaCorrection function on the
% path that actually does the forward gammaCorrection.
gammaLength = size(theDisplay.gamma,1);
stimulusImageLinear1 = ieLUTDigital(round((gammaLength-1)*stimulusImageRGB1), theDisplay.gamma);

% Check that things we think should match actually do.  Don't expect exact
% agreement because gamma table is discrete.
examplePixelStimulusRGB = squeeze(stimulusImageRGB(50,50,:));
if (any(stimRGB ~= examplePixelStimulusRGB))
    error('Stimulus not built the way we expect');
end
examplePixelStimulusRGB1 = squeeze(stimulusImageRGB1(50,50,:));
if (max(abs(examplePixelStimulusRGB1-stimRGB))/mean(stimRGB) > 1e-2)
    error('Explict gamma correction from linear values doesn''t match RGB input');
end
examplePixelLinear = squeeze(stimulusImageLinear(50,50,:));
examplePixelLinear1 = squeeze(stimulusImageLinear1(50,50,:));
if (max(abs(examplePixelLinear-examplePixelLinear1))/mean(examplePixelLinear) > 1e-3)
    error('Explict gamma correction from linear values doesn''t match RGB input');
end

% Grab the primaries of the display
B_primary = theDisplay.spd;
wls = theDisplay.wave;

% Build a mosaic object and pull out parts we need.
theConeMosaic = ConeResponseCmosaic(0, 0, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', true);

% Get lens transmittance.  Should really to this using a get
theOI = theConeMosaic.PSF;
theOptics = oiGet(theOI,'optics');
theLens = opticsGet(theOptics,'lens');
lensTransmittance = theLens.transmittance';

theMosaic = theConeMosaic.Mosaic;
coneWls = theMosaic.wave;
coneQENoLens = theMosaic.qe';
coneQE = coneQENoLens .* lensTransmittance(ones(size(coneQENoLens,1),1),:);
coneQESpline = SplineCmf(coneWls,coneQE,wls);
coneFundamentals = EnergyToQuanta(wls,coneQESpline')';
figure; clf; hold on; plot(wls,coneFundamentals','r');

% Compute cone excitations from RGB values
coneExcitations = (coneFundamentals*B_primary)*stimRGB;

% Perturb M cones
coneExcitations1 = coneExcitations + [0 0.25*coneExcitations(2) 0]';
stimPerturbMLinear = inv(coneFundamentals*B_primary)*coneExcitations1;
stimPerturbMRGB = gammaCorrection(stimPerturbMLinear, theDisplay);
