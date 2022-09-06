%% Intialize
close all; clear all

%% Establish sizes
nPixels = 100;
fieldSizeDegs = 0.5;

% Create stimulus values and full image
stimRVal = 0.80;  
stimGVal = 0.65; 
stimBVal = 0.10;
stimRGB = [stimRVal; stimGVal; stimBVal];
stimulusImageRGB = ones(nPixels, nPixels, 3);
stimulusImageRGB(:, :, 1) = stimRGB(1);
stimulusImageRGB(:, :, 2) = stimRGB(2);
stimulusImageRGB(:, :, 3) = stimRGB(3);

%% Create the display for processing
%
% Ideally would use a conventional/CRT12BitDisplay here instead of mono
displayName = 'mono';
displayFieldName = 'monoDisplay';
aoReconDir = getpref('ISETImagePipeline','aoReconDir'); helpDir = '/helperFiles';
theDisplayLoad = load(fullfile(aoReconDir,helpDir,[displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
wls = theDisplay.wave;

%% Adjustments if using mono displayName
% 
% Overwrite Gamma portion that was included in the building of mono cone
% mosaics when using a mono display
displayGammaBits = 12;
displayGammaGamma = 2;
if strcmp(displayName, 'mono')
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end

%% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);

%% Explicitly do gamma correction and check below.  
% 
% We'll need to be able to do this lower down.
stimulusImageRGB1 = gammaCorrection(stimulusImageLinear, theDisplay);

% And because it's here, do inverse gamma correction the ISETBio way.
% Note that currently there is an invGammaCorrection function on the
% path that actually does the forward gammaCorrection.
gammaLength = size(theDisplay.gamma,1);
stimulusImageLinear1 = ieLUTDigital(round((gammaLength-1)*stimulusImageRGB1), theDisplay.gamma);

% Check that things we think should match actually do.  Don't expect exact
% agreement because gamma table is discrete. 
% 
% This also pulls out a single pixel's linear values, which we will need
% when we compute metamers below
centerPixel = round(nPixels/2);
examplePixelStimulusRGB = squeeze(stimulusImageRGB(centerPixel,centerPixel,:));
if (any(stimRGB ~= examplePixelStimulusRGB))
    error('Stimulus not built the way we expect');
end
examplePixelStimulusRGB1 = squeeze(stimulusImageRGB1(centerPixel,centerPixel,:));
if (max(abs(examplePixelStimulusRGB1-stimRGB))/mean(stimRGB) > 1e-2)
    error('Explict gamma correction from linear values doesn''t match RGB input');
end
stimLinear = squeeze(stimulusImageLinear(centerPixel,centerPixel,:));
stimLinear1 = squeeze(stimulusImageLinear1(centerPixel,centerPixel,:));
if (max(abs(stimLinear-stimLinear1))/mean(stimLinear) > 1e-3)
    error('Explict gamma correction from linear values doesn''t match RGB input');
end

%% Build a mosaic object and pull out parts we need.
theConeMosaic = ConeResponseCmosaic(0, 0, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', true);
theOI = theConeMosaic.PSF;
theMosaic = theConeMosaic.Mosaic;

%% Get cone fundamentals
% 
% Get lens transmittance from OI
theOptics = oiGet(theOI,'optics');
theLens = opticsGet(theOptics,'lens');
lensTransmittance = theLens.transmittance';

% Combine QE from mosaic with lens transmittance
coneWls = theMosaic.wave;
coneQENoLens = theMosaic.qe';
% The code below adjusts for size by putting lens transmission into 3 x 31
coneQE = coneQENoLens .* ...
    lensTransmittance(ones(size(coneQENoLens,1),1),:);
coneQESpline = SplineCmf(coneWls,coneQE,wls);
coneFundamentals = EnergyToQuanta(wls,coneQESpline')';
figure; clf; hold on; plot(wls,coneFundamentals','r');

%% Make all the mosaic M cones into L cones
coneTypes = theMosaic.coneTypes;
coneIndex = find(coneTypes == cMosaic.MCONE_ID);
theMosaic.reassignTypeOfCones(coneIndex, cMosaic.LCONE_ID);
coneTypes = theMosaic.coneTypes;
coneIndex = find(coneTypes == cMosaic.MCONE_ID);
if (~isempty(coneIndex))
    error('Did not actually change mosaic cone types');
end

%% Compute cone mosaic responses to original image 
% 
% Checking the OI computation since the documentation says this should be
% flipped
% theOI = oiCompute(theOI,stimulusScene);
theOI = oiCompute(stimulusScene, theOI);
origMosaicExcitations = theMosaic.compute(theOI);

%% Compute cone excitations directly from linear stimulus RGB values
%
% These are not scaled right, but we don't care about that as
% all we need is the stimulus direction that is silent to a cone class.
B_primary = theDisplay.spd;
stimDirectExcitations = (coneFundamentals*B_primary)*stimLinear;

%% Perturb M cone component of the directl computed cone excitations
perturbAmount = 0.15;
perturbDirectExcitations = stimDirectExcitations + [0 perturbAmount*stimDirectExcitations(2) 0]';
perturbDirectLinear = inv(coneFundamentals*B_primary)*perturbDirectExcitations;
perturbRGB = gammaCorrection(perturbDirectLinear, theDisplay);
perturbImageRGB = ones(nPixels, nPixels, 3);
perturbImageRGB(:, :, 1) = perturbRGB(1);
perturbImageRGB(:, :, 2) = perturbRGB(2);
perturbImageRGB(:, :, 3) = perturbRGB(3);

%% Create perturbed scene and mosaic excitations
[perturbScene, ~, perturbImageLinear] = sceneFromFile(perturbImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
perturbScene = sceneSet(perturbScene, 'fov', fieldSizeDegs);
perturbOI = oiCompute(theOI,perturbScene);
perturbMosaicExcitations = theMosaic.compute(perturbOI);

%% Figure compares mosaic excitations for the two different stimuli
figure; clf; hold on;
maxVal = 1500;
plot(origMosaicExcitations(:),perturbMosaicExcitations(:),'ro','MarkerFaceColor','r','MarkerSize',8);
plot([0 maxVal],[0 maxVal],'k');
xlim([0 maxVal]); ylim([0 maxVal]);
axis('square');

%% Compute and report fraction difference
fracDiff = max(abs(origMosaicExcitations(:)-perturbMosaicExcitations(:)))/mean(origMosaicExcitations(:));
fprintf('Maximum fractional difference in mosaic excitation = %0.4g\n',fracDiff);
