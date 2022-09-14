%% Intialize
% close all; clear

%% Establish sizes
nPixels = 100;
fieldSizeDegs = 0.5;

% Create stimulus values and full image
stimRVal = 0.1620;  
stimGVal = 0.8461; 
stimBVal = 0.9490;
stimRGB = [stimRVal; stimGVal; stimBVal];
stimulusImageRGB = ones(nPixels, nPixels, 3);
stimulusImageRGB(:, :, 1) = stimRGB(1);
stimulusImageRGB(:, :, 2) = stimRGB(2);
stimulusImageRGB(:, :, 3) = stimRGB(3);

%% Create the display for processing
%
% Ideally would use a conventional/CRT12BitDisplay here instead of mono
displayName = 'conventional';
displayFieldName = 'CRT12BitDisplay';
aoReconDir = getpref('ISETImagePipeline','aoReconDir'); helpDir = '/helperFiles';
theDisplayLoad = load(fullfile(aoReconDir,helpDir,[displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
wlsDisplayOrig = theDisplay.wave;

% Spline underlying display wavelength down to the wavelength sampling we
% will eventually use in the calculations.
wls = (380:10:780)';
theDisplay = displaySet(theDisplay,'wave',wls);
figure; plot(wls,theDisplay.spd);

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
theConeMosaic = ConeResponseCmosaic(2.0, 0, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', true, 'defocusDiopters', 0, ...
        'wave', wls);
% theConeMosaic = ConeResponseCmosaic(2.0, 0, ...
%     'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', true, ...
%     'defocusDiopters', 0);
theOI = theConeMosaic.PSF;
theMosaic = theConeMosaic.Mosaic;

%% Get cone fundamentals
% 
% Get lens transmittance from OI
theOptics = oiGet(theOI,'optics');
theLens = opticsGet(theOptics,'lens');
lensTransmittance = theLens.transmittance';
theMacular = theMosaic.macular;
macTransmittance = theMacular.transmittance';

% Combine QE from mosaic with lens transmittance
coneWls = theMosaic.wave;
coneQENoLens = theMosaic.qe';

% The code below adjusts for size by putting lens transmission into QE.
% We don't think we need to put in the macular pigment, but optionally can.
useMac = false;
if (useMac)
    coneQE = coneQENoLens .* ...
        lensTransmittance(ones(size(coneQENoLens,1),1),:) .* ...
        macTransmittance(ones(size(coneQENoLens,1),1),:);
else
    coneQE = coneQENoLens .* ...
        lensTransmittance(ones(size(coneQENoLens,1),1),:);
end
if (any(coneWls ~= wls))
    error('All our work getting wavelength support consistent failed.');
end
%coneQESpline = SplineCmf(coneWls,coneQE,wls);
coneFundamentals = EnergyToQuanta(wls,coneQE')';

%% Plot cone fundamentals obtained various ways
LconeIndices = find(theMosaic.coneTypes == cMosaic.LCONE_ID);
MconeIndices = find(theMosaic.coneTypes == cMosaic.MCONE_ID);
SconeIndices = find(theMosaic.coneTypes == cMosaic.SCONE_ID);
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
perturbAmount = 0.30;
% perturbAmount2 = perturbAmount1 - (stimDirectExcitations(1) / stimDirectExcitations(2)) + 1;
% 
% c1 = stimDirectExcitations(2) - stimDirectExcitations(1) + (perturbAmount2 * stimDirectExcitations(2));
% c2 = stimDirectExcitations(1) - stimDirectExcitations(2) + (perturbAmount1 * stimDirectExcitations(2));
% 
perturbDirectExcitations = stimDirectExcitations - [0 (perturbAmount * stimDirectExcitations(2)) 0]';
% perturbDirectExcitations = (stimDirectExcitations + [c1 (c2) 0]') .* [1 0.4484 1]';
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
figure; clf;
subplot(1,3,1); hold on;
minVal = 1500;
maxVal = 3000;
origTemp = origMosaicExcitations(:,:,LconeIndices); 
perturbTemp = perturbMosaicExcitations(:,:,LconeIndices);
fracDiff(1) = max(abs(origTemp(:)-perturbTemp(:)))/mean(origMosaicExcitations(:));
slope(:,1) = [origTemp(:) ones(size(origTemp(:)))]\perturbTemp(:);
meanExcitations(1) = mean(origTemp(:));
plot(origTemp(:),perturbTemp(:),'ro','MarkerFaceColor','r','MarkerSize',8);
plot([minVal maxVal],[minVal maxVal],'k');
xlim([minVal maxVal]); ylim([minVal maxVal]);
title("Originally L cones, now M");
axis('square');

subplot(1,3,2); hold on;
minVal = 1500;
maxVal = 3000;
origTemp = origMosaicExcitations(:,:,MconeIndices); 
perturbTemp = perturbMosaicExcitations(:,:,MconeIndices);
fracDiff(2) = max(abs(origTemp(:)-perturbTemp(:)))/mean(origMosaicExcitations(:));
slope(:,2) = [origTemp(:) ones(size(origTemp(:)))]\perturbTemp(:);
meanExcitations(2) = mean(origTemp(:));
plot(origTemp(:),perturbTemp(:),'go','MarkerFaceColor','g','MarkerSize',8);
plot([minVal maxVal],[minVal maxVal],'k');
xlim([minVal maxVal]); ylim([minVal maxVal]);
title("Originally M cones, still M");
axis('square');

subplot(1,3,3); hold on;
minVal = 200;
maxVal = 1000;
origTemp = origMosaicExcitations(:,:,SconeIndices); 
perturbTemp = perturbMosaicExcitations(:,:,SconeIndices);
fracDiff(3) = max(abs(origTemp(:)-perturbTemp(:)))/mean(origMosaicExcitations(:));
slope(:,3) = [origTemp(:) ones(size(origTemp(:)))]\perturbTemp(:);
meanExcitations(3) = mean(origTemp(:));
plot(origTemp(:),perturbTemp(:),'bo','MarkerFaceColor','b','MarkerSize',8);
plot([minVal maxVal],[minVal maxVal],'k');
xlim([minVal maxVal]); ylim([minVal maxVal]);
title("Originally S cones, still S");
axis('square');

%% Compute and report fraction difference
fprintf('Maximum fractional difference in L mosaic excitation = %0.4g\n',fracDiff(1));
fprintf('Maximum fractional difference in M mosaic excitation = %0.4g\n',fracDiff(2));
fprintf('Maximum fractional difference in S mosaic excitation = %0.4g\n',fracDiff(3));
fprintf('Slope for L mosaic excitation = %0.4g, intercept %0.3g\n',slope(1,1),slope(2,1));
fprintf('Slope for M mosaic excitation = %0.4g, intercept %0.3g\n',slope(1,2),slope(2,2));
fprintf('Slope for S mosaic excitation = %0.4g, intercept %0.3g\n',slope(1,3),slope(2,3));
fprintf('Mean excitations %0.4g, %0.4g, %0.4g\n',meanExcitations(1),meanExcitations(2),meanExcitations(3)); 