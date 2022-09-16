%% Intialize
close all; clear

%% Parameters
% 
% Establish sizes
nPixels = 100;
fieldSizeDegs = 0.5;

% Randomize intensities across image?
%
% Seemed clever at one point to better constrain the
% regression, but chromatic aberation
% messes up the metamerism for an image with spatial
% structure so isn't all that useful in the end.
randomizeImageIntensities = false;

% Force all mosaic params to foveal values
forceFovealValues = true;

% Use cone fundamentals obtained through simulated measurement?
useFundamentalsBySimulation = true;

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
% Use conventional/CRT12BitDisplay
displayName = 'conventional';
displayFieldName = 'CRT12BitDisplay';
aoReconDir = getpref('ISETImagePipeline','aoReconDir'); helpDir = '/helperFiles';
theDisplayLoad = load(fullfile(aoReconDir,helpDir,[displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
wlsDisplayOrig = theDisplay.wave;

% Spline underlying display wavelength down to the wavelength sampling we
% will eventually use in the calculations.
wls = (380:5:780)';
theDisplay = displaySet(theDisplay,'wave',wls);

%% Adjustments if using mono displayName
% 
% Overwrite gamma function that was included in the building of mono cone
% mosaics when using a mono display
displayGammaBits = 12;
displayGammaGamma = 2;
if strcmp(displayName, 'mono')
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end

% Make sure display has no ambient
if (any(theDisplay.ambient ~= 0))
    error('Code assumes display ambient is zero, but it is not');
end

%% Show the stimulus by creating an ISETBio scene
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);

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

%% Randomize image intensities?
%
% This scales all of the pixels in the uniform field by a random
% number between 0 and 1, in the linear image domain. The reason
% for implementing this was to better constrain the regression comparing
% excitations for metameric pairs, since any scaling of a metameric pair
% is another metameric pair.  But this approach hits the cold hard reality
% of chromatic aberration - since this produces spatial structure in the
% image, aberration will break the metamerism.  Leaving the code here in
% case we ever want to explore this effect in more detail.
%
% Note that when this is true, we rewrite stimulusScene and
% stimulusImageLinear so they no longer match their values above.
%
% Also note that we'll use the same randVal below to apply the same random
% scaling to the perturbed scene/image. So don't rewrite that variable
% between here and below.
%
% Finally note that to preserve metamerism, the scaling must be applied to
% the linear image, not the gamma corrected RGB values.
if (randomizeImageIntensities)
    randVals = rand(size(squeeze(stimulusImageLinear1(:,:,1))));
    for cc = 1:size(stimulusImageLinear,3)
        stimulusImageLinear(:,:,cc) = stimulusImageLinear(:,:,cc) .* randVals;
    end
    stimulusImageRGB = gammaCorrection(stimulusImageLinear, theDisplay);

    [stimulusScene, ~, stimulusImageLinear] = sceneFromFile(stimulusImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
    stimulusScene = sceneSet(stimulusScene, 'fov', fieldSizeDegs);
end

%% Visualize scene
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);

%% Build a mosaic object and pull out parts we need.
% 
% For this excercise, we want a simple mosaic so we take explicit
% control of all the fancy eccentricity varying parameters.
theConeMosaic = ConeResponseCmosaic(2.0, 0, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', 3, 'useRandomSeed', false, 'defocusDiopters', 0, ...
        'wave', wls, ...
        'rodIntrusionAdjustedConeAperture', false, ...
        'eccVaryingConeAperture', false, ...
        'eccVaryingConeBlur', false, ...
        'eccVaryingOuterSegmentLength', false, ...
        'eccVaryingMacularPigmentDensity', false, ...
        'eccVaryingMacularPigmentDensityDynamic', false, ...
        'anchorAllEccVaryingParamsToTheirFovealValues', forceFovealValues);
theOI = theConeMosaic.PSF;
theMosaic = theConeMosaic.Mosaic;

%% Get cone fundamentals out of optics and mosaic
% 
% Get lens transmittance from OI
theOptics = oiGet(theOI,'optics');
theLens = opticsGet(theOptics,'lens');
lensTransmittance = theLens.transmittance';

% Combine QE from mosaic with lens transmittance
coneWls = theMosaic.wave;
coneQENoLens = theMosaic.qe';

% The code below adjusts for size by putting lens transmission into 3 x 31
coneQEFromObjects = coneQENoLens .* ...
    lensTransmittance(ones(size(coneQENoLens,1),1),:);
if (any(coneWls ~= wls))
    error('All our work getting wavelength support consistent failed.');
end
coneFundamentalsFromObjects = EnergyToQuanta(wls,coneQEFromObjects')';

%% Another way to get the fundamentals
%
% Build a set of monochromatic images at each sample wavelength.  The photon
% level is scene radiance in photons/sr-m2-nm-sec.
photonsPerSrM2NMSec = 1e25;
for ww = 1:length(wls)
    % Set up a dummy scene. Spatially uniform with a black body
    % spectrum of 5000 degK.  We'll replace the scene contents just below.
    scene{ww} = sceneCreate('uniform bb',nPixels,5000,wls);

    % Use small field of view to minimize effects of eccentricity, and also
    % so we don't need too many pixels (for efficiency in this demo).
    scene{ww} = sceneSet(scene{ww},'fov',fieldSizeDegs);

    % Get photons and rewrite to be monochromatic constant power in
    % photons/sec-nm.
    photons = sceneGet(scene{ww},'photons');
    photons = zeros(size(photons));
    photons(:,:,ww) = photonsPerSrM2NMSec*ones(nPixels,nPixels);
    scene{ww} = sceneSet(scene{ww},'photons',photons);
end

% Compute the retinal image and cone excitations for each scene
%
% Use this to get the cone fundamental that ISETBio is effectively
% using.
LConeIndices = find(theMosaic.coneTypes == cMosaic.LCONE_ID);
MConeIndices = find(theMosaic.coneTypes == cMosaic.MCONE_ID);
SConeIndices = find(theMosaic.coneTypes == cMosaic.SCONE_ID);
coneQEBySimulation = zeros(size(coneQEFromObjects));
fprintf('Finding ')
for ww = 1:length(wls)
    % Compute retinal image
    oiComputed{ww} = oiCompute(scene{ww}, theOI);

    % Compute noise free cone excitations
    coneExcitationsToMono{ww} = theMosaic.compute(oiComputed{ww},'nTrials', 1);

    % Find L, M and S cone excitations at this wavelength.  This is
    % accomplished by extracting the mean response of each cone type from
    % the excitations computed just above and diviting by the input power
    % in the scene.
    %
    % We expect the answer to be proportional to the fundamentals we
    % computed above, because we have not (yet) accounted for the geometry
    % between radiance in the scene and retinal irradiance, nor the cone
    % integration time.
    excitationsConeISETBioMono(1,ww) = mean(coneExcitationsToMono{ww}(LConeIndices));
    excitationsConeISETBioMono(2,ww) = mean(coneExcitationsToMono{ww}(MConeIndices));
    excitationsConeISETBioMono(3,ww) = mean(coneExcitationsToMono{ww}(SConeIndices));
    coneQEBySimulation(1,ww) = excitationsConeISETBioMono(1,ww)/photonsPerSrM2NMSec;
    coneQEBySimulation(2,ww) = excitationsConeISETBioMono(2,ww)/photonsPerSrM2NMSec;
    coneQEBySimulation(3,ww) = excitationsConeISETBioMono(3,ww)/photonsPerSrM2NMSec;
end

% Match scale to the coneQE obtained directly
%
% Do this jointly for the three cone classes.  We want to keep a consistent
% relative scale.
coneScaleFactor = coneQEBySimulation(:)\coneQEFromObjects(:);
for cc = 1:size(coneQEFromObjects,1)
    coneQEBySimulation(cc,:) = coneScaleFactor*coneQEBySimulation(cc,:);
end
coneFundamentalsBySimulation = EnergyToQuanta(wls,coneQEBySimulation')';

%% Choose which fundamentals
if (useFundamentalsBySimulation)
    coneFundamentals = coneFundamentalsBySimulation;
else
    coneFundamentals = coneFundamentalsFromObjects;
end

%% Plot cone fundamentals obtained various ways
figure; clf; hold on;
plot(wls,coneFundamentalsFromObjects','r','LineWidth',4);
plot(wls,coneFundamentalsBySimulation','y-','LineWidth',2);
xlabel('Wavelength'); ylabel('Fundamental');
legend({'From Object', 'By Simulation'})

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
perturbAmount = 0.20;
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
if (randomizeImageIntensities)
    for cc = 1:size(perturbImageLinear,3)
        perturbImageLinear(:,:,cc) = perturbImageLinear(:,:,cc) .* randVals;
    end
    perturbImageRGB = gammaCorrection(perturbImageLinear, theDisplay);

    [perturbScene, ~, perturbImageLinear] = sceneFromFile(perturbImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
    perturbScene = sceneSet(perturbScene, 'fov', fieldSizeDegs);
end
visualizeScene(perturbScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
perturbOI = oiCompute(theOI,perturbScene);
perturbMosaicExcitations = theMosaic.compute(perturbOI);

%% Figure compares mosaic excitations for the two different stimuli
figure; clf;
subplot(1,3,1); hold on;
origTemp = origMosaicExcitations(:,:,LConeIndices); 
perturbTemp = perturbMosaicExcitations(:,:,LConeIndices);
fracDiff(1) = max(abs(origTemp(:)-perturbTemp(:)))/mean(origMosaicExcitations(:));
slope(:,1) = [origTemp(:) ones(size(origTemp(:)))]\perturbTemp(:);
meanExcitation(1) = mean(origTemp(:));
relStddevExcitation(1,1) = std(origTemp(:))/meanExcitation(1);
relStddevExcitation(1,2) = std(perturbTemp(:))/meanExcitation(1);
minVal = 0.9*min(origTemp(:));
maxVal = 1.1*max(origTemp(:));
plot(origTemp(:),perturbTemp(:),'ro','MarkerFaceColor','r','MarkerSize',8);
plot([minVal maxVal],[minVal maxVal],'k');
xlim([minVal maxVal]); ylim([minVal maxVal]);
title("Originally L cones");
axis('square');

subplot(1,3,2); hold on;
origTemp = origMosaicExcitations(:,:,MConeIndices); 
perturbTemp = perturbMosaicExcitations(:,:,MConeIndices);
fracDiff(2) = max(abs(origTemp(:)-perturbTemp(:)))/mean(origMosaicExcitations(:));
slope(:,2) = [origTemp(:) ones(size(origTemp(:)))]\perturbTemp(:);
meanExcitation(2) = mean(origTemp(:));
relStddevExcitation(2,1) = std(origTemp(:))/meanExcitation(2);
relStddevExcitation(2,2) = std(perturbTemp(:))/meanExcitation(2);
minVal = 0.9*min(origTemp(:));
maxVal = 1.1*max(origTemp(:));
plot(origTemp(:),perturbTemp(:),'go','MarkerFaceColor','g','MarkerSize',8);
plot([minVal maxVal],[minVal maxVal],'k');
xlim([minVal maxVal]); ylim([minVal maxVal]);
title("Originally M cones");
axis('square');

subplot(1,3,3); hold on;
origTemp = origMosaicExcitations(:,:,SConeIndices); 
perturbTemp = perturbMosaicExcitations(:,:,SConeIndices);
fracDiff(3) = max(abs(origTemp(:)-perturbTemp(:)))/mean(origMosaicExcitations(:));
slope(:,3) = [origTemp(:) ones(size(origTemp(:)))]\perturbTemp(:);
meanExcitation(3) = mean(origTemp(:));
relStddevExcitation(3,1) = std(origTemp(:))/meanExcitation(3);
relStddevExcitation(3,2) = std(perturbTemp(:))/meanExcitation(3);
minVal = 0.9*min(origTemp(:));
maxVal = 1.1*max(origTemp(:));
plot(origTemp(:),perturbTemp(:),'bo','MarkerFaceColor','b','MarkerSize',8);
plot([minVal maxVal],[minVal maxVal],'k');
xlim([minVal maxVal]); ylim([minVal maxVal]);
title("Originally S cones");
axis('square');

%% Compute and report fraction difference
fprintf('Maximum fractional difference in L mosaic excitation = %0.4g\n',fracDiff(1));
fprintf('Maximum fractional difference in M mosaic excitation = %0.4g\n',fracDiff(2));
fprintf('Maximum fractional difference in S mosaic excitation = %0.4g\n',fracDiff(3));
fprintf('Relative stdevs of L mosaic excitations = %0.4g, %0.4g\n',relStddevExcitation(1,1),relStddevExcitation(2,1));
fprintf('Relative stdevs of M mosaic excitations = %0.4g, %0.4g\n',relStddevExcitation(2,1),relStddevExcitation(2,2));
fprintf('Relative stdevs of M mosaic excitations = %0.4g, %0.4g\n',relStddevExcitation(3,1),relStddevExcitation(3,2));
fprintf('Slope for L mosaic excitation = %0.4g, intercept = %0.4g\n',slope(1,1),slope(2,1));
fprintf('Slope for M mosaic excitation = %0.4g, intercept = %0.4g\n',slope(1,2),slope(2,2));
fprintf('Slope for S mosaic excitation = %0.4g, intercept = %0.4g\n',slope(1,3),slope(2,3));
fprintf('Mean L %0.4g\n',meanExcitation(1));
fprintf('Mean M %0.4g\n',meanExcitation(2));
fprintf('Mean S %0.4g\n',meanExcitation(3));
