% t_monoDisplay
%
% Create a display with monochromatic primaries, and demonstrate
% transformation from linear coordinates in a more standard diplsay
% to this one.

% History:
%  10/12/22  Fix primary ordering, relative intensities

%% Clear
clear; close all;

%% Load up a standard display
%displayFile = 'CRT12BitDisplay.mat';
%display = load(fullfile(dataBaseDir, displayFile));
origDisplayFile = 'LCD-Apple';
origDisplay = displayCreate(origDisplayFile);
wls = displayGet(origDisplay,'wave');
S = WlsToS(wls);
if (S(2) > 4)
    error('Want finer wavelength spacing to represent monochromatic primaries');
end

%% Load in cone spectral sensitivities
coneData = load('T_cones_ss2');
T_cones = SplineCmf(coneData.S_cones_ss2,coneData.T_cones_ss2,S);

%% Load in an RGB image
imageSettings = im2double(imread('eagle.jpg'));
[origScene,imageSettingsCheck,imagePrimary] = sceneFromFile(imageSettings, 'rgb', [], origDisplay);
if (max(abs(imageSettings(:)-imageSettingsCheck(:))) ~= 0)
    error('Scene from file not behaving as expected');
end
vcAddObject(origScene);
sceneWindow;

%% The linearized image is what we call a primary image in PTB
[imagePrimaryCal,imageM,imageN] = ImageToCalFormat(imagePrimary);

%% Get the display calibration information in a form DHB understands
%
% The ptb. package routine handles the power unit convention difference
% as it goes between ISETBio to PTB.
gammaMethod = 1;
origCalStruct = ptb.GeneratePTCalStructFromIsetbioDisplayObject(origDisplay);
origCalObj = ObjectToHandleCalOrCalStruct(origCalStruct);
SetSensorColorSpace(origCalObj,T_cones,S);
SetGammaMethod(origCalObj,gammaMethod);
if (any(origCalObj.get('P_ambient') ~= 0))
    error('Below we assume display ambient is zero.');
end

%% Get radiance using PTB knowledge
origPDevice = origCalObj.get('P_device');
origRadiancePTBCal = origPDevice*imagePrimaryCal;
origExcitationsPTBCal = T_cones*origRadiancePTBCal;

%% Get the radiance image out of the ISETBio scene, and into PTB convention.
%
% Verify this works as expected.  Note that we multiply by S(2) when we get
% the energy out of ISETBio, to handle the ISETBio/PTB convention
% difference.
origRadianceFromScence = sceneGet(origScene,'energy')*S(2);
origRadianceFromScenceCal = ImageToCalFormat(origRadianceFromScence);
origExcitationsFromSceneCal = T_cones*origRadianceFromScenceCal;
if (max(abs(origRadiancePTBCal(:)-origRadianceFromScenceCal(:))/mean(origRadianceFromScenceCal(:))) > 1e-6)
    error('Cannot connect PTB and ISETBio worlds the way we should be able to');
end

%% Make the monochromatic primaries
%
% Need a little trial and error to get reasonable relative power for the
% monochromatic primaries.  That's done here by the tweak factors, which
% operate after power of each monochromatic primary has been normalized to
% that of corresponding original primary.
redPowerUW = 6.087;
greenPowerUW = 0.1131;
spotWavelengthsNm = [680 543 440];
spotFWHMsNm = [5 5 5];
primaryTweakFactors = [10 0.7 0.8];
for ww = 1:length(spotWavelengthsNm)
    monoPDeviceRaw(:,ww) = normpdf(wls,spotWavelengthsNm(ww),FWHMToStd(spotFWHMsNm(ww)));
    monoPDevice(:,ww) = primaryTweakFactors(ww)*sum(origPDevice(:,ww))*monoPDeviceRaw(:,ww)/sum(monoPDeviceRaw(:,ww));
end
fprintf('Primary power\n');
for pp = 1:size(monoPDevice,2)
    fprintf('Primary %d, orig: %0.3f; mono: %0.3f\n',pp,sum(origPDevice(:,pp)),sum(monoPDevice(:,pp)));
end

%% Create calibration file for new display
monoCalStruct = origCalStruct;
monoCalStruct.P_device = monoPDevice;
monCalObj = ObjectToHandleCalOrCalStruct(origCalStruct);

%% Create matrix that transforms from original display primaries to mono display primaries
%
% The matrix we want is M_origToMono.  This operates on primaries for the
% original display to produce primaries for a monochromatic display.
M1 = M_PToT(origPDevice,T_cones);
M2 = M_TToP(T_cones,monoPDevice);
M_origToMono = M2*M1;

%% Create primary representation for mono device from that for original device
%
% Note that we divide by S(2) when we push the energy image back into the
% scene, to handle the PTB/ISETBio power convention difference.
monoImagePrimaryCal = M_origToMono*imagePrimaryCal;
monoRadiancePTBCal = monoPDevice*monoImagePrimaryCal;
monoRadiancePTBImage = CalFormatToImage(monoRadiancePTBCal,imageM,imageN);
monoScene = origScene;
monoScene = sceneSet(origScene,'energy',monoRadiancePTBImage/S(2));
monoScene = sceneSet(monoScene,'name','Mono Primaries');
fprintf('Max primaries of original/monochromatic image\n');
for pp = 1:size(imagePrimaryCal,1)
    fprintf('Primary %d, orig: %0.3f; mono: %0.3f\n',pp,max(imagePrimaryCal(pp,:)),max(monoImagePrimaryCal(pp,:)));
end
vcAddObject(monoScene);
sceneWindow;

%% Maybe you want an ISETBio display object for the monochromatic display.
%
% There might be some parameters that don't make it back through, but this
% should work for the spectral ones.
%
% The ptb. package routine handles the power unit convention difference
% as it goes from PTB to ISETBio.
extraCalData = ptb.ExtraCalData;
extraCalData.distance = displayGet(origDisplay,'distance');
monoDisplay = ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct('MonoPrimaries', monoCalStruct, extraCalData, false);
monoDisplay  = rmfield(monoDisplay,'dixel');

%% Visualize some images on display
%
% Use this to fuss with primary intensities so that white point is about
% white, and equal mixture of red and green about yellow.
nPixels = 100;
factor = 1.5;
theImageRGB(:,:,1) = sqrt(factor*0.16)*ones(nPixels,nPixels);
theImageRGB(:,:,2) = sqrt(factor*0.1225)*ones(nPixels,nPixels);
theImageRGB(:,:,3) = sqrt(factor*0.0025)*ones(nPixels,nPixels);
meanLuminanceCdPerM2 = [];
[stimulusScene, ~, stimulusImageLinear] = sceneFromFile(theImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, monoDisplay);
stimulusScene = sceneSet(stimulusScene, 'fov', 0.5);
visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
% save monoDisplay monoDisplay
figure; imshow(theImageRGB);
