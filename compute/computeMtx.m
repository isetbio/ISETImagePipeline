%% FOV mtx with optics
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

fovDegs = 0.20;
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

retina.resetSCone();
retina.reassignSCone(0.1);
retina.PSF = oiCreate('human', pupilSize);

renderFovOptic = retina.forwardRender(imageSize, false, true, false);
save('renderFovOptic.mat', 'renderFovOptic', '-v7.3');
disp('FOV w/ Optics');

%% FOV mtx without optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

renderFovDflmt = retina.forwardRender(imageSize, false, true, false);
save('renderFovDflmt.mat', 'renderFovDflmt', '-v7.3');
disp('FOV without Optics');

%% Peripheral with optics
% Setup a cone mosaic
imageSize = [128, 128, 3];
display = displayCreate('CRT12BitDisplay');
prior   = load('sparsePrior.mat');

% Generate cone mosaic - [eccX, eccY] deg ecc
eccX = 18.0; eccY = 18.0;
retina = ConeResponseCmosaic...
    (eccX, eccY, 'fovealDegree', 1.0, 'pupilSize', 3.0, 'subjectID', 6);

renderPeriOptic = retina.forwardRender(imageSize, false, false);
save('renderPeriOptic.mat', 'renderPeriOptic', '-v7.3');
disp('Peri w/ Optics');

%% Turn off optics
% Use a large pupil size to reduce the effect of diffraction
diffPupil = 10.0;
retina.PSF = ConeResponse.psfDiffLmt(diffPupil);

renderPeriDflmt = retina.forwardRender(imageSize, false, false);
save('renderPeriDflmt.mat', 'renderPeriDflmt', '-v7.3');
disp('Peri without Optics');

%% Create cone mosaic
fovDegs = 1.0;
pupilSize = 3.0;
retina = ConeResponse('eccBasedConeDensity', false, 'eccBasedConeQuantal', false, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

retina.resetSCone();
retina.reassignSCone(0.1);

%% Change all M cone to L cone
retina.reassignCone(0.0, retina.M_Cone_Idx, retina.L_Cone_Idx, false);

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

renderDeuteranopia = retina.forwardRender(imageSize, false, true, false);
save('renderDeuteranopia.mat', 'renderDeuteranopia', '-v7.3');
disp('Deuteranopia');

%% Change all L cone to M cone
retina.reassignCone(0.0, retina.L_Cone_Idx, retina.M_Cone_Idx, false);

[L, M, S] = retina.coneCount();
fprintf('Number of cones: L - %d, M - %d, S - %d \n', L, M, S);

renderProtanopia = retina.forwardRender(imageSize, false, true, false);
save('renderProtanopia.mat', 'renderProtanopia', '-v7.3');
disp('Protanopia');

%% Dichromacy without S cone
retina = ConeResponse('eccBasedConeDensity', false, 'eccBasedConeQuantal', false, ...
    'fovealDegree', fovDegs, 'display', display, 'pupilSize', pupilSize, ...
    'integrationTime', 1.0);

retina.resetSCone();
retina.reassignSCone(0.0);

renderTritanopia = retina.forwardRender(imageSize, false, true, false);
save('renderTritanopia.mat', 'renderTritanopia', '-v7.3');
disp('Tritanopia');
