%% Show how to manipulate different components
% that influence and contribute to cone sensitivity function

% Description:
%    In the first part of the tutorial, 
%    we show that cone can get cone fundamentals from the product of:
%     - Lens Transmittance
%     - Macular Transmittance
%     - Photopigment Absorbance
% 
%    In the second part of the tutorial,
%    we get rid of the differential effect of the lens and
%    macular transmittance on short wavelength light

%% Figure format
try 
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
        'figureWidthInches', 15, ...
        'figureHeightInches', 10);
catch EXP
    fprintf('plotlab not available, use default MATLAB style \n');
end

%% Create a uniform scene (with equal photon count across the entire spectrum)
uniformScene = sceneCreate('uniform equal photon', 128);
uniformScene = sceneSet(uniformScene, 'wAngular', 1.0);
uniformScene = sceneSet(uniformScene, 'distance', 0.57);

%% Generate a cone mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.4, 'display', displayCreate('CRT12BitDisplay'), 'pupilSize', 2.0, ...
    'integrationTime', 1.0);

retina.reassignSCone(0.25);
% retina.visualizeMosaic();

%% Plot cone fundamental
transPlot = figure();
subplot(1, 2, 1);

% get lens transmittance
lens = oiGet(retina.PSF, 'lens');
lensTrans = lens.transmittance;
plot(lens.wave, lensTrans, 'LineWidth', 2); hold on;

% get macular transmittance
macular = retina.Mosaic.macular;
macuTrans = macular.transmittance;
plot(macular.wave, macuTrans, '--', 'LineWidth', 2);
legend('lens', 'macular');

xlabel('Wavelength');
ylabel('Transmittance');
grid off; box off;

% get photopigment object
pigment = retina.Mosaic.pigment;

conePlot = figure();
subplot(1, 2, 1); hold on;
colorSeq = ['r', 'g', 'b'];
for idx = 1 : length(colorSeq)
    % plot absorbance against wavelength
    plot(pigment.wave, log(pigment.absorbance(:, idx) .* lensTrans .* macuTrans), ...
        'LineWidth', 2, 'Color', colorSeq(idx));
end

xlabel('Wavelength');
ylabel('Absorption'); ylim([-15, 0]);
grid off; box off;

%% Change the optical density of the lens
optics = retina.PSF;
lens0 = oiGet(optics, 'lens');
wls = lens0.wave;

lensUnitDensity1 = lens0.unitDensity;
lensUnitDensity1 = zeros(size(lensUnitDensity1));
lensPeakDensity1 = 1;

lens1 = Lens('wave', wls, ...
    'unitDensity', lensUnitDensity1, 'density', lensPeakDensity1);
retina.PSF = oiSet(optics, 'lens', lens1);

%% Change the macular density
macular0 = retina.Mosaic.macular;
wls = macular0.wave;

macularUnitDensity1 = zeros(size(macular0.unitDensity));
macularDensity1 = macular0.density;

macular1 = Macular('wave', wls, 'unitDensity', macularUnitDensity1, 'density', macularDensity1);
retina.Mosaic.macular = macular1;

%% Plot Cone fundamental after changing the optical and macular density
figure(transPlot);
subplot(1, 2, 2);

% get lens transmittance
lens = oiGet(retina.PSF, 'lens');
lensTrans = lens.transmittance;
plot(lens.wave, lensTrans, 'LineWidth', 2); hold on;

% get macular transmittance
macular = retina.Mosaic.macular;
macuTrans = macular.transmittance;
plot(macular.wave, macuTrans, '--', 'LineWidth', 2);
legend('lens', 'macular');

xlabel('Wavelength');
ylabel('Transmittance');
grid off; box off;

% get photopigment object
pigment = retina.Mosaic.pigment;

figure(conePlot);
subplot(1, 2, 2); hold on;
colorSeq = ['r', 'g', 'b'];
for idx = 1 : length(colorSeq)
    % plot absorbance against wavelength
    plot(pigment.wave, log(pigment.absorbance(:, idx) .* lensTrans .* macuTrans), ...
        'LineWidth', 2, 'Color', colorSeq(idx));
end

xlabel('Wavelength');
ylabel('Absorption'); ylim([-15, 0]);
grid off; box off;

%% Compute cone isomerization with uniform scene
[~, ~, L, M, S] = retina.computeWithScene(uniformScene);

% retina.visualizeExcitation();
% visualizeOpticalImage(retina.LastOI, 'displayRadianceMaps', true, ...
%                 'displayRetinalContrastProfiles', true);
