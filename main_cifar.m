% Example shows how to use ISETPipeineToolbox/ConeResponse
% Create cone mosaic with default (simple) parameter and human psf (optics)
retina = ConeResponse();

%% cone mosaic
retina.visualizeCone()

%% Response to image
imgIdx = 1;
load('./image_all.mat');
image = reshape(image_all(imgIdx, :, :), [32, 32, 3]);

% Compute response
[excitation, oi, l, m, s] = retina.compute(image);

%% Visualize OI
retina.visualizeOI();

%% Visualize cone excitation
retina.visualizeExcitation();

%% previous code
fovDeg = 1;
theMosaic = coneMosaicHex(5, ...                % hex lattice sampling factor
   'fovDegs', fovDeg, ...                       % match mosaic width to stimulus size
   'eccBasedConeDensity', false, ...            % cone density varies with eccentricity
   'eccBasedConeQuantalEfficiency', false, ...  % cone quantal efficiency varies with eccentricity
   'integrationTime', 0.1, ...                  % 0.1s integration time
   'maxGridAdjustmentIterations', 50);          % terminate iterative lattice adjustment after 50 iterations

theMosaic.visualizeGrid(...
    'backgroundColor', [1 1 1], ...
    'ticksInVisualDegs', true);

% Poisson noise model, mean response
theMosaic.noiseFlag = 'none';

%% Compute response
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.57);
psf = oiCreate('wvf human');

%% Response to image
load('./image_all.mat');
imgIdx = 1;
image = reshape(image_all(imgIdx, :, :), [32, 32, 3]);
[excitation, io, L, M, S] = computeResponse(presentationDisplay, fovDeg, psf, theMosaic, image);

%% Visualization
visualizeOpticalImage(io, 'displayRetinalContrastProfiles', true);
visualizeConeMosaicResponses(theMosaic, excitation, 'R*/cone/tau');
