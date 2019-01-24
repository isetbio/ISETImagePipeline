%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
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
