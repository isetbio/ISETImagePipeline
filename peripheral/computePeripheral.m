%% Load display setup & constant
dataBaseDir = getpref('ISETImagePipeline', 'dataDir');
display = load(fullfile(dataBaseDir, 'CRT12BitDisplay.mat'));

imageSize = [128, 128, 3];

thisImageSet = 'ILSVRC';
imageName    = 'ILSVRC2017_test_00000021.JPEG';

fileDir = fullfile(dataBaseDir, thisImageSet, imageName);
image   = imresize(im2double(imread(fileDir)), 0.4);

patch = sampleImage(image, imageSize(1));
patch(patch < 0) = 0;
patch(patch > 1) = 1;
imshow(patch, 'InitialMagnification', 400);

%% Generate cone mosaic - Foveal
eccX = 0.5; eccY = 0;
retinaFov = ConeResponsePeripheral(eccX, eccY, 'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

%% Visualization
[~, ~, linear, coneVec] = retinaFov.compute(patch);

retinaFov.visualizeOI();
retinaFov.visualizeMosaic();
retinaFov.visualizeExcitation();

%% Generate cone mosaic - 5 deg ecc
eccX = 5; eccY = 0;
retina5 = ConeResponsePeripheral(eccX, eccY, 'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

%% Visualization
[~, ~, linear, coneVec] = retina5.compute(patch);

retina5.visualizeOI();
retina5.visualizeMosaic();
retina5.visualizeExcitation();

%% Generate cone mosaic - 10 deg ecc
eccX = 10; eccY = 0;
retina10 = ConeResponsePeripheral(eccX, eccY, 'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'pupilSize', 2.0);

%% Visualization
[~, ~, linear, coneVec] = retina10.compute(patch);

retina10.visualizeOI();
retina10.visualizeMosaic();
retina10.visualizeExcitation();

%% Reconstruction Computation - render matrix
renderFov = retinaFov.forwardRender(imageSize);
render5   = retina5.forwardRender(imageSize);
render10  = retina10.forwardRender(imageSize); 