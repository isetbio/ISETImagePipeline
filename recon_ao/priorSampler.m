%% Clear
clear; close all;
display = displayCreate('CRT12BitDisplay');

%% Sample images from the sparse prior
prior = load('sparsePrior.mat');

imSize = [100, 100, 3];
image = sparseSampler(prior, imSize);

figure();
imshow(gammaCorrection(image, display), ...
       'InitialMagnification', 1000);

%% Sample images from a 1/f specturm prior
imSize = [100, 100, 3];
image = spectrumSampler(imSize);

figure();
imshow(gammaCorrection(image, display), ...
       'InitialMagnification', 1000);

%% Sample from specturm prior + nonlinear transformation
imSize = [100, 100, 3];

gain = 20; 
image = spectrumSampler(imSize, gain);

figure();
imshow(gammaCorrection(image, display), ...
       'InitialMagnification', 1000);
