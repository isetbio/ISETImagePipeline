% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

% M1
load('retina_render_M1.mat');

regPara = 5e-4;
outputArray_M1 = reconArray(inputLinear, renderArray, regBasis, MU, regPara, imageSize);

% M2
load('retina_render_M2.mat');

regPara = 5e-4;
outputArray_M2 = reconArray(inputLinear, renderArray, regBasis, MU, regPara, imageSize);