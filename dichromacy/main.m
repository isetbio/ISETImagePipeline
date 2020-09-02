% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

% % M1
% load('retina_render_M1.mat');
% 
% regPara = 5e-4;
% outputArray_M1 = reconArray(inputLinear, renderArray, regBasis, MU, regPara, imageSize);
% 
% save('output_m1.mat', 'outputArray_M1', '-v7.3');
% 
% % M2
% load('retina_render_M2.mat');
% 
% regPara = 5e-4;
% outputArray_M2 = reconArray(inputLinear, renderArray, regBasis, MU, regPara, imageSize);
% 
% save('output_m2.mat', 'outputArray_M2', '-v7.3');

load('render_anomalous.mat');
regPara = 0.05;
outputArray = reconArray(inputLinear, {renderAnomalous}, regBasis, MU, regPara, imageSize);

save('output_noise.mat', 'outputArray', '-v7.3');