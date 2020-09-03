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

load('render_normal.mat');
% load('render_anomalous.mat');

% 0.01 (dim light), 0.05, 
% renderNormal = 0.2 * double(renderNormal);
% renderAnomalous = 0.2 * double(renderAnomalous);

regPara = 0.05;
renderNormal = double(renderNormal);
outputArray = reconArray(inputLinear, {renderNormal}, regBasis, MU, regPara, imageSize);

save('output_noise_low_normal.mat', 'outputArray', '-v7.3');