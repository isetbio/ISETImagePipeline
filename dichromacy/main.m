% Saved variable
load('sparsePrior.mat');
load('inputImage_128.mat');

load('render_normal.mat');
load('render_anomalous.mat');

% render mtx for dim light
renderNormal = 0.01 * double(renderNormal);
renderAnomalous = 0.01 * double(renderAnomalous);

% 0.001 for night light (0.01 render), 0.0075 for dim light (0.2 render),
% 0.01 for regular light (render), 0.05 for bright light (5 render)
regPara = 0.001;
outputArray = reconArray(inputLinear, {renderNormal, renderAnomalous}, regBasis, MU, regPara, imageSize);

save('output_noise_nightlight.mat', 'outputArray', '-v7.3');

%% archive
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