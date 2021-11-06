function anomalousRun(lightLevel, outName)

imageSize = [128, 128, 3];

% load variables
prior = load('sparsePrior.mat');
regBasis = prior.regBasis;
mu = prior.MU;

input = load('inputImage_128.mat');
inputLinear = input.inputLinear;

renderNormal = load('render_normal.mat');
renderNormal = lightLevel * renderNormal.renderNormal;

renderAnoma = load('render_anomalous.mat');
renderAnoma = lightLevel * renderAnoma.renderAnomalous;

% 0.001 for night light (0.01 render), 0.0075 for dim light (0.2 render),
% 0.01 for regular light (render), 0.05 for bright light (5 render)
regPara = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.05, 0.075, 0.1];
outputArray = cell(1, length(regPara));

for idx = 1:length(regPara)
    outputArray{idx} = reconArray(inputLinear, {renderNormal, renderAnoma}, ...
                                    regBasis, mu, regPara(idx), imageSize);
end

save(outName, 'outputArray', '-v7.3');

end

