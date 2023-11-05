function anomalousRunGray(lightLevel, outName)
% light level:
% night light (0.01 render), dim light (0.2 render),
% regular light (render), bright light (5 render)

imageSize = [128, 128, 3];

% load variables
prior = load('sparsePrior.mat');
regBasis = prior.regBasis;
mu = prior.MU;

input = load('inputImage_128.mat');
inputLinear = input.inputLinear;

renderNormal = load('render_normal.mat');
renderNormal = lightLevel * renderNormal.renderNormal;

regPara = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.05, 0.075, 0.1];
outputArray = cell(1, length(regPara));

for idx = 1:length(regPara)
    fprintf('\n \n Iteration %d \n \n', idx);
    outputArray{idx} = reconArrayGray(inputLinear, {renderNormal}, ...
                                    regBasis, mu, regPara(idx), imageSize);
end

save(outName, 'outputArray', '-v7.3');

end

