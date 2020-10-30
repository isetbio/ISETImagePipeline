function dataOut = reconNeuralEngine(coneRespObj, ~, ~, sceneSequence, ~, instancesNum, varargin)

if (nargin == 0)
    dataOut = struct();
    return;
end

p = inputParser;
p.addParameter('noiseFlags', {'random'});
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

noiseFlags = p.Results.noiseFlags;

[~, coneVec] = coneRespObj.computeWithScene(sceneSequence{:});
coneVec = repmat(coneVec, 1, instancesNum);

coneResponses = containers.Map();
for idx = 1:length(noiseFlags)
    if strcmp(noiseFlags{idx}, 'none')
        coneResponses(noiseFlags{idx}) = coneVec;
    elseif strcmp(noiseFlags{idx}, 'random')
        coneResponses(noiseFlags{idx}) = poissrnd(coneVec);
    end
end

dataOut = struct(...
    'neuralResponses', coneResponses, ...
    'temporalSupport', []);

end