function dataOut = reconNeuralEngine(reconObj, coneRespObj, dispOption, ~, ~, sceneSequence, ~, instancesNum, varargin)

if (nargin == 0)
    dataOut = struct();
    return;
end

p = inputParser;
p.addParameter('noiseFlags', {'random'});
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

noiseFlags = p.Results.noiseFlags;
imageSize = reconObj.Size;
mask = makeMast(imageSize);

[~, coneVec] = coneRespObj.computeWithScene(sceneSequence{:});

nIter = 600;
reconResponses = containers.Map();
for idx = 1:length(noiseFlags)
    if strcmp(noiseFlags{idx}, 'none')
        
        % reset random number generator for template
        rng('default');
        init = rand([prod(imageSize), 1]);
        recon = (reconObj.estimate(coneVec, nIter, init, true, 1.0, dispOption)) .* mask;
        reconResponses(noiseFlags{idx}) = recon(:);
        
    elseif strcmp(noiseFlags{idx}, 'random')
        
        recon = zeros([prod(imageSize), instancesNum]);
        parfor itr = 1:instancesNum
            init = rand([prod(imageSize), 1]);
            output = (reconObj.estimate(poissrnd(coneVec), nIter, init, true, 1.0, dispOption)) .* mask;
            recon(:, itr) = output(:);
        end
        reconResponses(noiseFlags{idx}) = recon;
        
    end
end
% coneVec = repmat(coneVec, 1, instancesNum);
dataOut = struct(...
    'neuralResponses', reconResponses, ...
    'temporalSupport', []);

end

function maskImage = makeMast(imageSize)

maskImage = ones(imageSize);

center = imageSize(1) * 0.5;
radius = imageSize(1) * 0.425;

for i = 1:imageSize(1)
    for j = 1:imageSize(1)
        dist = sqrt((i - center) ^ 2 + (j - center) ^ 2);
        if (dist > radius)
            maskImage(i, j, :) = 0;
        end
    end
end

end
