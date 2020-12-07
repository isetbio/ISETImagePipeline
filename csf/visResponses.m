%% Vis responses for one contrast level
oneSided = false;
[nullTemp, stimTemp, nullNoise, stimNoise] = showRecon(responses, oneSided);
[DV_CR, DV_IR, DV] = decisionVar(stimTemp, nullTemp, stimNoise, nullNoise, oneSided);

%% Show images
oneSided = false;
responses = responseObj{1}{1};
[nullTemp, stimTemp, nullNoise, stimNoise] = showRecon(responses, oneSided);
[DV_CR, DV_IR, DV] = decisionVar(stimTemp, nullTemp, stimNoise, nullNoise, oneSided);

%% Helper function
function [nullTemp, stimTemp, nullNoise, stimNoise] = showRecon(responses, oneSided)
imageSize = [50, 50, 3];

nullTemp = responses.inSampleNullStimResponses('none');
stimTemp = responses.inSampleTestStimResponses('none');

nullNoise = responses.outOfSampleNullStimResponses('random');
stimNoise = responses.outOfSampleTestStimResponses('random');

figure();
subplot(1, 2, 1);
imshow(reshape(nullTemp, imageSize));

subplot(1, 2, 2);
imshow(reshape(stimTemp, imageSize));

if ~ oneSided
    figure();
    for idx = 1:9
        subplot(3, 3, idx);
        imshow(reshape(nullNoise(:, idx), imageSize));
    end
end

figure();
for idx = 1:9
    subplot(3, 3, idx);
    imshow(reshape(stimNoise(:, idx), imageSize));
end

end

function [DV_CR, DV_IR, DV] = decisionVar(stimTemp, nullTemp, stimNoise, nullNoise, oneSided)

DV_CR = zeros(1, size(stimNoise, 2));
DV_IR = zeros(1, size(stimNoise, 2));

for idx = 1:size(stimNoise, 2)
        
    if oneSided 
        DV_CR(idx) = norm(stimNoise(:, idx) - stimTemp);
        DV_IR(idx) = norm(stimNoise(:, idx) - nullTemp);
    else
        DV_CR(idx) = norm([stimNoise(:, idx) - stimTemp; nullNoise(:, idx) - nullTemp]);
        DV_IR(idx) = norm([stimNoise(:, idx) - nullTemp; nullNoise(:, idx) - stimTemp]);
    end            
end

DV = DV_CR - DV_IR;

end
