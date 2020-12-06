%% Vis responses for one contrast level
[nullTemp, stimTemp, nullNoise, stimNoise] = showRecon(responses);
[DV_CR, DV_IR, DV] = decisionVar(stimTemp, nullTemp, stimNoise, nullNoise);

%% Show images
responses = responseObj{1}{1};
[nullTemp, stimTemp, nullNoise, stimNoise] = showRecon(responses);
[DV_CR, DV_IR, DV] = decisionVar(stimTemp, nullTemp, stimNoise, nullNoise);

%% Helper function
function [nullTemp, stimTemp, nullNoise, stimNoise] = showRecon(responses)
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

figure();
for idx = 1:9
    subplot(3, 3, idx);
    imshow(reshape(nullNoise(:, idx), imageSize));
end

figure();
for idx = 1:9
    subplot(3, 3, idx);
    imshow(reshape(stimNoise(:, idx), imageSize));
end
end

function [DV_CR, DV_IR, DV] = decisionVar(stimTemp, nullTemp, stimNoise, nullNoise)

DV_CR = zeros(1, size(nullNoise, 2));
DV_IR = zeros(1, size(nullNoise, 2));
for idx = 1:size(nullNoise, 2)
    DV_CR(idx) = norm([stimNoise(:, idx) - stimTemp; nullNoise(:, idx) - nullTemp]);
    DV_IR(idx) = norm([stimNoise(:, idx) - nullTemp; nullNoise(:, idx) - stimTemp]);
    DV = DV_CR - DV_IR;
end

end
