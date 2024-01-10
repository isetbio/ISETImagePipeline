%% Generate cone mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, 'fovealDegree', 1.0, 'pupilSize', 2.0);

imageSize = [128, 128, 3];
input = ones(imageSize) * 1;
[~, ~, ~, coneVec] = retina.compute(input);

%% Linearity test
% Run N images with only one pixel set to a value of 1
nInput = 10;
coneResp = zeros(length(coneVec), nInput);

for idx = 1 : nInput
    input = zeros(imageSize);
    input(idx) = 1;
    
    [~, ~, ~, coneVec] = retina.compute(input);
    
    % Save cone response for each image
    coneResp(:, idx) = coneVec;
end

% Run 1 image with N pixel set to the value of 1
input = zeros(imageSize);
input(1:nInput) = 1;
[~, ~, ~, coneVec] = retina.compute(input);

%% The sum of the cone response should be equal to the cone response to the sum of images
% Scatter plot
figure(); hold on;
scatter(coneVec, sum(coneResp, 2));

refPoint = [0, 20];
plot(refPoint, refPoint);
axis square;
xlim(refPoint);
ylim(refPoint);
