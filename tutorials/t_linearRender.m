% Demonstrate how to compute the linear render matrix approximation of the
% ISETBio pipeline.
%
% Description:
%    This tutorial shows how we calculate a linear render matrix
%    approximation of the ISETBio, with random (white noise) images.

%% Generate desired cone mosaic
% ConeResponse is a wrapper class of isetbio calculation for RGB image
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.2);

% We mainly need the last two variable, linearImage and coneVec for our
% calculation. linearImage is the input image after gamma function, since
% the linear calculation starts after that; coneVec is the pattern of cone
% excitation but expressed simply as a vector (instead of a two-dimensional
% matrix).

% Use smaller image so this example won't take forever to run, but this
% calculation can be done for very large images (e.g., 200*200*3).
imageSize = [16, 16, 3];
testImage = rand(imageSize);
[excitation, theOI, linearImage, coneVec] = retina.compute(testImage);

%% Generate a "training set" for regression
% Since the dimensionality of input image is 16*16*3, we want the
% size of the trianing set to be at least larger than that.
nTrain = 1e3;
allConeVec   = zeros(nTrain, length(coneVec));
allLinearImg = zeros(nTrain, length(linearImage(:)));

parfor idx = 1:nTrain
    % In practice, we found that white noise images produce better results
    % and are easier to work with (since we can simply generate them on the
    % fly).
    input = rand(imageSize);
    [~, ~, linearImage, coneVec] = retina.compute(input);
    
    allConeVec(idx, :)   = coneVec;
    allLinearImg(idx, :) = linearImage(:);
end

%% Regression
% Calculate the linear regression. RegressionEstimator is a wrapper class
% for a standard SVD based least square estimation.
regEstimator = RegressionEstimator(allLinearImg, allConeVec);
renderMatrix = regEstimator.W';

%% Test the render matrix
% Here we use a newly generated random (white noise) image as the test, but
% could (and should) also use natural images.
testImage = rand(imageSize);

% Generate coneVec with ISETBio
[~, ~, linearImage, coneVec] = retina.compute(testImage);

% Generate coneVec with our render matrix
coneVecRender = renderMatrix * linearImage(:);

% Check it is indeed perfect
figure(); axis equal;
scatter(coneVec, coneVecRender);
fprintf('Average Difference (in percentage) : %.4f%% \n', ...
    mean(abs(coneVec - coneVecRender) ./ coneVec) * 100);
