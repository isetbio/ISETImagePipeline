% Compute comparison methods for dichromatic visualization
%
% 

% History:
%   12/10/21  dhb  Read it and added comments.

%% Clear and close
clear; close all;

%% Compute the matrix from linear RGB -> LMS
display = displayCreate('CRT12BitDisplay');

% Allocate
matrix_LinRGBToLMS = zeros(3, 3);

% Loop over RGB and create a one pixel image that
% with R, G, or B set to max value, and pull out 
% the corresponding LMS coordinates. 
for channelIdx = 1:3
    oneChannelImg = zeros([1, 1, 3]);
    oneChannelImg(:, :, channelIdx) = 1;

    % Create a scene from the displayed iamge.  Although we think
    % this method gamma corrects, gamma correction maps 1 to 1.
    [scene, ~, ~] = sceneFromFile(oneChannelImg, 'rgb', [], display);

    % Get LMS values using sceneGet and fill in the column of the matrix
    val = sceneGet(scene, 'lms');
    matrix_LinRGBToLMS(:, channelIdx) = val;
end
wave = sceneGet(scene, 'wave');

% Let's make sure we understand what the above is doing, by
% computing the matrix directly from the display channel spd
% and the PTB Stockman-Sharpe 2-deg fundamentals.  Agreement
% is OK enough to make us think it's just a numerical choice
% somewhere, probably in the tabulated values of the fundamentals
% between PTB and ISET which are different.
CHECK = false;
if (CHECK)
    T_cones = WlsToS(wave);
    condData = load('T_cones_ss2.mat');
    T_cones = SplineCmf(condData.S_cones_ss2,condData.T_cones_ss2,T_cones);
    matrix_Check = T_cones*display.spd*T_cones(2);
    for ii = 1:3
        for jj = 1:3
            matrix_Check2(ii,jj) = trapz(display.wave,T_cones(ii,:)' .* display.spd(:,jj));
        end
    end
    (matrix_LinRGBToLMS-matrix_Check)./matrix_Check;
    (matrix_LinRGBToLMS-matrix_Check2)./matrix_Check;
    T_cones_iset = ieReadSpectra('stockman', wave)';
    figure; hold on;
    plot(wave,T_cones_iset,'r');
    plot(wave,T_cones,'k');
end

%% White point of the display in LMS
whiteImg = ones([1, 1, 3]);
[scene, ~, ~] = sceneFromFile(whiteImg, 'rgb', [], display);
whiteLMS = sceneGet(scene, 'lms');

%% Visulization, Brettel (aka 'brettel') method
imSize = [128, 128, 3];
load('inputImage_128.mat');

figure();
plotAxis = tight_subplot(3, size(inputLinear, 1), ...
    [.01 .01], [.01 .01], [.01 .01]);

method = 'brettel';
for dichromType = 1:3
    for imageIdx = 1:size(inputLinear, 1)
        % linear RGB
        image = reshape(inputLinear(imageIdx, ...
            :, :, :), imSize);

        % linear RGB -> LMS
        lms = linearTrans(image, imSize, matrix_LinRGBToLMS);

        % LMS -> dichromatic LMS, Brettel method
        lmsDichma = lms2lmsDichromat(lms, dichromType, method, whiteLMS);

        % LMS -> Linear RGB
        rgb = linearTrans(lmsDichma, imSize, inv(matrix_LinRGBToLMS));

        % plot image
        plotIdx = (dichromType - 1) * size(inputLinear, 1) + imageIdx;
        axes(plotAxis(plotIdx));
        imshow(gammaCorrection(rgb, display));
    end
end

%% Visulization, Jiang et al. (aka 'linear') method
imSize = [128, 128, 3];
load('inputImage_128.mat');

figure();
plotAxis = tight_subplot(3, size(inputLinear, 1), ...
    [.01 .01], [.01 .01], [.01 .01]);

method = 'linear';
for dichromType = 1:3
    for imageIdx = 1:size(inputLinear, 1)
        % linear RGB
        image = reshape(inputLinear(imageIdx, ...
            :, :, :), imSize);
    
        % linear RGB -> LMS
        lms = linearTrans(image, imSize, matrix_LinRGBToLMS);
        
        % LMS -> dichromatic LMS, linear method, Jiang et al., 2015
        lmsDichma = lms2lmsDichromat(lms, dichromType, method);

        % LMS -> Linear RGB
        rgb = linearTrans(lmsDichma, imSize, inv(matrix_LinRGBToLMS));

        % plot image
        plotIdx = (dichromType - 1) * size(inputLinear, 1) + imageIdx;
        axes(plotAxis(plotIdx));
        imshow(gammaCorrection(rgb, display));
    end
end

%% Helper function
function output = linearTrans(input, imSize, matrix)

input = reshape(input, [imSize(1) * imSize(2), 3]);
output = matrix * input';
output = reshape(output', imSize);

end
