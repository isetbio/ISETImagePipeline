%% Compute the matrix from RGB -> LMS
display = displayCreate('CRT12BitDisplay');

% Matrix that takes RGB to LMS
matrix = zeros(3, 3);

% R, G, B
for idx = 1:3
    image = zeros([1, 1, 3]);
    image(:, :, idx) = 1;

    [scene, ~, ~] = sceneFromFile(image, 'rgb', [], display);
    val = sceneGet(scene, 'lms');

    % fill in the column
    matrix(:, idx) = val;
end

% White point of the display
image = ones([1, 1, 3]);
[scene, ~, ~] = sceneFromFile(image, 'rgb', [], display);
whitePt = sceneGet(scene, 'lms');

%% Visulization, Brettel method
imSize = [128, 128, 3];
load('inputImage_128.mat');

figure();
plotAxis = tight_subplot(3, size(inputLinear, 1), ...
    [.01 .01], [.01 .01], [.01 .01]);

method = 'brettel';
for type = 1:3
    for idx = 1:size(inputLinear, 1)
        % linear RGB
        image = reshape(inputLinear(idx, ...
            :, :, :), imSize);

        % linear RGB -> LMS
        lms = linearTrans(image, imSize, matrix);

        % LMS -> dichromatic LMS, Brettel method
        lmsDichma = lms2lmsDichromat(lms, type, method, whitePt);

        % LMS -> Linear RGB
        rgb = linearTrans(lmsDichma, imSize, inv(matrix));

        % plot image
        plotIdx = (type - 1) * size(inputLinear, 1) + idx;
        axes(plotAxis(plotIdx));
        imshow(gammaCorrection(rgb, display));
    end
end

%% Visulization, linear method
imSize = [128, 128, 3];
load('inputImage_128.mat');

figure();
plotAxis = tight_subplot(3, size(inputLinear, 1), ...
    [.01 .01], [.01 .01], [.01 .01]);

method = 'linear';
for type = 1:3
    for idx = 1:size(inputLinear, 1)
        % linear RGB
        image = reshape(inputLinear(idx, ...
            :, :, :), imSize);
    
        % linear RGB -> LMS
        lms = linearTrans(image, imSize, matrix);
        
        % LMS -> dichromatic LMS, linear method, Jiang et al., 2015
        lmsDichma = lms2lmsDichromat(lms, type, method);

        % LMS -> Linear RGB
        rgb = linearTrans(lmsDichma, imSize, inv(matrix));

        % plot image
        plotIdx = (type - 1) * size(inputLinear, 1) + idx;
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
