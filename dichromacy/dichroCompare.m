%% Compute the matrix from RGB -> LMS
display = displayCreate('CRT12BitDisplay');

matrix = zeros(3, 3);
for idx = 1:3
    image = zeros([1, 1, 3]);
    image(:, :, idx) = 1;

    [scene, ~, ~] = sceneFromFile(image, 'rgb', [], display);
    val = sceneGet(scene, 'lms');

    matrix(:, idx) = val;
end

image = ones([1, 1, 3]);
[scene, ~, ~] = sceneFromFile(image, 'rgb', [], display);
whitePt = sceneGet(scene, 'lms');

%% Visulization, linear method
imSize = [128, 128, 3];
load('inputImage_128.mat');

figure();
plotAxis = tight_subplot(3, size(inputLinear, 1), ...
    [.01 .01], [.01 .01], [.01 .01]);

method = 'linear';
for type = 1:3
    for idx = 1:size(inputLinear, 1)
        image = reshape(inputLinear(idx, ...
            :, :, :), imSize);

        % linear RGB -> LMS
        image = reshape(image, [imSize(1) * imSize(2), 3]);
        lms = matrix * image';
        lms = reshape(lms', imSize);

        % LMS -> dichromatic LMS
        lmsDichma = lms2lmsDichromat(lms, type, method);

        % LMS -> Linear RGB
        lmsDichma = reshape(lmsDichma, [imSize(1) * imSize(2), 3]);
        rgb = matrix \ lmsDichma';
        rgb = reshape(rgb', imSize);

        plotIdx = (type - 1) * size(inputLinear, 1) + idx;
        axes(plotAxis(plotIdx));
        imshow(gammaCorrection(rgb, display));
    end
end

%% Visulization, Brettel method
imSize = [128, 128, 3];
load('inputImage_128.mat');

figure();
plotAxis = tight_subplot(3, size(inputLinear, 1), ...
    [.01 .01], [.01 .01], [.01 .01]);

method = 'brettel';
for type = 1:3
    for idx = 1:size(inputLinear, 1)
        image = reshape(inputLinear(idx, ...
            :, :, :), imSize);

        % linear RGB -> LMS
        image = reshape(image, [imSize(1) * imSize(2), 3]);
        lms = matrix * image';
        lms = reshape(lms', imSize);

        % LMS -> dichromatic LMS
        lmsDichma = lms2lmsDichromat(lms, type, method, whitePt);

        % LMS -> Linear RGB
        lmsDichma = reshape(lmsDichma, [imSize(1) * imSize(2), 3]);
        rgb = matrix \ lmsDichma';
        rgb = reshape(rgb', imSize);

        plotIdx = (type - 1) * size(inputLinear, 1) + idx;
        axes(plotAxis(plotIdx));
        imshow(gammaCorrection(rgb, display));
    end
end
