%% Image Preprocessing and Visualization

load('./cifar-all-mat/data_all.mat');
nimg = size(data_all); nimg = nimg(1);
image_all = zeros([nimg, 16*16]);

for idx = 1:100000
    img = reshape(data_all(idx, :), [32, 32, 3]);

    for dim = 1:3
        img(:, :, dim) = img(:, :, dim)';
    end

img = rgb2gray(imresize(im2double(img), 0.5));
image_all(idx, :) = reshape(img, [1, 16*16]);
end

%% Check 
load('./image_all.mat');
for loop = 1:1000
    idx = randi(100000, 1);
    img = reshape(image_all(idx, :), [32, 32, 3]);
    imshow(img, 'InitialMagnification', 500);
    pause(1);
end

%% Check 
load('./cifar-all-mat/image_all_16.mat');
for loop = 1:1000
    idx = randi(100000, 1);
    img = reshape(image_all(idx, :), [16, 16, 3]);
    
    imshow(img, 'InitialMagnification', 300);
    pause(0.1);
end

%% Check 
load('./cifar-all-mat/image_all_16_gray.mat');
for loop = 1:1000
    idx = randi(100000, 1);
    img = reshape(image_all(idx, :), [16, 16]);
    
    imshow(img, 'InitialMagnification', 300);
    pause(0.1);
end