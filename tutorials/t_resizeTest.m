

prBase.nPixels = 128;
prBase.imageName = 'puppy';
prBase.imageType = 'jpeg';

% JPEG options include 'butterfly', 'dragonfly', 'panda', 'pizza'
% 'puppy', and 'sunflower'.  Not all of these images are square, so
% we select out the largest square in the center before resizing.
theImageRGB = imread(fullfile(prBase.aoReconDir,'images',[prBase.imageName '.jpeg']),'jpeg');
[m,n,k] = size(theImageRGB);
minDim = min([m,n]); 
mSpace = minDim/2; nSpace = minDim/2;
lowM = round(m/2-mSpace)+1; highM = lowM+minDim-1; lowN = round(n/2-nSpace)+1; highN = lowN+minDim-1;
fullImage = theImageRGB(lowM:highM,lowN:highN,:);
smallImage = imresize(fullImage,'OutputSize',[prBase.nPixels prBase.nPixels]);
figure; imshow(fullImage); title('Full Square');
f = figure; imshow(smallImage); title('Resized');
imwrite(smallImage,'/Users/dhb/Desktop/test.tiff','tiff');
saveas(f,'/Users/dhb/Desktop/testfig.tiff','tiff');