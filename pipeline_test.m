% Setup input and output directory
toolboxName = 'ISETImagePipeline';
dataBaseDir = getpref(toolboxName,'dataDir');
dataDirIn   = strcat(dataBaseDir, '/CIFAR_all/image_cifar_all.mat');
dataDirOut  = strcat(dataBaseDir, '/CIFAR_all/');

% Load data and ConeResponse object
load(dataDirIn);
retina = ConeResponse();

nImage     = 100;
excitaions = zeros(nImage, 740, 740);

% Compute mosaic excitation pattern and optical image
for imgIdx = 1 : nImage
    input = reshape(image_all(imgIdx, :, :), [32, 32, 3]);
    [exci, oi] = retina.compute(input);
    
    excitaions(imgIdx, :, :) = exci;    
    if imgIdx == 1
        optical_img = oi;
        optical_img = repmat(optical_img, 1, nImage);
    else
        optical_img(imgIdx) = oi;
    end        
end

% Write to output directory
save(strcat(dataDirOut, 'cifar_exci_false_1.mat'), 'excitaions');
save(strcat(dataDirOut, 'cifar_oi_false_1.mat'),   'optical_img');
