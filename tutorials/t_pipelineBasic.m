% Demonstrate basic computation of retinal image and mosaic responses from RGB
%
% Description:
%    This tutorial shows how we process a set of RGB images through ISETBIO
%    to give us a dataset that has the RGB image, the retinal image, and
%    the mosaic responses.  We can then use this dataset in various machine
%    learning applications, or other fun stuff.
%
%    The work is done by the ConeReponse object, which itself lives in the
%    ISETPipelineToolbox.
%
%    This tutorial uses the CIFAR image data set, which we have already
%    stored in .mat file format in the project data directory.  This dataset
%    is described here:
%      https://www.cs.toronto.edu/~kriz/cifar.html
%    We have combined the CIFAR-10 and CIFAR-100 datasets and stored them 
%    as a 1000000 by 3072 matrix, where each row represents one 32 by 32 by
%    3 image.  Unpack a row to an image using
%      image = reshape(row, [32, 32, 3])

% History:
%   01/30/19  lz       Wrote it.
%   01/31/19  lz, dhb  Comments, clean up etc.
%   02/04/19  lz       Update with cone excitation vector added etc.

%% Clear etc
clear; close all;
rng(1);

%% Setup input and output directory
%
% We use the dataBaseDir root set by the project preference, and then work
% down a directory tree from there.
projectName = 'ISETImagePipeline';
thisImageSet = 'CIFAR_all';
theDir = 'true_1_tutorial';
dataBaseDir = getpref(projectName, 'dataDir');
dataFileIn = fullfile(dataBaseDir, thisImageSet, 'image_cifar_all.mat');
dataDirOut = fullfile(dataBaseDir, thisImageSet, theDir);
if (~exist(dataDirOut, 'dir'))
    mkdir(dataDirOut);
end

%% Load in the dataset data and create ConeResponse object
if (~exist(dataFileIn, 'file'))
    error('Need to get CIFAR dataset in .mat format into expected location before this will work');
end
load(dataFileIn);

%% Create the compute object, show cone mosaic being used
fprintf('Constructing ConeResponse object ... \n');  
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true);
fprintf('Finish Constructing ConeResponse object. \n');
retina.visualizeMosaic();

%% Define number of images
nImages = 5;

%% Compute optical image and mosaic excitation pattern for each input image
for imgIdx = 1:nImages
    % The CIFAR dataset comes in with each image as a row vector.  Here we
    % pull out just one image. We rely on the fact that we know the size of
    % the image in each row. THis knowledge is specific to the CIFAR
    % dataset as we downloaded it.
    inputImage = reshape(image_all(imgIdx, :), [32, 32, 3]);
    
    % The ConeResponse object does all the work, and gives as the result
    [excitationImage, oi, allCone, L, M, S] = retina.compute(inputImage);
    
    % Write to output directory
    %
    % First time through, save full optical image structure.  This is
    % constant across images except for the photon data.  To save space, we
    % only store one of these, and then store the photon data separately
    % for each input image
    if (imgIdx == 1)
        % Save optical image structure
        oistructFile = 'oiStruct.mat';
        save(fullfile(dataDirOut, oistructFile), 'oi');
        
        % Save one copy of the mosaic object that was used
        % to compute everything
        mosaic = retina.getMosaic();
        mosaicFile = 'mosaicStruct.mat';
        save(fullfile(dataDirOut, mosaicFile), 'mosaic');
        
        % Show input image, optical image, and cone mosaic excitation        
        retina.visualizeOI();
        retina.visualizeExcitation();
    end
    
    % For optical image, just store the photon map (isomerizations)
    oiFile = sprintf('oiPhotons_%d.mat', imgIdx);
    photons = oi.data.photons;
    save(fullfile(dataDirOut, oiFile), 'photons');
    
    % Save excitation image
    excitationFile = sprintf('excitations_%d.mat', imgIdx);
    save(fullfile(dataDirOut, excitationFile), 'excitationImage');
      
    % Save excitations as a vector with length number of cones. This might 
    % be useful, for example, if we are using them as likelihood function, 
    % or in regression settings    
    coneVectorFile = sprintf('coneVector_%d.mat', imgIdx);
    save(fullfile(dataDirOut, coneVectorFile), 'allCone');    
    
end


