%% Load dataset
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

figure(); 
subplot(1,3,1); plot(display.CRT12BitDisplay.gamma(:,1), 'r-'); 
subplot(1,3,2); plot(display.CRT12BitDisplay.gamma(:,2), 'g-'); 
subplot(1,3,3); plot(display.CRT12BitDisplay.gamma(:,3), 'b-');

retina = ConeResponse('eccBasedConeDensity', false, 'eccBasedConeQuantal', false, ...
                      'display', display.CRT12BitDisplay);

%% Compute test image
imageSize = [32, 32, 3];
testImage1 = rand(imageSize) * 0.6;
[~, ~, testLinearImage1, testConeVec1] = retina.compute(testImage1);

%% Visulization
figure;
scatter(testImage1(:), testLinearImage1(:));

figure;
invTestImage1 = invGammaCorrection(testLinearImage1, display.CRT12BitDisplay);
scatter(testImage1(:), invTestImage1(:));

%% Linear test
testImage2 = rand(imageSize) * 0.6;
[~, ~, testLinearImage2, testConeVec2] = retina.compute(testImage2);

%% Sum Image
sumLinear = testLinearImage1 + testLinearImage2;
sumImage  = invGammaCorrection(sumLinear, display.CRT12BitDisplay);
[~, ~, testLinearImage3, testConeVec3] = retina.compute(sumImage);

%% Input check
figure;
scatter(sumLinear(:), testLinearImage3(:));

%% Visulization
figure;
sumConeVec = testConeVec1 + testConeVec2;
scatter(sumConeVec(:), testConeVec3(:));
xlim([0, 7000]); ylim([0, 7000]);
axis square;