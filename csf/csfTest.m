%% TODO:
% QUEST+ routine for threshold measurement

%% Setup variables
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Create a mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'integrationTime', 0.05);
retina.visualizeMosaic();

%% Stimulus generation test
background = stimulusCSF('L+M+S', 0, 1);
target = stimulusCSF('L+M+S', 0.05, 1);

background = invGammaCorrection(background, display.CRT12BitDisplay);
target = invGammaCorrection(target, display.CRT12BitDisplay);

figure(); subplot(1, 2, 1);
imshow(background, 'initialMagnification', 200);

subplot(1, 2, 2);
imshow(target, 'initialMagnification', 200);

AFC(retina, background, target)

%% Plot CSF
crst = 1e-4 : 2e-4 : 5e-2;
pCorr1 = zeros(size(crst));
pCorr2 = zeros(size(crst));

background = stimulusCSF(0, 1);
for idx = 1:length(crst)
    fprintf('%d / %d \n', idx, length(crst));
    pCorr1(idx) = runAFC(retina, background, 'L+M+S', crst(idx), 1);
end

background = stimulusCSF(0, 5);
for idx = 1:length(crst)
    fprintf('%d / %d \n', idx, length(crst));
    pCorr2(idx) = runAFC(retina, background, 'L+M+S', crst(idx), 5);
end

figure();
plot(crst, pCorr1, 'LineWidth', 2); hold on;
plot(crst, pCorr2, 'LineWidth', 2);

%% Helper function
function pCorr = runAFC(retina, background, targetType, targetCrst, targetFreq)    
    target = stimulusCSF(targetType, targetCrst, targetFreq);
    pCorr = AFC(retina, background, target);
end


function pCorr = AFC(retina, background, target)

nTrial = 100;
[~, ~, ~, rateB] = retina.compute(background);
[~, ~, ~, rateT] = retina.compute(target);


sampleT = poissrnd(repmat(rateT, 1, nTrial));
llB = sum(log(poisspdf(sampleT, repmat(rateB, 1, nTrial))), 1);
llT = sum(log(poisspdf(sampleT, repmat(rateT, 1, nTrial))), 1);
logRatio = llT - llB;

pCorr = (sum(logRatio > 0) / nTrial);

end
