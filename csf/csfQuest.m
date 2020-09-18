%% Setup variables
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Create a mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
                      'fovealDegree', 0.5, 'display', display.CRT12BitDisplay, 'integrationTime', 0.05);
retina.visualizeMosaic();

%% CSF with QUEST procedure
observer = PoissonTemplateObserver(retina, display.CRT12BitDisplay, 'L+M+S', 1);
questData = qpInitialize('stimParamsDomainList', {-5 : 0.01 : -1}, ...
                         'psiParamsDomainList',  {-5 : 0.01 : -1, 2:5, 0.5, 0:0.01:0.05});

nTrial = 64;
for idx = 1 : nTrial
    % Get stimulus for this trial
    targetCrst = qpQuery(questData);
    
    % Simulate outcome
    outcome = observer.singleTrial(power(10, targetCrst)) + 1;
    
    % Update quest data structure
    questData = qpUpdate(questData, targetCrst, outcome); 
end

%% Find aximum likelihood fit.  Use psiParams from QUEST+ as the starting
psiParamsFit = qpFit(questData.trialData, questData.qpPF, psiParamsQuest, questData.nOutcomes,...
    'lowerBounds', [1e-4 1 0.5 0],'upperBounds',[1e-1 5 0.5 0.05]);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    psiParamsFit(1), psiParamsFit(2), psiParamsFit(3), psiParamsFit(4));

%% Plot of trial locations with maximum likelihood fit
figure; clf; hold on
stimCounts = qpCounts(qpData(questData.trialData),questData.nOutcomes);
stim = [stimCounts.stim];
stimFine = linspace(1e-4, 1e-1, 1000)';
plotProportionsFit = qpPFWeibull(stimFine,psiParamsFit);
for cc = 1:length(stimCounts)
    nTrials(cc) = sum(stimCounts(cc).outcomeCounts);
    pCorrect(cc) = stimCounts(cc).outcomeCounts(2)/nTrials(cc);
end
for cc = 1:length(stimCounts)
    h = scatter(stim(cc),pCorrect(cc),100,'o','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1],...
        'MarkerFaceAlpha',nTrials(cc)/max(nTrials),'MarkerEdgeAlpha',nTrials(cc)/max(nTrials));
end
plot(stimFine,plotProportionsFit(:,2),'-','Color',[1.0 0.2 0.0],'LineWidth',3);
xlabel('Stimulus Value');
ylabel('Proportion Correct');
xlim([1e-4 1e-1]); ylim([0 1]);
title({'Estimate Weibull threshold, slope, and lapse', ''});
drawnow;
