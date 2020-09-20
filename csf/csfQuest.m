%% Setup variables
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

%% Create a mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 0.5, 'display', display.CRT12BitDisplay, 'integrationTime', 0.05);
retina.visualizeMosaic();

%% CSF with QUEST single trial procedure
observer = PoissonTemplateObserver(retina, display.CRT12BitDisplay, 'L+M+S', 1);

crstDomain = -5 : 0.025 : 0;
slopeRange = 0.1 : 0.5 : 100;

questData = qpInitialize('stimParamsDomainList', {crstDomain}, ...
    'psiParamsDomainList',  {crstDomain,  slopeRange, 0.5, 0.0});

nTrial = 100;
for idx = 1 : nTrial
    fprintf('Trial: %d / %d \n', idx, nTrial);
    % Get stimulus for this trial
    targetCrst = qpQuery(questData);
    
    % Simulate outcome
    outcome = observer.singleTrial(10 ^ targetCrst) + 1;
    
    % Update quest data structure
    questData = qpUpdate(questData, targetCrst, outcome);
end

%% CSF with QUEST multi trial procedure
questData = qpInitialize('stimParamsDomainList', {crstDomain}, ...
    'psiParamsDomainList',  {crstDomain,  slopeRange, 0.5, 0.0});

nTrial = 32;
nRepeat = 20;
for idx = 1 : nTrial
    fprintf('Trial: %d / %d \n', idx, nTrial);
    % Get stimulus for this trial
    targetCrst = qpQuery(questData);
    
    % Simulate outcome
    [~, outcome] = observer.multiTrial(10 ^ targetCrst, nRepeat);
    outcome = outcome + 1;
    
    % Update quest data structure
    for trial = 1:length(outcome)
        questData = qpUpdate(questData, targetCrst, outcome(trial));
    end
    
end

%% Data
data = questData.trialData;
nTrial = length(data);

stim = zeros(1, nTrial);
response = zeros(1, nTrial);

for idx = 1:nTrial
    stim(idx) = data(idx).stim;
    response(idx) = data(idx).outcome - 1;
end

figure();
stimVal = unique(stim);
for idx = 1:length(stimVal)
    prop = response(stim == stimVal(idx));       
        
    scatter(stimVal(idx), sum(prop) / length(prop), 25 * 100 / nTrial * length(prop), ...
        'MarkerEdgeColor', zeros(1, 3), 'MarkerFaceColor', ones(1, 3) * 0.5, 'MarkerFaceAlpha', 0.5);
    hold on;
end

psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);

fprintf('Max posterior QUEST+ parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    psiParamsQuest(1), psiParamsQuest(2), psiParamsQuest(3), psiParamsQuest(4));

psiParamsFit = qpFit(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', [min(crstDomain) min(slopeRange) 0.5 0.0],'upperBounds',[max(crstDomain) max(slopeRange) 0.5 0.0]);

fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    psiParamsFit(1), psiParamsFit(2), psiParamsFit(3), psiParamsFit(4));

hold on;
stimSpace = crstDomain;

fitCurve = qpPFWeibull(stimSpace', psiParamsFit);
plot(stimSpace, fitCurve(:,2), '-', 'Color', [1.0 0.2 0.0], 'LineWidth',3);
