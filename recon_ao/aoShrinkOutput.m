
%% Get list of output files
% 
% Make pwd two levels above their location to find
fileList = dir('*/*/xRunOutput.mat');

for ii = 1:length(fileList)
    fprintf('Fixing %d of %d, %s\n',ii,length(fileList),fileList(ii).folder)

    fixOutput(fullfile(fileList(ii).folder,fileList(ii).name));
end

function fixOutput(name)


fprintf('\tLoading ... ')
load(name)
fprintf('loaded ... ')
clear estimator
clear reconImageLinearTemp psfSupportTemp initImageLinearTemp
clear tempFig theAxes theFIg axesHandle temp initSceneTemp
clear forwardConeMosaic reconConeMosaic
save(name,'-v7.3');
fprintf('done\n');
end