function status = ipipeRunExamplesAll
%% Run all the examples in the ISETImagePipeline tree
%
% Syntax:
%     ipipeRunExamplesAll
%
% Description:
%     Run all the examples in the ISETImagePipeline  tree,
%     excepthose that contain a line of the form
%     "% ETTBSkip"
%
% Inputs:
%    None.
%
% Outputs:
%    status    - 1 if all examples run OK, 0 otherwise.
%
% Optional key/value pairs:
%    None.
%

% History:
%   01/17/18  dhb  Wrote it.

[~, functionStatus] = ExecuteExamplesInDirectory(tbLocateProject('ISETImagePipeline'),'verbose',false);
status = all(functionStatus ~= -1);