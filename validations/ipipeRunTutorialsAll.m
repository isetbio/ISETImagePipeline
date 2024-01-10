function status = ipipeRunTutorialsAll
% Run all ISETImagePipeline tutorials and collect up status
%
% Syntax
%    status = ipipeRunTutorialsAll
%
% Description
%    Run all of the ISETImagePipeline  tutorials that we think should work, and print
%    out a report at the end as to whether they threw errors, or not.
%    Scripts inside of ISETImagePipeline/tutorials are run, except that
%    scripts within the directory 'underDevelopment' and 'support' are skipped.
%
% Inputs:
%    None
%
% Outputs:
%    status    - Set to true if all tutorials ran OK, false otherwise.
% 
% See also:
%

% History:
%   07/26/17  dhb  Wrote this, because we care.

% User/project specific preferences
p = struct(...
    'rootDirectory',            fileparts(which(mfilename())), ...
    'tutorialsSourceDir',       fullfile(tbLocateProject('ISETImagePipeline'), 'tutorials') ...                % local directory where tutorial scripts are located
    );

%% List of scripts to be skipped
%
% Anything with this in its path name is skipped.
scriptsToSkip = {...
    'underDevelopment' ...
    'support' ...
    };

%% Use UnitTestToolbox method to do this.
status = UnitTest.runProjectTutorials(p, scriptsToSkip, 'All');

end