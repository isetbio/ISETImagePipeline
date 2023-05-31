function prBase = prBaseDefaults;
% Set defaults for prBase parameters

%% Initialize
prBase = struct;

%% Point at directory with data files for this subproject
%
% This will allow us to load in project specific precomputed information.
% Also records initials of version editors, otherwise set to 'main'
prBase.aoReconDir = getpref('ISETImagePipeline','aoReconDir');