function [replaceCones, startCones, newCones] = assignCones(chrom)
% Syntax:
%    [replaceCones, startCones, newCones] = assignCones(chrom);
% Synopsis
%    Establish the chromatic assignments for cone mosaics
%
% Description:
%    Provide string for desired chromaticity and return a boolean to
%    indicate desired changes as well as the starting cone values and the
%    cone types replacing them. Established to be passed through aoStimRecon 
%    and used in buildRenderStruct
%    
% Inputs: 
%    chrom       - A string identifying chromaticity of built mosaic. In 
%                  the form of "chromNorm", "chromProt", "chromDeut", 
%                  "chromTrit", "chromAllL", "chromAllM", "chromAllS"
%
% Outputs: 
%   replaceCones - Boolean stating whether cones will be replaced or not
%   startCones   - Double representing value for initial cones in mosaic. 
%                  [1 2 3] corresponds to [L M S]
%   newCones     - Double representing value for desired change in mosaic. 
%                  [1 2 3] corresponds to [L M S]. Entered in the form 
%                  cMosaic.LCONE_ID, cMosaic.MCONE_ID, cMosaic.SCONE_ID
%
% Optional key/value pairs:
%    None
%
% See also: buildRenderStruct, aoStimRecon

% History
%
% 9/23/22  chr  Wrote initial version as callable function


if strcmp(chrom, 'chromProt')
    replaceCones = true;
    startCones = 1;
    newCones = cMosaic.MCONE_ID;
elseif strcmp(chrom, 'chromDeut')
    replaceCones = true;
    startCones = 2;
    newCones = cMosaic.LCONE_ID;
elseif strcmp(chrom, 'chromTrit')
    replaceCones = true;
    startCones = 3;
    newCones = cMosaic.LCONE_ID;
elseif strcmp(chrom, 'chromAllL')
    replaceCones = true;
    startCones = [2 3];
    newCones = cMosaic.LCONE_ID;
elseif strcmp(chrom, 'chromAllM')
    replaceCones = true;
    startCones = [1 3];
    newCones = cMosaic.MCONE_ID;
elseif strcmp(chrom, 'chromAllS')
    replaceCones = true; 
    startCones = [1 2]; 
    newCones = cMosaic.SCONE_ID;
elseif strcmp(chrom, 'chromNorm')
    replaceCones = false; 
    startCones = [];
    newCones = [];
else
    error('Unrecognized chromaticity input')
end

