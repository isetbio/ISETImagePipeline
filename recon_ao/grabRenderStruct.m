function grabRenderStruct(renderStructure, eccXDegs, eccYDegs, ...
    fieldSizeDegs, nPixels, pupilDiamMM, aoRender, defocusDiopters)
% Synopsis:
%    Load caches files and conduct safety checks
%
% Description:
%    Used to load existing cache files of recon data as checked. After
%    loading, confirms that cached parameters match current ones. 
% 
%    This script organizes and saves ts output in the directory hierarchy
%    set up by the local hook file.
%
% See also: aoStimRecon, aoStimReconRunMany

% History:
%   08/16/22  chr  Made the check a callable function
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options


if (renderStructure.eccX ~= eccXDegs || renderStructure.eccY ~= eccYDegs)
    error('Precomputed forward rendering matrix not computed for current eccentricity');
end
if (renderStructure.fieldSizeDegs ~= fieldSizeDegs)
    error('Precomputed forward rendering matrix not computed for current field size');
end
if (renderStructure.nPixels ~= nPixels)
    error('Precomputed forward rendering nPixels not equal to current value');
end
if (renderStructure.pupilDiamMM ~= pupilDiamMM)
    error('Precompued forward pupil size not equal to current value');
end
if (renderStructure.AORender ~= aoRender)
    error('Precompued forward AO state not equal to current value');
end
if (renderStructure.defocusDiopters ~= defocusDiopters)
    error('Precompued forward defocus diopters not equal to current value');
end
end


