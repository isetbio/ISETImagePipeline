function [equivWavelength,equivWavelengthCal] = RGBToEquivWavelength(inputImageRGB, theDisplay, varargin)
% Convert image RGB values to equivalent wavelength
%
% Synopsis:
%    [equivWavelength] = RGB2EquivWavelength(inputImageRGB, theDisplay)
%
% Description:
%    Input a gamma corrected RGB image, with respect to passed display.
%    Compute equivalent wavelength for each RGB triplet.  The image can be
%    just one pixel, in which case it should be of dimension 1 by 1 by 3.
%
%    The tutorial t_renderMonoDisplayImage shows how this routine is
%    called, and makes some plots that illustrate that the equivalent
%    wavelength returned is reasonable.
%
%    This works by converting the input image at each pixel to LMS cone
%    coordinates, and then to l,m chromaticity, where l = L/(L+M) and m =
%    M/(L+M).  It then finds the wavelength whose l,m coordinates are
%    closest to the l, m coordinates of each pixel, in a least squares
%    sense.
%
%    Currently the Stockman-Sharpe 2-deg fundamentals are used.  It would
%    be easy enough to add key-value pairs to use more customized (field
%    size, observer age, pupil size) fundamentals.
%
%    One could also consider doing a calculation with respect to a
%    specified white point, where one found where the line defined in chromaticity
%    from the white point through the stimulus intesected with the spectrum locus.
%    For cases where the stimuli are in the Rayleigh region and on the spectrum locus,
%    using different methods will make very little difference.  For other cases, the
%    chromaticity space and method used will matter some, but for those situations equivalent
%    wavelength itself becomes a more approximate concept and it is not really clear one 
%    wants to fuss too much.
%
%    The specification of min and max wavelength enforces that we only look
%    in the Rayleigh region of the spectrum locus.
%
% Inputs:
%    inputImageRGB     - Gamma corrected RGB input image, for theDisplay.
%                        Should be m by n by 3.
%    theDisplay        - ISETBio display structure corresponding to the input image 
%
% Outputs:
%    equivWavelength   - Equivalent wavelength corresponding to each pixel
%                        of input. in same image format as input but m by n
%                        by 1.
%    equivWavelengthCal - Equivalent wavelength in cal format, (m x m) by
%                        1.
%
% Optional key/value pairs
%    'linearInput'               - Input linear rgb rather than gamma corrected
%                                  RGB image. Default false.
%    'wls'                       - Spectral wavelength sampling, column
%                                  vector.  Default (380:1:780)'
%    'lowWl'                     - Lower bound on eq wl, nm. Default 540
%    'highWl'                    - Upper bound on eq wl, nm. Default 680
%    'verbose'                   - Print out diagnostic info. Default false.
%
% See also: t_renderMonoDisplayImage, aoStimRecon, aoStimReconRunMany

% History:
%   10/12/23  dhb  Wrote it as separate function.
      
% Parse key value pairs
p = inputParser;
p.addParameter('linearInput', false, @islogical);
p.addParameter('wls',(380:1:780)',@isnumeric);
p.addParameter('verbose',false,@islogical);
p.addParameter('lowWl',540,@isnumeric);
p.addParameter('highWl',680,@isnumeric);
parse(p, varargin{:});

% Set wavelength sampling.
wls = p.Results.wls;
S = WlsToS(wls);

% Get and spline display primaries.  Doing it way avoids NaN outside the
% range of the specified display primaries because SplineSpd allows an
% extrapolation method to be specified (see help SplineSpd), whereas if we
% set the wl sampling of the display followed by a get, we'll get back some
% NaN's in the primary spectra
displayPrimaries  = SplineSpd(displayGet(theDisplay,'wave'), ...
    displayGet(theDisplay,'spdprimaries'),S,2);

% Convert input to gamma corrected RGB if needed.
if (p.Results.linearInput)
    inputImageRGB = gammaCorrection(inputImageRGB,forwardDisplay);
end

% Get cone fundamentals.  Spline to desired wavelengths.
theConeData = load('T_cones_ss2');
T_cones = SplineCmf(theConeData.S_cones_ss2,theConeData.T_cones_ss2,S);

% Find l/m cone chromaticities for each wavelength.  We define
%    l = L/(L+M) and m = M/(L+M).
% These are normalized and thus may be compared with normalized versions
% for the various combinations of red and green primaries.
T_cones_LMChrom = T_cones(1:2,:) ./ (T_cones(1,:) + T_cones(2,:));

% We only search for equivalent wavelengths for wavelengths
% between the primaires (approx). Set these here and find
% wavelength index of each.
lowWl = p.Results.lowWl;
highWl = p.Results.highWl;
lowIndex = find(wls == lowWl);
highIndex = find(wls == highWl);

% Convert input into an ISETBio scene.  This gives us the spectrum
% at each location, although in the end we are really only using this
% routine to ungamma correct as we reconstruct the spectrum at each
% pixel below.  
meanLuminanceCdPerM2 = [];
[~, ~, inputImagergb] = sceneFromFile(inputImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, theDisplay);
if (p.Results.verbose)
    minr = min(min(inputImagergb(:,:,1)));
    ming = min(min(inputImagergb(:,:,2)));
    minb = min(min(inputImagergb(:,:,3)));
    maxr = max(max(inputImagergb(:,:,1)));
    maxg = max(max(inputImagergb(:,:,2)));
    maxb = max(max(inputImagergb(:,:,3)));
    fprintf('\trgb2EquivWavelength: min input linear: %0.4f, %0.4f, %0.4f\n',minr,ming,minb);
    fprintf('\trgb2EquivWavelength: max input linear: %0.4f, %0.4f, %0.4f\n',maxr,maxg,maxb);
end

% Get the spectrum at each pixel and compute cone coordinates.  The spectra
% are also sitting in the ISETBio scene and we could get them that way too.
[inputImagergbCal,m,n] = ImageToCalFormat(inputImagergb);
nValues = size(inputImagergbCal,2);
equivWavelengthCal = zeros(1,nValues);

% Loop over pixels.  Someone smart and motivated could probably vectorize
% the nested loops here, and make it go much faster.
for rr = 1:nValues
    % Spectrum is weighted sum of display primaries.
    theSpd = inputImagergbCal(1,rr)*displayPrimaries(:,1) + inputImagergbCal(2,rr)*displayPrimaries(:,2) +inputImagergbCal(3,rr)*displayPrimaries(:,3);
    
    % LM cone cooordinates
    theLM = T_cones(1:2,:)*theSpd;
    theLMChrom = theLM/(theLM(1)+theLM(2));

    % Find the monochromatic l/m chromaticity that is closest by looping
    % over all of them.
    nearestDist = Inf;
    for ww = lowIndex:highIndex
        % Get the distance beween chromaticities for this wavelength
        diff = theLMChrom - T_cones_LMChrom(:,ww);
        theDist = norm(diff);

        % Keep track of nearest
        if (theDist < nearestDist)
            nearestWl = ww;
            nearestDist = theDist;
        end
    end

    % Get equivalent wavelength as that which had smallest distance.
    equivWavelengthCal(rr) = wls(nearestWl);
    fprintf('Equivalent wavelength for rgb = %0.2f,%0.2f,%0.2f: %d\n',inputImagergbCal(1,rr),inputImagergbCal(2,rr),inputImagergbCal(3,rr),equivWavelengthCal(rr));
end

% Convert to image format
equivWavelength = CalFormatToImage(equivWavelengthCal,m,n);

end
