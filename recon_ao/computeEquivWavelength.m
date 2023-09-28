% computeEquivWavlength
%
% Find equivalent wavelengths for various linear combinations of r and g in
% our monochromatic monitor.
%
% Code at the end shows how to get equiv wavelength 
%   1) Specified list of linear rgb triplets
%   2) For a given r/(r+g) value with b == 0
%
% Which of the above you get depends on the value of boolean COMPUTEFORRGB
%   If this is true, you get 1.  If false you get 2.
% See the relevant code sections for the specification of the values to be
% converted.
%
% For the r/(r+g) computation, there is also a table at the end with
% r/(r+g) in the first column, and equivalent wavelength in second.  For
% this option there is also example code that converts a given r/(r+g)
% value to equivalent wavelength and vice-versa.  The conversion of equiv
% wavelength to an rgb triplet is not unique, although you could add code
% to do it if you specified b.
%
% Note that r, g and b here are linear r and g values, not gamma corrected
% R, G and B values.
%
% History:
%   2023/09/26  dhb  Wrote it.
%   2023/09/28  dhb  Added rgb options.

% Setup
clear; close all;
project = 'ISETImagePipeline';
displayDir = fullfile(getpref(project,'aoReconDir'),'displays');
displayFile = fullfile(displayDir,'monoDisplay');

% When true, compute for a list of specified linear rgb values.
% When false, compute for a set of r/(r+g) on assumption that b == 0.
COMPUTEFORRGB = true;

% Work on 1 nm spacing
S = [380 1 401];
wls = SToWls(S);

% Load display and get r and g primaires out.  Spline to desired
% wavelengths.
theDisplayData = load(displayFile);
primaries = SplineSpd(displayGet(theDisplayData.monoDisplay,'wave'),displayGet(theDisplayData.monoDisplay,'spdprimaries'),S,2);

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
lowWl = 540;
highWl = 680;
lowIndex = find(wls == lowWl);
highIndex = find(wls == highWl);

% Decide which specification to use
if (COMPUTEFORRGB)
    % Make a list of linear rgb values to compute for.  You can have as
    % many columns as you like
    rgbValues = [ [0.5 0.5 0.2]' [1 0 0.2]' [0 1 0.2]' [0.5 0.5 0]' [1 0 0]' [0 1 0]'];
    nValues = size(rgbValues,2);

    for rr = 1:nValues
        theSpd = rgbValues(1,rr)*primaries(:,1) + rgbValues(2,rr)*primaries(:,2) +rgbValues(3,rr)*primaries(:,3);
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
        theEquivWls(rr) = wls(nearestWl);
        fprintf('Equivalent wavelength for rgb = %0.2f,%0.2f,%0.2f: %d\n',rgbValues(1,rr),rgbValues(2,rr),rgbValues(3,rr),theEquivWls(rr));
    end

else

    % Set up r/r+g values to convert for.  Note that
    % if we range r from 0 to 1, then set g = 1 - r,
    % then r/(r+g) = r, because r + g = 1 by construction.
    r_Over_rPlusg = linspace(0,1,100);
    rVals = r_Over_rPlusg;
    gVals = 1-rVals;

  
    % Set DEBUGFIG = true to make some debugging figures as we go.
    DEBUGFIG = false;
    if (DEBUGFIG)
        debugFig = figure;
    end

    % Loop over r_Over_rPlusg values.  For each, find the l/m chromaticy,
    % and from there the monochromatic l/m chromaticity that is closest.
    % That corresponds to the equivalent wavelength.
    for rr = 1:length(r_Over_rPlusg)
        % Get spd, L/M coordinates, and l/m chromaticity for this r_Over_rPlusg
        theSpd = rVals(rr)*primaries(:,1) + gVals(rr)*primaries(:,2);
        theLM = T_cones(1:2,:)*theSpd;
        theLMChrom = theLM/(theLM(1)+theLM(2));

        % Pop the chromaticity onto the left panel of the debug figure
        if (DEBUGFIG)
            figure(debugFig); clf;
            subplot(1,2,1); hold on;
            plot(theLMChrom(1),theLMChrom(2),'bo','MarkerSize',18,'MarkerFaceColor','b')
        end

        % Find the monochromatic l/m chromaticity that is closest by looping
        % over all of them.
        nearestDist = Inf;
        for ww = lowIndex:highIndex
            % Add to debug figure
            if (DEBUGFIG)
                subplot(1,2,1); hold on;
                plot(T_cones_LMChrom(1,ww),T_cones_LMChrom(2,ww),'ro','MarkerSize',10,'MarkerFaceColor','r');
            end

            % Get the distance beween chromaticities for this wavelength
            diff = theLMChrom - T_cones_LMChrom(:,ww);
            theDist = norm(diff);

            % Keep track of nearest
            if (theDist < nearestDist)
                nearestWl = ww;
                nearestDist = theDist;
            end

            % More debug figure stuff. Ploting distances versus wavelength in
            % right panel.
            if (DEBUGFIG)
                subplot(1,2,2); hold on;
                plot(wls(ww),theDist,'ro','MarkerSize',10,'MarkerFaceColor','r');
            end
        end

        % Get equivalent wavelength as that which had smallest distance.
        theEquivWls(rr) = wls(nearestWl);

        % Finish up debug figures
        if (DEBUGFIG)
            subplot(1,2,1); hold on;
            plot(T_cones_LMChrom(1,nearestWl),T_cones_LMChrom(2,nearestWl),'go','MarkerSize',10,'MarkerFaceColor','g');
            title(sprintf('r = %0.2f, equiv wl = %d',rVals(rr),theEquivWls(rr)));
            xlabel('L chrom'); ylabel('M chrom');

            subplot(1,2,2);
            xlabel('Wavelength (nm)'); ylabel('Distance');
            pause;
        end
    end

    % A figure with equivalent wavelenght versus r/(r+g)
    figure; clf; hold on;
    plot(r_Over_rPlusg,theEquivWls,'ro','MarkerSize',10,'MarkerFaceColor','r');
    plot(r_Over_rPlusg,theEquivWls,'r');
    xlabel('r/(r+g)');
    ylabel('Equiv Wl (nm)')

    % If you want to find equivalent wavelength for a given r/(r+g), do this.
    target_r_Over_rPlusg = 0.5;
    [~,minIndex] = min(abs(r_Over_rPlusg-target_r_Over_rPlusg));
    fprintf('Equiv wavelength matche to r/(r+g) of %0.2f: %d\n',target_r_Over_rPlusg,theEquivWls(minIndex(1)));

    % If you want to find r/(r+g) that corresponds to a target equivalent
    % wavelength, do this.
    targetEquivWl = 580;
    [~,minIndex] = min(abs(theEquivWls-targetEquivWl ));
    fprintf('r/(r+g) matched to equiv wl of %d nm: %0.2f\n',targetEquivWl,r_Over_rPlusg(minIndex(1)));

    % Table of r/(r+g) in first column, and equivalent wavelength in second
    theTable = [r_Over_rPlusg' theEquivWls'];

end





