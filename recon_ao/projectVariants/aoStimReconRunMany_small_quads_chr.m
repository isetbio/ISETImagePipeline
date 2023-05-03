%% aoStimReconMany
%
% Description:
%    Call aoStimRecon with many combinations of parameters.
%
% See also: aoStimRecon

% History:
%   08/15/22  dhb  Wrote after converting aoStimRecon to a function
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options
%   09/22/22  chr  Convert to its own dichrom file
%   09/27/22  chr  Incorporate inputs for stimulus centering position
%   10/05/22  dhb  Lots of changes for parallel
%   02/20/23  chr  Updates based on Tuten Meeting

%% Clear
clear; close all;

%% Control size of parpool, otherwise may crush memory
% thePool = gcp('nocreate');
% if (isempty(thePool))
%     parpool(5);
% end

%% Portion for ReRunning Figures from previous simulations after updates
rrf = struct;
rrf.rerunImages = true;
rrf.slidePlots = false;
rrf.statPlots = false;
rrf.trueDisplayName = 'mono';
rrf.viewingDisplayName = 'conventional';
rrf.stimDispScale = 3;
rrf.reconDispScale = 1;

if (rrf.slidePlots)
    rows = 12;
    colms = 13;
    mosaicSpread =  fliplr([0:2:20] / 20); %fliplr([0:8] / 8);
    sizes = [10.0 2.0 3.5 5.0 6.5]; % [10.5 2.5 3.5 4.5 6.5]; 
    fig1Legend = {'2.5 arcmin', '3.5 arcmin', '4.5 arcmin', '6.5 arcmin', '10.5 arcmin'};
    fig2Legend = {'90 %L', '80 %L', '70 %L', '60 %L', '50 %L', '40 %L', '30 %L', '20 %L', '10 %L'};
    % {'1.00 L:M', '0.88 L:M', '0.75 L:M', '0.63 L:M', '0.50 L:M','0.38 L:M', '0.25 L:M', '0.13 L:M', '0.00 L:M'};
    zoomLim = 14;
    cellRecons = cell(rows-1, colms);
    cellStatsRecon = cell(size(cellRecons));

    stimGrab = [1:colms] *  (rows-1);
    cellStim = cell(1,colms);
    cellStatsStim = cell(size(cellStim));
    
end

if (rrf.rerunImages)
    rrf.aoReconDir = getpref('ISETImagePipeline','aoReconDir');
    rrf.versEditor = 'small_quads_chr';

    rrf.wrapDir = fullfile(rrf.aoReconDir , rrf.versEditor, '/StimSize_v6/Rerun');
    rrf.wrapDirInfo = dir(rrf.wrapDir);

    for i = 4:length(rrf.wrapDirInfo)
        counter = 1;
        rrf.mainDir = fullfile(rrf.wrapDir, rrf.wrapDirInfo(i).name);
        rrf.mainDirInfo = dir(rrf.mainDir);

        for j = 4:length(rrf.mainDirInfo)
            if rrf.mainDirInfo(j).isdir
                rrf.outputDir = fullfile(rrf.mainDir, rrf.mainDirInfo(j).name);

                if(rrf.slidePlots)
                    load(fullfile(rrf.outputDir, 'xRunOutput.mat'), ...
                        "reconRGBDispCorrectedBoost", "stimRGBDispCorrectedBoost", ...
                        "rgbStatsStim", "rgbStatsRecon")

                    cellRecons{counter} = reconRGBDispCorrectedBoost(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);
                    cellStatsRecon{counter} = rgbStatsRecon;

                    if ismember(counter, stimGrab)
                        cellStim{1, counter/(rows-1)} = stimRGBDispCorrectedBoost(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);
                        cellStatsStim{1, counter/(rows-1)} = rgbStatsStim;
                    end

                    counter = counter + 1;

                else
                    load(fullfile(rrf.outputDir, 'xRunOutput.mat'), "pr", "cnv")
                    rrf.renderDir = fullfile(rrf.aoReconDir, rrf.versEditor);
                    aoStimRecon(pr, cnv, rrf)

                end
            end
        end

        rrgValsStim = ones(1,colms);
        for k = 1:colms
            cellStats{1,k} = cellStatsStim(k);
            cellStats{2,k} = cellStatsRecon(:,k);
            rrgValsRecons = [];

            for q = 1:(rows-1)
                rrgValsRecons = [rrgValsRecons cellStatsRecon{q,k}{1,5}];
            end
            rrgValsStim(k) = cellStatsStim{1,k}{1,5};
            cellStats{3,k} = (rrgValsRecons);
            cellStats{4,k} = mean(rrgValsRecons);
            cellStats{5,k} = std(rrgValsRecons);
            cellStats{6,k} = mosaicSpread;
        end


        save(fullfile(rrf.mainDir, "cellStats.mat"), "cellStats")

        
        if (rrf.slidePlots)
            figure()
%             cellRecons = cellRecons(2:end-1,:);
            cellFull = [cellStim; cellRecons];
%             figFull = imtile(cellFull, 'GridSize', [colms, (rows-2)]);
            figFull = imtile(cellFull, 'GridSize', [colms, rows]);
            imshow(imrotate(figFull, -90))
            saveas(gcf,fullfile(rrf.mainDir,'reconSlidePlot.tiff'),'tiff');
        end
    end

    if (rrf.statPlots)
        cellStatsAll = cell(1,(length(rrf.wrapDirInfo)-3));
        for i = 4:length(rrf.wrapDirInfo)
            rrf.mainDir = fullfile(rrf.wrapDir, rrf.wrapDirInfo(i).name);
            load(fullfile(rrf.mainDir, 'cellStats.mat'));
            cellStatsAll{i-3} = cellStats;
        end
        
        cellStatsAll = [num2cell(sizes);cellStatsAll];


        % Temporary reorganizing since the above cycles through the directory
        % names starting with "10" before "2"
        holder=cellStatsAll(:,1);
        cellStatsAll=cellStatsAll(:,2:end);
        cellStatsAll(:,end+1)=holder;


%         % Figure 1
%         figure()
%         subplot(3,1,1)
%         for j=1:length(cellStatsAll)
%             plot(cellStatsAll{2,j}{6,7}, cellStatsAll{2,j}{3,7}, '-o'); hold on;
%         end
%         set(gca, 'XDir', 'reverse')
%         xlabel('L:M Ratio', 'FontSize', 14);
%         ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
%         title('Recon Primary Proportions over Cone Proportions', 'FontSize', 16);
%         legend(fig1Legend);
% 
%         % Figure 2
%         subplot(3,1,2)
%         for i = 1:(rows-1)
%             tempStore = ones(1,length(cellStatsAll));
%             for j=1:length(cellStatsAll)
%                 tempStore(j) = cellStatsAll{2,j}{3,7}(i);
%             end
%             plot(cell2mat(cellStatsAll(1,:)), tempStore, '-o'); hold on;
%         end
%         xlabel('Stimulus Size (arcmin)', 'FontSize', 14);
%         ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
%         title('Recon Primary Proportions over Stim Size',  'FontSize', 16);
%         legend(fig2Legend, 'NumColumns', 2);
% 
%         % Figure 3
%         subplot(3,1,3)
%         tempStore = ones(1,length(cellStatsAll));
%         for j=1:length(cellStatsAll)
%             tempStore(j) = cellStatsAll{2,j}{5,7};
%         end
%         plot(cell2mat(cellStatsAll(1,:)), tempStore, '-o'); hold on;
%         xlabel('Stimulus Size (arcmin)', 'FontSize', 14);
%         ylabel('Yellow Recon Variation', 'FontSize', 14);
%         title('Recon Variation over Stim Size', 'FontSize', 16);
%         set(gcf, 'Position', [335     1   564   976]);
%         saveas(gcf,fullfile(rrf.wrapDir,'sumStatPlotsYellow.tiff'),'tiff')



%         Figure for the stimuli color progression
        figure()
        plot(fliplr([1:colms]), rrgValsStim, '-o')
        ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
        xlabel('Stimulus Number', 'FontSize', 14);
        title('Change in Stimulus Appearance', 'FontSize', 16);
        set(gcf, 'Position', [119   321   661   518]);
        saveas(gcf,fullfile(rrf.wrapDir,'stimAppearance.tiff'),'tiff');

        % For this, UY corresponds to the stimulus at index 6 (looks
        % yellow)

        % Figure for the recon color progression
        for k = 1:size(cellStatsAll,2)
            figure()
            cellStatsAll{3,k} = ones((rows-4),colms);
            for i = 1:(rows-1)
                for j = 1:colms
                    cellStatsAll{3,k}(i,j) = cellStatsAll{2,k}{3,j}(i);
                end
            end
            plot(fliplr([1:colms]), cellStatsAll{3,k}, '-o'); hold on
            ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
            xlabel('Stimulus Number', 'FontSize', 14);
            legend(fig2Legend, 'NumColumns',2)   
            title(['Change in Recon Appearance: ', num2str(cellStatsAll{1,k}), ' arcmin'], 'FontSize', 16)
            set(gcf, 'Position', [119   321   661   518]);
%             saveas(gcf,fullfile(rrf.wrapDir,['reconAppearance', num2str(cellStatsAll{1,k}), 'arcmin.tiff']),'tiff');
        end
        
        % For matlab plots 0.1300    0.1100    0.7750    0.8150

        % Figure for shift in UY
        figure()
        for k = 1:1:size(cellStatsAll,2)

            cellStatsAll{4,k} = ones((rows-1),2);
            for i = 1:(rows-1)
                cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
                cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'spline');
            end
            plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
            xlim([0.1 0.9]); ylim([0 1]);
            set(gca, 'XDir', 'reverse')
            ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
            xlabel('%L Cones', 'FontSize', 14);
            title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
        end

        legend(fig1Legend, 'Location', 'southeast')
        set(gcf, 'Position', [119   321   661   518]);
        saveas(gcf,fullfile(rrf.wrapDir,['reconUYStim.tiff']),'tiff');
        %             usableLims = cellStatsAll{4,k}(:,2) > 0 & cellStatsAll{4,k}(:,2) < 1;
        %             usableLims(1) = false; usableLims(end) = false;
        %             dummyVar = cellStatsAll{4,k}(:,2);
        %             usablePoints = dummyVar(usableLims);
        %             slopeVals = polyfit(1:length(usablePoints), usablePoints, 1);
        %             cellStatsAll{5,k} = slopeVals(1);
        %             cellStatsAll{6,k} = std(usablePoints);




        %             y2 = polyval(slopeVals,cellStatsAll{2,1}{6,7});
        %             plot(cellStatsAll{2,1}{6,7}, y2);
        %             dim = [.2 .5 .3 .3];
        %             linFit = -1 * (polyfit(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), 1));
        %             str = ['Linear slope: ', num2str(linFit(1))];
        %             annotation('textbox',dim,'String',str,'FitBoxToText','on');


        %             set(gcf, 'Position', [119   321   661   518]);
        %             saveas(gcf,fullfile(rrf.wrapDir,['reconAppearance', num2str(cellStatsAll{1,k}), 'arcmin.tiff']),'tiff');
        save(fullfile(rrf.wrapDir, "cellStatsAll.mat"), "cellStatsAll")
    end





elseif ~(rrf.rerunImages)
    %% Set defaults in prBase
    prBase = prBaseDefaults;

    %% Version editor string
    %
    % Helps us keep different calcs separate
    prBase.versEditor = 'small_quads_chr';

    %% Parameters
    %
    % Display, options are:
    %    'conventional'    - A conventional display
    %    'mono'            - A display with monochromatic primaries
    prBase.displayName = 'mono';
    prBase.displayGammaBits = 12;
    prBase.displayGammaGamma = 2;
    displayScaleFactorList = [1];

    %% Spatial parameters
    %
    % Common to forward and recon models
    prBase.nPixels = 50;
    prBase.trueCenter = round(prBase.nPixels/2);

    %% Mosaic parameters
    prBase.fieldSizeMinutes = 30;
    prBase.eccXDegs = 2.0;
    prBase.eccYDegs = 0.0;
    prBase.forwardRandSeed = false;
    prBase.reconRandSeed = false;
    prBase.forwardEccVars = false;
    prBase.reconEccVars = false;
    prBase.reconstructfromRenderMatrix = true;
    prBase.addPoissonNoise = false;

    % Mosaic chromatic type, options are:
    %    "chromNorm", "chromProt", "chromDeut", "chromTrit",
    %    "chromAllL", "chromAllM", "chromAllS", "quadSeq" and number
    %    Currently established quadSeq1 - quadSeq56
    forwardChromList = ["quadSeq78" "quadSeq78" "quadSeq79" "quadSeq80" "quadSeq81" "quadSeq82" "quadSeq83" "quadSeq84" "quadSeq85" "quadSeq86" "quadSeq87"]; % Don't forget to run QS34 on 4@0.5
    reconChromList =   ["quadSeq78" "quadSeq78" "quadSeq79" "quadSeq80" "quadSeq81" "quadSeq82" "quadSeq83" "quadSeq84" "quadSeq85" "quadSeq86" "quadSeq87"]; % 36, 38, 40, 42, 44

    % Build new sequence by
    prBase.quads(1).name = 'useQuadSeq';
    prBase.quads(1).value = true;

    % If want to apply percentages to the full mosaic instead of a
    % quadrant, set to true and use percentages corresponding to Quadrant 4
    % (prBase.quads(5))
    fullMosaicPercent = true;

    if(prBase.quads(1).value)
        % Initialize storage structure with information on each quadrant
        prBase.quads(2).name = 'Quad1';
        prBase.quads(3).name = 'Quad2';
        prBase.quads(4).name = 'Quad3';
        prBase.quads(5).name = 'Quad4';

        % Enter desired percent as decimal of L cones per region across
        % quadrants. The remaining percent will be made of M cones. Entries
        % should start with outermost regions first and progress inward
        prBase.quads(2).percentL = [0.27];
        prBase.quads(3).percentL = [0.53];
        prBase.quads(4).percentL = [0.71];
        prBase.quads(5).percentL = [0.70];

        % Enter desired percent as decimal of S cones per region across
        % quadrants. Follows same form as above
        prBase.quads(2).percentS = [0.05];
        prBase.quads(3).percentS = [0.05];
        prBase.quads(4).percentS = [0.05];
        prBase.quads(5).percentS = [0.1];

        % Establish initial region boundaries in the x and y direction for all
        % four quadrants based on FOV
        prBase.quads(2).xbounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccXDegs;
        prBase.quads(3).xbounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccXDegs;
        prBase.quads(4).xbounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccXDegs;
        prBase.quads(5).xbounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccXDegs;
        prBase.quads(2).ybounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccYDegs;
        prBase.quads(3).ybounds = [0 prBase.fieldSizeMinutes/60/2] + prBase.eccYDegs;
        prBase.quads(4).ybounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccYDegs;
        prBase.quads(5).ybounds = [-prBase.fieldSizeMinutes/60/2 0] + prBase.eccYDegs;

        if (fullMosaicPercent)
            prBase.quads(5).xbounds = [-prBase.fieldSizeMinutes/60/2 prBase.fieldSizeMinutes/60/2] + prBase.eccXDegs;
            prBase.quads(5).ybounds = [-prBase.fieldSizeMinutes/60/2 prBase.fieldSizeMinutes/60/2] + prBase.eccYDegs;
        end
    end

    prBase.quads(6).name = 'overrideQuadSeq';
    prBase.quads(6).value = true;

    % Add indices of cones to be silenced.
    prBase.kConeIndices = [];

    % Select which quadrants from the above to activate, of the order:
    % Q1 Q2 Q3 Q4 Origin
    quadSelectList = [[false false false false true]]';%...

    % Force build and save of render structures.  This
    % only affects this script, and will typically be false.
    buildNewForward = false;
    buildNewRecon = false;

    %% Stimulus parameters.
    %
    % Size list parameter in degs, expressed as min/60 (because 60 min/deg)
    stimSizeDegsList = [2.0 3.5 5.0 6.5 10.5] / 60;

    % RGB values (before gamma correction)
    prBase.stimBgVal = 1;% [0.1054 0.1832 0.1189]
    stimRValList = [1];% 1.0 0.0];
    stimGValList = [1];% 0.0 1.0];
    stimBValList = [0];% 0.0 0.0];

    % Overwrite stim values to make isoluminant colors on the RG channel.
    % Allow for offset from the true isoLum values based on variability in
    % staircase procedure where (-) is more red and (+) is more green,
    % (Intervals of 50s? 100s? 1000s?)
    isoLumRG = true;
    colorStepRG = [-1729 -1280 -880 -480 0 240 480 720 960 1200 1440 1680 1729];

    if (isoLumRG)
        % Load the appropriate display
        theDisplayLoad = load(fullfile(prBase.aoReconDir, 'displays', [prBase.displayName 'Display.mat']));
        switch (prBase.displayName)
            case 'conventional'
                displayFieldName = 'CRT12BitDisplay';
                prBase.stimBgVal = [0.1054 0.1832 0.1189]/10;
            case 'mono'
                displayFieldName = 'monoDisplay';
                prBase.stimBgVal = [0.1054 0.1832 0.1189]/10;
            otherwise
                error('Unknown display specified');
        end

        % Grab the display primary xyz values and pull out luminance column
        eval(['primariesXYZ = displayGet(theDisplayLoad.' displayFieldName ', ''primaries xyz'');']);
        lumRGB = primariesXYZ(:,2);

        % Create offset vectors for R and G-raw values, then multiply raw value
        % by luminance ratio between channels. Find index where R and G
        % channels are most similar and add step if desired
        isoLumR = fliplr(0:0.0001:1); isoLumGRaw = 1 - isoLumR;
        ratioRG = lumRGB(1) / lumRGB(2); alpha = lumRGB(2) / lumRGB(1);
        isoLumG = ratioRG .* isoLumGRaw;

        gMax1 = (1-prBase.stimBgVal(1))/alpha; rMax1 = gMax1 * alpha;
        gMax2 = (1-prBase.stimBgVal(2)); rMax2 = gMax2 * alpha;

        if gMax1 > 1 || rMax1 > 1
            rBound = find(isoLumR < rMax2);
            isoLumR = isoLumR(min(rBound):end);
            isoLumG = isoLumG(min(rBound):end);
        elseif gMax2 > 1 || rMax2 > 1
            rBound = find(isoLumR < rMax1);
            isoLumR = isoLumR(min(rBound):end);
            isoLumG = isoLumG(min(rBound):end);
        else
            error('Both calculations with background exceed limits')
        end

        difRG = abs(isoLumR - isoLumG)';
        indRG = find(difRG == min(difRG)) + colorStepRG;

        % Overwrite R and G val list based on isolum conditions. To all three
        % channels add background stim value.
        stimRValList = isoLumR(indRG) + prBase.stimBgVal(1);
        stimGValList = isoLumG(indRG) + prBase.stimBgVal(2);
        stimBValList = zeros(1,length(stimRValList)) + prBase.stimBgVal(3);

        % Clean workspace
        clear theDisplayLoad; clear displayFieldName; clear primariesXYZ;
        clear isoLumR; clear isoLumG; clear isoLumGRaw; clear difRG; clear indRG;
    end

    % Check that all channels receive same number of inputs
    if (length(stimGValList) ~= length(stimRValList) || length(stimBValList) ~= length(stimRValList))
        error('Stimulus value lists must have same length');
    end

    %% Positions
    %
    % Input desired x and y position for stimulus to be centered over. Function
    % will end if values exceed pixel limits.
    %
    % Position specified in pixels, could consider specifying in minutes.
    pixelsPerMinute = prBase.nPixels/prBase.fieldSizeMinutes;
    shiftInMinutesListX = [0];
    shiftInMinutesListY = [0];
    fullSquareShift = false;

    % Convert the shifts to pixel positions
    shiftInPixelsListX = round(pixelsPerMinute*shiftInMinutesListX);
    shiftInPixelsListY = round(pixelsPerMinute*shiftInMinutesListY);

    % Consider including the term "-(prBase.nPixels * 0.1)" to better center
    % the quad stim at larger values so not pushing against the outer edge.
    quadCenters = round((prBase.nPixels) / 4);
    centerXPosition = prBase.trueCenter + quadCenters + shiftInPixelsListX;
    centerYPosition = prBase.trueCenter + quadCenters + shiftInPixelsListY;
    prBase.stimCenter = [centerXPosition ; centerYPosition];

    % Loop through created pixel positions if want to create a square grid of
    % movement instead of default horizontal shift
    if (fullSquareShift)
        centerXPosition = repelem(prBase.stimCenter(1,:), length(shiftInMinutesList));
        centerYPosition = repmat(prBase.stimCenter(1,:), [1,length(shiftInMinutesList)]);
        prBase.stimCenter = [centerXPosition; centerYPosition];
    end
    deltaCenterList = [prBase.stimCenter - prBase.trueCenter];

    %% Prior parameters
    %
    % conventionalSparsePrior - from the paper, images analyzed on conventional display.
    prBase.sparsePriorStr = 'conventional';

    %% Reconstruction parameters
    %
    % Should cycle through a few of these regs to optimize for 58x58 pixels
    % Previous pairs: 100x100 at 5e-3, 128x128 at 1e-2
    regParaList = 0.1;
    prBase.stride = 2;
    prBase.maxReconIterations = 2000;
    prBase.whiteNoiseStarts = 0;
    prBase.pinkNoiseStarts = 1;
    prBase.sparsePriorPatchStarts = 0;
    prBase.stimulusStart = false;
    % Should note, the line below only works if have >= 3 starting fields
    prBase.uniformStartVals = []; % [[0.5 0.5 0.5]' [0 0.5 0]' [0 0 0.5]'];
    prBase.boundedSearch = false;

    % Use AO in forward rendering? And determine optics pupil size
    prBase.forwardAORender = true;
    prBase.forwardNoLCA = true;
    prBase.reconAORender = false;
    prBase.reconNoLCA = false;
    reconPupilDiamListMM =  [2];
    forwardPupilDiamListMM = [7];

    % Define optics.  Subject only matters if we use a database.  Ignored for
    % Marimont and Wandell.  For database, subjectID of 0 means diffraction
    % limited.
    %
    % Databases are 'Polans2015' and 'Artal2012'
    prBase.forwardSubjectID = 6;
    prBase.forwardZernikeDataBase = 'Polans2015';
    prBase.reconSubjectID = 6;
    prBase.reconZernikeDataBase = 'Polans2015';

    % Residual defocus for forward and recon rendering, of equal sizes
    forwardDefocusDioptersList = [0.05];
    reconDefocusDioptersList = [0.0];

    %% Set up list conditions
    runIndex = 1;
    for ss = 1:length(stimSizeDegsList)
        for cc = 1:length(stimRValList)
            for yy = 1:size(deltaCenterList,2)
                for ff = 1:length(forwardDefocusDioptersList)
                    for rr = 1:length(regParaList)
                        for dd = 1:length(forwardChromList)
                            for pp = 1:length(forwardPupilDiamListMM)
                                for qq = 1:length(quadSelectList(1,:))
                                    for dsf = 1:length(displayScaleFactorList)

                                        stimSizeDegs(runIndex) = stimSizeDegsList(ss);

                                        stimRVal(runIndex) = stimRValList(cc);
                                        stimGVal(runIndex) = stimGValList(cc);
                                        stimBVal(runIndex) = stimBValList(cc);

                                        stimCenter(:,runIndex) = deltaCenterList(:,yy);

                                        forwardDefocusDiopters(runIndex) = forwardDefocusDioptersList(ff);
                                        reconDefocusDiopters(runIndex) = reconDefocusDioptersList(ff);

                                        regPara(runIndex) = regParaList(rr);

                                        forwardChrom(runIndex) = forwardChromList(dd);
                                        reconChrom(runIndex) = reconChromList(dd);

                                        forwardPupilDiamMM(runIndex) = forwardPupilDiamListMM(pp);
                                        reconPupilDiamMM(runIndex) = reconPupilDiamListMM (pp);

                                        quadSelect(:,runIndex) = quadSelectList(:,qq);

                                        displayScaleFactor(runIndex) = displayScaleFactorList(dsf);

                                        runIndex = runIndex + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    %% Build render structures we need if they are not cached
    for pp = 1:length(regPara)

        % Set up paramters structure for this loop, filling in fields that come
        % out of lists precreated above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor);
        pr.quadSelect = quadSelect(:,pp);
        % Compute convenience parameters
        cnv = computeConvenienceParams(pr);

        % Build foward cone mosaic and render matrix if needed
        if (buildNewForward || ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'file'))
            renderStructure = buildRenderStruct(pr, cnv, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardNoLCA, ...
                pr.forwardDefocusDiopters, pr.forwardRandSeed, cnv.replaceForwardCones, cnv.forwardStartCones, ...
                cnv.forwardNewCones, pr.forwardEccVars, pr.forwardSubjectID, pr.forwardZernikeDataBase, pr.forwardChrom);
            save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'renderStructure','-v7.3');
            forwardRenderStructure = renderStructure; clear renderStructure;
        end

        % Build recon cone mosaic and render structure if needed
        if (buildNewRecon || ~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'file'))
            renderStructure = buildRenderStruct(pr, cnv, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconNoLCA, ...
                pr.reconDefocusDiopters, pr.reconRandSeed, cnv.replaceReconCones, cnv.reconStartCones, ...
                cnv.reconNewCones, pr.reconEccVars, pr.reconSubjectID, pr.reconZernikeDataBase, pr.reconChrom);
            save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'renderStructure','-v7.3');
            reconRenderStructure = renderStructure; clear renderStructure;
        end
    end

    % THIS SHOULD BE A PARFOR AFTERWARDS DON'T FORGET
    parfor pp = 1:length(regPara)

        % Set up paramters structure for this loop, filling in fields that come
        % out of lists above.
        pr = prFromBase(prBase,pp,stimSizeDegs,stimRVal,stimGVal,stimBVal, ...
            stimCenter,forwardDefocusDiopters,reconDefocusDiopters,regPara, ...
            forwardChrom,reconChrom,forwardPupilDiamMM,reconPupilDiamMM,displayScaleFactor);
        pr.quadSelect = quadSelect(:,pp);

        % Compute convenience parameters
        cnv = computeConvenienceParams(pr);

        % Call the driving function
        aoStimRecon(pr,cnv, rrf);
    end
end
