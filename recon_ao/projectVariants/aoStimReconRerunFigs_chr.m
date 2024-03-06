%% aoStimReconRerunFigs
%
% Description:
%    Rerun aoStimRecon to apply new updates or pull out stat info/figures.
%    Bypasses the actual reconstruction calculation.
%
% See also: aoStimRecon, aoStimReconRunMany, correctForViewing

% Notes:
%    Cone ratios are such that the first number in a series (i.e. 110) is 90% L and then
%    systematically goes down towards 10%L for a given chunk.
%
%    Also I can double check as for the >136 I may not have pushed those edits to the git. 137-140 are positional changes
%    in same proportionality, 141 is all L surround, 142 is all M surround

% History:
%   06/1/23  chr  Organize into its own file

%% Initialize
clear; close all;

%% Set variables for the rrf sequence

% Establish the ReRunFigs struct
rrf = struct;

% Indicate where we want to go today
%
%   rerunImages: Call the full aoStimRecon script after updates
%   montage: Pull recons from files for presentable montage
%   dispStim: Display a montage of only the stimuli images
%   statPlots: Create quantification plots for presentations
rrf.rerunImages = false;   % DHB: Setting this to true seems to be the thing that doesn't rerun the images.
rrf.montage = true;
rrf.dispStim = false;
rrf.statPlots = true;

% Select monitor display arrangements for correctForViewing.m procedure
rrf.startDisplayName = 'mono';
rrf.viewingDisplayName = 'conventional';

% Select scaling arrangements for correctForViewing
rrf.stimDispScale = 1;
rrf.reconDispScale = 1;

%% Call the montage portion

% Get preference and version editor for project
%
% Assign the proper rerun wrapper directory with the associated
% version and pull out corresponding information.
% Each subdirectory should correspond to a different stimulus size
rrf.aoReconDir = getpref('ISETImagePipeline','aoReconDir');
rrf.versEditor = 'small_quads_chr';
rrf.wrapDir = fullfile(rrf.aoReconDir , rrf.versEditor, '/StimSize_v8/Rerun');
rrf.wrapDirInfo = dir(rrf.wrapDir);

% When reading directories there are initial filler entries (., ..,
% DS_Store). The script should ignore these fillers and start with actual
% output directories, but the number changes between mac/megalodon.
if contains(rrf.aoReconDir, 'megalodon')
    firstEntry = 3;
else
    firstEntry = 4;
end

% Cycle through each stim size directory
for i = firstEntry:length(rrf.wrapDirInfo)
    if rrf.wrapDirInfo(i).isdir

        % Setup some base variables for the montage
        if (rrf.montage)

            % Establish montage dimensions and whether or not you would like to
            % disregard the first and final mosaics as extrema (i.e. if have a 100%
            % and 0% L mosaic but don't want to include due to unrealistic nature)
            numStim = 11;
            numMosaics = 9;
            trimExtrema = false;

            % Input sizes for stimuli used and progression of L proportionality
            % across mosaics tested.
            sizes = [3.5];
            mosaicSpread =  fliplr(0.1:0.1:0.9);

            % Impose some pixel limitation for the image zoom on the montage when
            % arranging for presentations
            zoomLim = 14;

            % Establish the cells that will be used to hold the all the images and
            % associated statistics for reconstructions and stimuli. Also pinpoint
            % a list of when to actually grab stimuli info instead of doing it on
            % every pass through (limit redundancy)
            cellRecons = cell(numMosaics, numStim);
            cellStatsRecon = cell(size(cellRecons));
            cellWaveRecon = cell(size(cellRecons));
            cellStim = cell(1,numStim);
            cellStatsStim = cell(size(cellStim));
            cellWaveStim = cell(size(cellStim));
            stimGrab = [1:numStim] *  (numMosaics);
        end

        % Set up counter to position images along the montage
        counter = 1;

        % Set the stim size directory name as the main directory and pull out
        % corresponding information on output directories
        rrf.mainDir = fullfile(rrf.wrapDir, rrf.wrapDirInfo(i).name);
        rrf.mainDirInfo = dir(rrf.mainDir);

        % Cycle through each output directory
        for j = firstEntry:length(rrf.mainDirInfo)

            % If the indexed mainDirInfo value is a valid directory
            if rrf.mainDirInfo(j).isdir

                % Grab the fullfile as the output directory for saving
                % information.
                rrf.outputDir = fullfile(rrf.mainDir, rrf.mainDirInfo(j).name);

                % If the goal is to create a montage
                if(rrf.montage)

                    % Load the pertinent variables from the output directory
                    load(fullfile(rrf.outputDir, 'xRunOutput.mat'), ...
                        "cfvRecon", "cfvStim", "ii", ...
                        "idxXRange", "idxYRange", ...
                        "stimulusImageLinear", "reconImageLinear")



                    % "rgbStatsStim", "rgbStatsRecon", deleted since not
                    % being used

                    cellWaveRecon{counter} = compareRenderingEW(cfvStim.stimulusRGBScaled{1}, ...
                        cfvRecon.reconScaledRGB{1}, stimulusImageLinear, reconImageLinear, ...
                        rrf.startDisplayName, rrf.viewingDisplayName, idxXRange, ...
                        'inwardMove', false, 'showFigs', false, 'scaleToMax', false);

                    % Pull the corrected and scaled recon from the xRunOutput,
                    % trim the edges based on the zoomLim value given above
                    % (effectively magnifying the image for presentation), and
                    % save in the desired index for the cell holding recon images
                    %                     cellRecons{counter} = cfvRecon.reconScaledRGB{ii}(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);
                    cellRecons{counter} = cellWaveRecon{counter}.reconImageRGB(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);




                    % Pull associated Stat information and also place in a
                    % corresponding cell
                    cellStatsRecon{counter} = cfvRecon.rgbStats;

                    % Repeat the above procedure for the stimulus when the
                    % counter value falls inside the StimGrab list
                    if ismember(counter, stimGrab)
                        %                         cellStim{1, counter/(numMosaics)} = cfvStim.stimulusRGBScaled{ii}(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);
                        cellWaveStim{1, counter/(numMosaics)} = cellWaveRecon{counter};    % Don't forget to clean this part too, just a renaming here
                        cellStim{1, counter/(numMosaics)} = cellWaveStim{1, counter/(numMosaics)}.stimImageRGB(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);

                        cellStatsStim{1, counter/(numMosaics)} = cfvStim.rgbStats;
                    end

                    % Add to the counter
                    counter = counter + 1;

                else

                    % Otherwise load the xRunOutput.mat file from the output
                    % directory, establish the new renderDir and rerun the
                end
            end
            if rrf.rerunImages
                % aoStimRecon
                load(fullfile(rrf.outputDir, 'xRunOutput.mat'), "pr", "cnv")
                rrf.renderDir = fullfile(rrf.aoReconDir, rrf.versEditor);
                aoStimRecon(pr, cnv, rrf)
            end

        end

        % Initialize a vector to hold the r/(r+g) Vals
        rrgValsStim = ones(1,numStim);

        % Remove excess mosaic proportionalities if trimming the bounds and
        % set as new variable to avoid overwriting. Otherwise, set new
        % variable equal to the full cell.
        if trimExtrema
            numMosaicsTrim = numMosaics - 2;
            cellRecons = cellRecons(2:end-1, :);
            cellStatsRecon = cellStatsRecon(2:end-1, :);
            cellWaveRecon = cellWaveRecon(2:end-1, :);
        else
            numMosaicsTrim = numMosaics;
        end

        % For each stimulus, pull out pertinent stats information for plotting
        for k = 1:numStim

            % Pull out the stored Stats cells for stimuli and recons, then
            % initialize vector
            cellStats{1,k} = cellStatsStim(k);
            cellStats{2,k} = cellStatsRecon(:,k);
            rrgValsRecons = [];

            %             % For each stimulus, append the r/(r+g) value across different mosaics
            %             for q = 1:(numMosaicsTrim)
            %                 %%% insert the call function
            % %                 avgReconLinearrgb(1,1) = {cellWaveRecon{q,k}(idxYRange,idxXRange, 1)};
            % %                 avgReconLinearrgb(2,1) = {cellWaveRecon{q,k}(idxYRange,idxXRange, 2)};
            % %                 avgReconLinearrgb(3,1) = {cellWaveRecon{q,k}(idxYRange,idxXRange, 3)};
            % %                 avgReconLinearrgb(1,2) = {mean(avgReconLinearrgb{1,1}(:))};
            % %                 avgReconLinearrgb(2,2) = {mean(avgReconLinearrgb{2,1}(:))};
            % %                 avgReconLinearrgb(3,2) = {mean(avgReconLinearrgb{3,1}(:))};
            % %                 eqWave = computeEquivWavelength([avgReconLinearrgb{1,2} ...
            % %                     avgReconLinearrgb{2,2} avgReconLinearrgb{3,2}]');
            %
            %                 eqWave = computeEquivWavelength([cellStatsRecon{q,k}{1,2} ...
            %                     cellStatsRecon{q,k}{2,2} cellStatsRecon{q,k}{3,2}]');
            %                 rrgValsRecons = [rrgValsRecons eqWave];
            %                 eqWave = [];
            % %                 rrgValsRecons = [rrgValsRecons cellStatsRecon{q,k}{1,5}];
            %             end

            % Store above values in a cell, then apply the same procedure to the
            % stimuli
            cellStats{3,k} = (rrgValsRecons);
            %             avgStimLinearrgb = squeeze(cellWaveRecon{1,k}(24, 24, :));
            %             avgStimLinearrgb = squeeze(cellStatsStim{1,k}(24, 24, :));
            %             rrgValsStim(k) = computeEquivWavelength([cellStatsStim{1,k}{1,2} ...
            %                     cellStatsStim{1,k}{2,2} cellStatsStim{1,k}{3,2}]');
            %             rrgValsStim(k) = cellStatsStim{1,k}{1,5};!!!!!

            % Calculate the mean and standard deviation across all average
            % reconstruction values. Also include the distribution of mosaic
            % proportionalities calculated above
            cellStats{4,k} = mean(rrgValsRecons);
            cellStats{5,k} = std(rrgValsRecons);
            cellStats{6,k} = mosaicSpread;
        end

        % Save cellStats as its own file in the main directory
        save(fullfile(rrf.mainDir, "cellStats.mat"), "cellStats")

        % If building the montage
        if (rrf.montage)
            figure()

            % Combine the image matrices for stimuli and recons into one big matrix
            cellFull = [cellStim; cellRecons];

            % Create a tiled montage of the images, rotate so the stimuli are
            % visualized on top (accomodating default order of imtile), and save
            figFull = imtile(cellFull, 'GridSize', [numStim, numMosaicsTrim+1]);
            imshow(imrotate(figFull, -90))
            saveas(gcf,fullfile(rrf.mainDir,'reconSlidePlot.tiff'),'tiff');
        end

        % If presenting the stimuli sequence alone
        if (rrf.dispStim)
            figure()

            % Combine the image matrices for stimuli and recons into one big matrix
            cellFull = [cellStim; cellRecons];

            % Create a tiled montage of the images, rotate so the stimuli are
            % visualized on top (accomodating default order of imtile), and save
            figFull = imtile(cellStim, 'GridSize', [numStim, 1]);
            imshow(imrotate(figFull, -90))
            saveas(gcf,fullfile(rrf.mainDir,'stimSlidePlot.tiff'),'tiff');
        end
    end
end


%% Create quantification plots


% If goal is to produce statistics plots
if (rrf.statPlots)

    % Initialize a cell with size based on number of stimuli sizes
    cellStatsAll = cell(1,length(sizes));

    % For each size
    for i = firstEntry:length(rrf.wrapDirInfo)

        % Set the name of the main directory as above
        rrf.mainDir = fullfile(rrf.wrapDir, rrf.wrapDirInfo(i).name);

        % If the indexed mainDirInfo value is a valid directory
        if isdir(rrf.mainDir)

            % Load the pertinent cellStats file
            load(fullfile(rrf.mainDir, 'cellStats.mat'));
            cellStatsAll{i-3} = cellStats;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % Reorganization since directories are read in lexicographical order.
    % Not confident this is a stable patch, should investigate this more.

    cellStatsAll = [cellStatsAll(:,end) cellStatsAll(:,1:end-1)];
    cellStatsAll = [num2cell(sizes);cellStatsAll];
    %%%%%%%%%%%%%%%%%%%%%%



    %  Begin the portion to create the plot of stimulus vs reconstruction
    %  wavelength

    % Create the figure 1 legends based on input
    legend1End = repmat(" arcmin", 1, length(sizes));
    legend1Str = append(string(sizes), legend1End);
    fig1Legend = num2cell(legend1Str);

    % Create the figure 2 legends based on input
    legend2End = repmat(" %L", 1, length(mosaicSpread));
    legend2Str = append(string(mosaicSpread*100), legend2End);
    fig2Legend = num2cell(legend2Str);

    plotColors = [1 - (0:1/(length(mosaicSpread)-1):1); ...
        (0:1/(length(mosaicSpread)-1):1); ...
        zeros(1, length(mosaicSpread))]';

    plotColorsScaled = plotColors ./ max(plotColors, [], 2);

    rrgValsStim = [];
    for w = 1:numStim
        rrgValsStim = [rrgValsStim int64(mean(cellWaveStim{w}.stimEW, 'all'))];
    end

    for k = 1:size(cellStatsAll,2)
        cellStatsAll{3,k} = ones((numMosaicsTrim-3),numStim);
        figure()

        % For each mosaic proportionality
        for i = 1:(numMosaicsTrim)

            % For each stimulus presented, pull out the reconstruction
            % r/(r+g) value
            for j = 1:numStim
                %                 cellStatsAll{3,k}(i,j) = cellStatsAll{2,k}{3,j}(i+1);
                cellStatsAll{3,k}(i,j) = int64(mean(cellWaveRecon{i,j}.reconEW, 'all'));
            end
            % Plot recon r/(r+g) value vs stimulus r/(r+g) value
            plot((rrgValsStim), cellStatsAll{3,k}(i, :), '-o', 'Color', plotColorsScaled(i,:), 'LineWidth', 3); hold on

        end

        % Label plot and format
        xlabel('Stim Wavelength', 'FontSize', 40);
        ylabel('Recon Wavelength', 'FontSize', 40);
        title(['Stim/Recon Comparison: ' num2str(cellStatsAll{1,k}) ' arcmin'], 'FontSize', 26)
        xlim([540 680])
        ylim([540 680])
        set(gcf, 'Position', [119   321   661   518]);
        box off
        axis square
        %     plot(0:0.1:1, 0:0.1:1, '--k', 'LineWidth', 3); hold on;
        yline(580, '--k', 'LineWidth', 3);
        legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')

        % Save output as image and eps file for easier formatting in Adobe
        % Illustrator
        saveas(gcf,fullfile(rrf.wrapDir,['reconQuant' num2str(cellStatsAll{1,k}) 'Arcmin.tiff']),'tiff');
        saveas(gcf,fullfile(rrf.wrapDir,['reconQuant' num2str(cellStatsAll{1,k}) 'Arcmin.eps']),'epsc');

        % Table with Wavelength Values
        allWavelengths = [rrgValsStim; cellStatsAll{3,k}];
        proportionL = ["Stimulus"; num2str(mosaicSpread')];

        fullTable = [proportionL allWavelengths];
        T = array2table(fullTable(2:end,:), 'VariableNames', fullTable(1,:));
        disp(T);

    end








  % Shift in unique yellow plot
    figure()

    % For each stimulus arcmin size
    for k = 1:1:size(cellStatsAll,2)

        % Interpolate the point at which the above line plots intersect the
        % stimulus r/(r+g) unique yellow value
        cellStatsAll{4,k} = ones(numMosaicsTrim,1);



        for i = 1:(numMosaicsTrim)
        % Approximations in the event that there are repeated wavelength 
        % values to maintain interpolation
        
            reconEW = cellStatsAll{3,k}(i,:);
%             reconEW = [1 2 3 4 5 6 7 7 9 10 11 12 13]
            reconEWUnique = unique(reconEW(:).');
    
            reconEWUniqueInd = zeros(1, length(reconEWUnique));
            for w = 1:length(reconEWUnique)
                reconEWUniqueInd(w) = find(reconEW == reconEWUnique(w), 1);
            end
            
            reconEWInterp = reconEW(reconEWUniqueInd);
            rrgValsStimInterp = double(rrgValsStim(reconEWUniqueInd));

%             reconEWUnique = unique(reconEW, 'stable')
%             reconEWUniqueMat = reconEW == reconEWUnique(:)
%             [~, reconEWUniqueInd] = find(diff(reconEWUniqueMat) == 1)

            cellStatsAll{4,k}(i,1) = interp1(reconEWInterp, rrgValsStimInterp, 580, 'spline');

            % Plot the values at which the stimulus r/(r+g) produces a Unique
            % Yellow reconstruction across each mosaic proportionality
            plot(cellStatsAll{2,1}{6,7}(i), cellStatsAll{4,k}(i,1), 'o', 'Color', plotColorsScaled(i,:), 'Linewidth', 5); hold on;
        end
        legend(fig2Legend, 'NumColumns',2, 'Location', 'northeast')
        plot(cellStatsAll{2,1}{6,7}(:), cellStatsAll{4,k}(:,1), '-k', 'Linewidth', 3, 'HandleVisibility', 'off'); hold on;


%         % Label plot and format
%         xlim([0.1 0.9]); ylim([0 1]);
%         set(gca, 'XDir', 'reverse')
        ylim([550 650])
        axis("square")
        xlabel('Local Proportion L', 'FontSize', 40);
        ylabel('Unique Yellow (nm)', 'FontSize', 40);
        title('Shift in Unique Yellow', 'FontSize', 44)
    end









    % Label plot and format
%     set(gcf, 'Position', [119   321   661   518]);
%     box off
%     axis square
%     legend(fig1Legend, 'Location', 'northwest')

%     Save output as image and eps file for easier formatting in Adobe
%     Illustrator
    saveas(gcf,fullfile(rrf.wrapDir,['shiftUY' num2str(cellStatsAll{1,k}) '.tiff']),'tiff');
    saveas(gcf,fullfile(rrf.wrapDir,['shiftUY' num2str(cellStatsAll{1,k}) '.eps']),'epsc');
end


















%
%
%
%
%

%     % % Stimulus quantification plots
%     figure()
%
%     % Plot stimulus r/(r+g) value vs stimulus number
%     plot(fliplr(1:numStim), rrgValsStim, '-o', 'LineWidth',3); hold on
%
%     % Plot lines based on the intersection when the stimulus looks most
%     % yellow (neither red nor green)
% %     plot([6 6], [0 rrgValsStim(end-6+1)], '--k', 'LineWidth', 3);
% %     plot([0 6], [rrgValsStim(end-6+1) rrgValsStim(end-6+1)], '--k', 'LineWidth', 3);
%
%     % Label plot and format
%     xlabel('Stimulus Number', 'FontSize', 24, 'FontName', 'Arial');
%     ylabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44, 'FontName', 'Arial');
%     title('Stimulus Quantification', 'FontSize', 26, 'FontName', 'Arial');
%     set(gcf, 'Position', [119   321   661   518]);
%     axis square
%     box off
%
%     % Save output as image and eps file for easier formatting in Adobe
%     % Illustrator
%     saveas(gcf,fullfile(rrf.wrapDir,'stimQuant.tiff'),'tiff');
%     saveas(gcf,fullfile(rrf.wrapDir,'stimQuant.eps'),'epsc');
%
%
%     % % Recon quantification plots
%     % For each stimulus arcmin size
%     for k = 1:size(cellStatsAll,2)
%         cellStatsAll{3,k} = ones((numMosaicsTrim-3),numStim);
%
%         % For each mosaic proportionality
%         for i = 1:(numMosaicsTrim)
%
%             % For each stimulus presented, pull out the reconstruction
%             % r/(r+g) value
%             for j = 1:numStim
% %                 cellStatsAll{3,k}(i,j) = cellStatsAll{2,k}{3,j}(i+1);
% %                 cellStatsAll{3,k}(i,j) =
%             end
%         end
%
%         figure()
%
%         % Plot recon r/(r+g) value vs stimulus r/(r+g) value
%         plot((rrgValsStim), cellStatsAll{3,k}(:, :), '-o', 'LineWidth', 3); hold on
%
%         % Label plot and format
%         xlabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
%         ylabel('Reconstruction $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
%         title(['Reconstruction Quantification: ' num2str(cellStatsAll{1,k}) ' arcmin'], 'FontSize', 26)
%         xlim([0 1])
%         set(gcf, 'Position', [119   321   661   518]);
%         box off
%         axis square
%         plot(0:0.1:1, 0:0.1:1, '--k', 'LineWidth', 3); hold on;
%         yline(rrgValsStim(8), '--k', 'LineWidth', 3);
%         legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')
%
%         % Save output as image and eps file for easier formatting in Adobe
%         % Illustrator
%         saveas(gcf,fullfile(rrf.wrapDir,['reconQuant' num2str(cellStatsAll{1,k}) 'Arcmin.tiff']),'tiff');
%         saveas(gcf,fullfile(rrf.wrapDir,['reconQuant' num2str(cellStatsAll{1,k}) 'Arcmin.eps']),'epsc');
%     end
%
%
% end
% %
% %     % % Shift in unique yellow plot
% %     figure()
% %
% %     % For each stimulus arcmin size
% %     for k = 1:1:size(cellStatsAll,2)
% %
% %         % Interpolate the point at which the above line plots intersect the
% %         % stimulus r/(r+g) unique yellow value
% %         cellStatsAll{4,k} = ones(numMosaicsTrim,1);
% %         for i = 1:(numMosaicsTrim)
% %             cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), rrgValsStim, rrgValsStim(8), 'spline');
% %         end
% %
% %         % Plot the values at which the stimulus r/(r+g) produces a Unique
% %         % Yellow reconstruction across each mosaic proportionality
% %         plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o', 'Linewidth', 3); hold on;
% %
% %         % Label plot and format
% %         xlim([0.1 0.9]); ylim([0 1]);
% %         set(gca, 'XDir', 'reverse')
% %         xlabel('%L Cones', 'FontSize', 24);
% %         ylabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
% %         title(['Stimulus Appearance for UY Recon'], 'FontSize', 26)
% %     end
% %
% %     % Label plot and format
% %     set(gcf, 'Position', [119   321   661   518]);
% %     box off
% %     axis square
% %     legend(fig1Legend, 'Location', 'northwest')
% %
% %     % Save output as image and eps file for easier formatting in Adobe
% %     % Illustrator
% %     saveas(gcf,fullfile(rrf.wrapDir,'shiftUY.tiff'),'tiff');
% %     saveas(gcf,fullfile(rrf.wrapDir,'shiftUY.eps'),'epsc');
% % end
% %
% %
% %
% %
