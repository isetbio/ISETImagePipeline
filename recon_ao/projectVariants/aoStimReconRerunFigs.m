%% aoStimReconRerunFigs
%
% Description:
%    Rerun aoStimRecon to apply new updates or pull out stat info/figures.
%    Bypasses the actual reconstruction calculation.
%
% See also: aoStimRecon, aoStimReconRunMany, correctForViewing

% History:
%   06/1/23  chr  Organize into its own file


%% Set variables for the rrf sequence

% Establish the ReRunFigs struct
rrf = struct;

% State the purpose of the rerun:
%
% rerunImages: Call the full aoStimRecon script after updates
% montage: Pull recons from files for presentable montage
% dispStim: Display a montage of only the stimuli images
% statPlots: Create quantification plots for presentations
rrf.rerunImages = true;
rrf.montage = true;
rrf.dispStim = true;
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
% Assign the proper Rerun wrapper directory with the associated StimSize
% version and pull out corresponding information.
% Each subdirectory should correspond to a different stimulus size
rrf.aoReconDir = getpref('ISETImagePipeline','aoReconDir');
rrf.versEditor = 'small_quads_chr';
rrf.wrapDir = fullfile(rrf.aoReconDir , rrf.versEditor, '/StimSize_v7/Rerun');
rrf.wrapDirInfo = dir(rrf.wrapDir);


% When reading directories there are initial filler entries (., ..,
% DS_Store). The script should ignore these fillers and start with actual
% output directories, but the number changes between mac/megalodon.
if contains(rrf.aoReconDir, 'megalodon')
    firstEntry = 3;
else
    firstEntry = 3;
end


% Cycle through each stim size directory
for i = firstEntry:length(rrf.wrapDirInfo)
    if rrf.wrapDirInfo(i).isdir

        % Setup some base variables for the montage
        if (rrf.montage)

            % Establish montage dimensions and whether or not you would like to
            % disregard the first and final mosaics as extrema (i.e. if have a 100%
            % and 0% L mosaic but don't want to include due to unrealistic nature)
            numStim = 13;
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
            cellStim = cell(1,numStim);
            cellStatsStim = cell(size(cellStim));
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
                        "rgbStatsStim", "rgbStatsRecon")

                    % Pull the corrected and scaled recon from the xRunOutput,
                    % trim the edges based on the zoomLim value given above
                    % (effectively magnifying the image for presentation), and
                    % save in the desired index for the cell holding recon images
                    cellRecons{counter} = cfvRecon.reconScaledRGB{ii}(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);

                    % Pull associated Stat information and also place in a
                    % corresponding cell
                    cellStatsRecon{counter} = cfvRecon.rgbStats;

                    % Repeat the above procedure for the stimulus when the
                    % counter value falls inside the StimGrab list
                    if ismember(counter, stimGrab)
                        cellStim{1, counter/(numMosaics)} = cfvStim.stimulusRGBScaled{ii}(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);
                        cellStatsStim{1, counter/(numMosaics)} = cfvStim.rgbStats;
                    end

                    % Add to the counter
                    counter = counter + 1;

                else

                    % Otherwise load the xRunOutput.mat file from the output
                    % directory, establish the new renderDir and rerun the
                    % aoStimRecon
                    load(fullfile(rrf.outputDir, 'xRunOutput.mat'), "pr", "cnv")
                    rrf.renderDir = fullfile(rrf.aoReconDir, rrf.versEditor);
                    aoStimRecon(pr, cnv, rrf)

                end
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

            % For each stimulus, append the r/(r+g) value across different mosaics
            for q = 1:(numMosaicsTrim)
                rrgValsRecons = [rrgValsRecons cellStatsRecon{q,k}{1,5}];
            end

            % Store above values in a cell, then apply the same procedure to the
            % stimuli
            cellStats{3,k} = (rrgValsRecons);
            rrgValsStim(k) = cellStatsStim{1,k}{1,5};

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


    % Create the figure 1 legends based on input
    legend1End = repmat(" arcmin", 1, length(sizes));
    legend1Str = append(string(sizes), legend1End);
    fig1Legend = num2cell(legend1Str);

    % Create the figure 2 legends based on input
    legend2End = repmat(" %L", 1, length(mosaicSpread));
    legend2Str = append(string(mosaicSpread*100), legend2End);
    fig2Legend = num2cell(legend2Str);


    % % Stimulus quantification plots
    figure()

    % Plot stimulus r/(r+g) value vs stimulus number
    plot(fliplr(1:numStim), rrgValsStim, '-o', 'LineWidth',3); hold on

    % Plot lines based on the intersection when the stimulus looks most
    % yellow (neither red nor green)
    plot([6 6], [0 rrgValsStim(end-6+1)], '--k', 'LineWidth', 3);
    plot([0 6], [rrgValsStim(end-6+1) rrgValsStim(end-6+1)], '--k', 'LineWidth', 3);

    % Label plot and format
    xlabel('Stimulus Number', 'FontSize', 24, 'FontName', 'Arial');
    ylabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44, 'FontName', 'Arial');
    title('Stimulus Quantification', 'FontSize', 26, 'FontName', 'Arial');
    set(gcf, 'Position', [119   321   661   518]);
    axis square
    box off

    % Save output as image and eps file for easier formatting in Adobe
    % Illustrator 
    saveas(gcf,fullfile(rrf.wrapDir,'stimQuant.tiff'),'tiff');
    saveas(gcf,fullfile(rrf.wrapDir,'stimQuant.eps'),'epsc');


    % % Recon quantification plots
    % For each stimulus arcmin size
    for k = 1:size(cellStatsAll,2)
        cellStatsAll{3,k} = ones((numMosaicsTrim-3),numStim);

        % For each mosaic proportionality 
        for i = 1:(numMosaicsTrim)

            % For each stimulus presented, pull out the reconstruction
            % r/(r+g) value
            for j = 1:numStim
                cellStatsAll{3,k}(i,j) = cellStatsAll{2,k}{3,j}(i+1);
            end
        end

        figure()

        % Plot recon r/(r+g) value vs stimulus r/(r+g) value
        plot((rrgValsStim), cellStatsAll{3,k}(:, :), '-o', 'LineWidth', 3); hold on

        % Label plot and format
        xlabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
        ylabel('Reconstruction $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
        title(['Reconstruction Quantification: ' num2str(cellStatsAll{1,k}) ' arcmin'], 'FontSize', 26)
        xlim([0 1])
        set(gcf, 'Position', [119   321   661   518]);
        box off
        axis square
        plot(0:0.1:1, 0:0.1:1, '--k', 'LineWidth', 3); hold on; 
        yline(rrgValsStim(8), '--k', 'LineWidth', 3);
        legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')

        % Save output as image and eps file for easier formatting in Adobe
        % Illustrator
        saveas(gcf,fullfile(rrf.wrapDir,['reconQuant' num2str(cellStatsAll{1,k}) 'Arcmin.tiff']),'tiff');
        saveas(gcf,fullfile(rrf.wrapDir,['reconQuant' num2str(cellStatsAll{1,k}) 'Arcmin.eps']),'epsc');
    end




    % % Shift in unique yellow plot
    figure()

    % For each stimulus arcmin size
    for k = 1:1:size(cellStatsAll,2)

        % Interpolate the point at which the above line plots intersect the
        % stimulus r/(r+g) unique yellow value
        cellStatsAll{4,k} = ones(numMosaicsTrim,1);
        for i = 1:(numMosaicsTrim)
            cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), rrgValsStim, rrgValsStim(8), 'spline');
        end

        % Plot the values at which the stimulus r/(r+g) produces a Unique
        % Yellow reconstruction across each mosaic proportionality
        plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o', 'Linewidth', 3); hold on;

        % Label plot and format
        xlim([0.1 0.9]); ylim([0 1]);
        set(gca, 'XDir', 'reverse')
        xlabel('%L Cones', 'FontSize', 24);
        ylabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
        title(['Stimulus Appearance for UY Recon'], 'FontSize', 26)
    end

    % Label plot and format
    set(gcf, 'Position', [119   321   661   518]);
    box off
    axis square
    legend(fig1Legend, 'Location', 'northwest')

    % Save output as image and eps file for easier formatting in Adobe
    % Illustrator 
    saveas(gcf,fullfile(rrf.wrapDir,'shiftUY.tiff'),'tiff');
    saveas(gcf,fullfile(rrf.wrapDir,'shiftUY.eps'),'epsc');
end




