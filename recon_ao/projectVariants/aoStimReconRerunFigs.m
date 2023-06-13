%% aoStimReconRerunFigs
%
% Description:
%    Rerun aoStimRecon to apply new updates or pull out stat info/figures.
%    Bypasses the actual reconstruction calculation.
%
% See also: aoStimRecon, aoStimReconRunMany, correctForViewing

% History:
%   06/1/23  chr  Organize into its own file


%% Base variables for the rrf sequence

% Establish the ReRunFigs struct
rrf = struct;

% State the purpose of the rerun:
%
% rerunImages: Call the full aoStimRecon script after updates
% montage: Pull recons from files for presentable montage
% statPlots: Create quantification plots for presentations
rrf.rerunImages = true;
rrf.montage = false;
rrf.statPlots = false;

% Select display arrangements for correctForViewing
rrf.startDisplayName = 'mono';
rrf.viewingDisplayName = 'conventional';

% Select scaling arrangements for correctForViewing (Should be 1!!!!!!!!!!!!!!!!!!!!!!!!!!)
rrf.stimDispScale = 1;
rrf.reconDispScale = 1;


%% Call the montage portion

% Setup some base variables for the montage
if (rrf.montage)

    % Establish montage dimensions
    numStim = 13;
    numMosaics = 11;

    % Input sizes for stimuli used and progression of L proportionality
    % across mosaics tested.
    sizes = [3.5 10.5];
    mosaicSpread =  fliplr([0:2:20] / 20); %fliplr([0:8] / 8);

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

    % Create the figure 1 legends based on input
    legend1End = repmat(" arcmin", 1, length(sizes));
    legend1Str = append(string(sizes), legend1End);
    fig1Legend = num2cell(legend1Str);

    % Create the figure 2 legends based on input
    legend2End = repmat(" %L", 1, length(mosaicSpread));
    legend2Str = append(string(mosaicSpread*100), legend2End);
    fig2Legend = num2cell(legend2Str);
end


% Get preference and version editor for project
rrf.aoReconDir = getpref('ISETImagePipeline','aoReconDir');
rrf.versEditor = 'small_quads_chr';

% When reading directories there are initial filler entries (., ..,
% DS_Store). The script should ignore these fillers and start with actual
% output directories, but the number changes between mac/megalodon. 
if contains(rrf.aoReconDir, 'megalodon')
    firstEntry = 3;
else
    firstEntry = 4;
end


% Assign the proper Rerun wrapper directory with the associated StimSize
% version and pull out corresponding information. 
% Each subdirectory should correspond to a different stimulus size
rrf.wrapDir = fullfile(rrf.aoReconDir , rrf.versEditor, '/StimSize_v6/Rerun');
rrf.wrapDirInfo = dir(rrf.wrapDir);

% Cycle through each stim size directory 
for i = firstEntry:length(rrf.wrapDirInfo)

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
                    "reconRGBDispCorrectedBoost", "stimRGBDispCorrectedBoost", ...
                    "rgbStatsStim", "rgbStatsRecon")
                
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % Pull the corrected and scaled recon from the xRunOutput,
                % trim the edges based on the zoomLim value given above
                % (effectively magnifying the image for presentation), and
                % save in the desired index for the cell holding recon images
                cellRecons{counter} = reconRGBDispCorrectedBoost(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);

                % Pull associated Stat information and also place in a
                % corresponding cell
                cellStatsRecon{counter} = rgbStatsRecon;
                
                % Repeat the above procedure for the stimulus when the
                % counter value falls inside the StimGrab list
                if ismember(counter, stimGrab)
                    cellStim{1, counter/(numMosaics)} = stimRGBDispCorrectedBoost(zoomLim:end-zoomLim, zoomLim:end-zoomLim, :);
                    cellStatsStim{1, counter/(numMosaics)} = rgbStatsStim;
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
end

%% Create quantification plots 

% Initialize a vector to hold the r/(r+g) Vals 
rrgValsStim = ones(1,numStim);

% For each stimulus, pull out pertinent stats information for plotting
for k = 1:numStim

    % Pull out the stored Stats cells for stimuli and recons, then
    % initialize vector 
    cellStats{1,k} = cellStatsStim(k);
    cellStats{2,k} = cellStatsRecon(:,k);
    rrgValsRecons = [];
    
    % For each stimulus, append the r/(r+g) value across different mosaics 
    for q = 1:(numMosaics)
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
% save(fullfile(rrf.mainDir, "cellStats.mat"), "cellStats")

% If building the montage
if (rrf.montage)
    figure()
    %             cellRecons = cellRecons(2:end-1,:);

    % Combine the image matrices for stimuli and recons into one big matrix
    cellFull = [cellStim; cellRecons];
    % figFull = imtile(cellFull, 'GridSize', [numStim, (numMosaics-1)]);

    % Create a tiled montage of the images, rotate so the stimuli are
    % visualized on top (accomodating default order of imtile), and save
    figFull = imtile(cellFull, 'GridSize', [numStim, numMosaics+1]);
    imshow(imrotate(figFull, -90))
%     saveas(gcf,fullfile(rrf.mainDir,'reconSlidePlot.tiff'),'tiff');
end
%     end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If goal is to produce statistics plots
if (rrf.statPlots)

    % Initialize a cell with size based on number of stimuli sizes
    cellStatsAll = cell(1,(length(rrf.wrapDirInfo)-3));

    % For each size
    for i = 4:length(rrf.wrapDirInfo)


        rrf.mainDir = fullfile(rrf.wrapDir, rrf.wrapDirInfo(i).name);
        load(fullfile(rrf.mainDir, 'cellStats.mat'));
        cellStatsAll{i-3} = cellStats;
    end

    % Rearrange the stimulus size order to match the lexicographical
    % sorting used when reading in through the directory.
    sizesRearrange = str2double(sort(string(sizes)));
    cellStatsAll = [num2cell(sizesRearrange);cellStatsAll];


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
    %         for i = 1:(numMosaics)
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
    plot(fliplr([1:numStim]), rrgValsStim, '-o')
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
        cellStatsAll{3,k} = ones((numMosaics-3),numStim);
        for i = 1:(numMosaics)
            for j = 1:numStim
                cellStatsAll{3,k}(i,j) = cellStatsAll{2,k}{3,j}(i);
            end
        end
        plot(fliplr([1:numStim]), cellStatsAll{3,k}, '-o'); hold on
        ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
        xlabel('Stimulus Number', 'FontSize', 14);
        legend(fig2Legend(2:end-1), 'NumColumns',2)
        title(['Change in Recon Appearance: ', num2str(cellStatsAll{1,k}), ' arcmin'], 'FontSize', 16)
        set(gcf, 'Position', [119   321   661   518]);
        %             saveas(gcf,fullfile(rrf.wrapDir,['reconAppearance', num2str(cellStatsAll{1,k}), 'arcmin.tiff']),'tiff');
    end

    % For matlab plots 0.1300    0.1100    0.7750    0.8150

    % Figure for shift in UY
    figure()
    for k = 1:1:size(cellStatsAll,2)

        cellStatsAll{4,k} = ones((numMosaics-2),2);
        for i = 1:(numMosaics-2)
            cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim, rrgValsStim(7), 'spline');
            %                 cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
        end
        plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
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

























































%% TEMP COMMAND WINDOW DUMP



set(gca,'YLabel',[], 'Yticklabel',[])
set(gca,'XLabel', [], 'Xticklabel',[])
set(gcf, 'Position', [119   321   661   518]);
xlim([(2 - (10/60)/2) (2 + (10/60)/2)])
ylim([(0 - (10/60)/2) (0 + (10/60)/2)])
open correctDispImage.m
open aoStimRecon.m
aoStimReconRunMany_small_quads_chr
load('/Volumes/ExtData/Data/Aguirre-Brainard Lab Dropbox/Carlos Rodriguez/IBIO_analysis/ISETImagePipeline/aoRecon/small_quads_chr/StimSize_v6/Rerun/bckgrn10arc3.5_small_quads_chr_AO7_NOAO2_Polans2015_6_0.05_0.00_50_30_mono_1.00_conventional/3.5_2.0_0.0_0.100000_2_0.0105_0.3563_0.1550_0.0119_renderMatrix_quadSeq67_0_1_quadSeq67_0_0_13_13_nonoise_0_5/xRunOutput.mat')
aoStimReconRunMany_small_quads_chr
load('/Volumes/ExtData/Data/Aguirre-Brainard Lab Dropbox/Carlos Rodriguez/IBIO_analysis/ISETImagePipeline/aoRecon/small_quads_chr/StimSize_v6/Rerun/bckgrn10arc10_small_quads_chr_AO7_NOAO2_Polans2015_6_0.05_0.00_50_30_mono_1.00_conventional/10.0_2.0_0.0_0.100000_2_0.0105_0.3563_0.1550_0.0119_renderMatrix_quadSeq67_0_1_quadSeq67_0_0_13_13_nonoise_0_5/xRunOutput.mat')
load('/Volumes/ExtData/Data/Aguirre-Brainard Lab Dropbox/Carlos Rodriguez/IBIO_analysis/ISETImagePipeline/aoRecon/small_quads_chr/StimSize_v6/Rerun/bckgrn10arc6.5_small_quads_chr_AO7_NOAO2_Polans2015_6_0.05_0.00_50_30_mono_1.00_conventional/6.5_2.0_0.0_0.100000_2_0.0105_0.3563_0.1550_0.0119_renderMatrix_quadSeq67_0_1_quadSeq67_0_0_13_13_nonoise_0_5/xRunOutput.mat')
clear all
load('/Volumes/ExtData/Data/Aguirre-Brainard Lab Dropbox/Carlos Rodriguez/IBIO_analysis/ISETImagePipeline/aoRecon/small_quads_chr/StimSize_v6/Rerun/bckgrn10arc6.5_small_quads_chr_AO7_NOAO2_Polans2015_6_0.05_0.00_50_30_mono_1.00_conventional/6.5_2.0_0.0_0.100000_2_0.0105_0.3563_0.1550_0.0119_renderMatrix_quadSeq67_0_1_quadSeq67_0_0_13_13_nonoise_0_5/xRunOutput.mat')
load('/Volumes/ExtData/Data/Aguirre-Brainard Lab Dropbox/Carlos Rodriguez/IBIO_analysis/ISETImagePipeline/aoRecon/small_quads_chr/StimSize_v6/Rerun/bckgrn10arc10_small_quads_chr_AO7_NOAO2_Polans2015_6_0.05_0.00_50_30_mono_1.00_conventional/10.0_2.0_0.0_0.100000_2_0.0105_0.0105_0.2273_0.0119_renderMatrix_quadSeq57_0_1_quadSeq57_0_0_13_13_nonoise_0_5/xRunOutput.mat')
aoStimReconRunMany_small_quads_chr
rrf.statPlots = true;
aoStimReconRunMany_small_quads_chr
cellStatsAll = cell(1,(length(rrf.wrapDirInfo)-3));
for i = 4:length(rrf.wrapDirInfo)
rrf.mainDir = fullfile(rrf.wrapDir, rrf.wrapDirInfo(i).name);
load(fullfile(rrf.mainDir, 'cellStats.mat'));
cellStatsAll{i-3} = cellStats;
end
cellStatsAll = [num2cell(sizes);cellStatsAll];
sizes = [10.0 3.5]; % [10.5 2.5 3.5 4.5 6.5];
cellStatsAll = [num2cell(sizes);cellStatsAll];
% names starting with "10" before "2"
holder=cellStatsAll(:,1);
cellStatsAll=cellStatsAll(:,2:end);
cellStatsAll(:,end+1)=holder;
figure()
plot(fliplr([1:colms]), rrgValsStim, '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
title('Change in Stimulus Appearance', 'FontSize', 16);
set(gcf, 'Position', [119   321   661   518]);
saveas(gcf,fullfile(rrf.wrapDir,'stimAppearance.tiff'),'tiff');
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
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-1),2);
for i = 1:(rows-1)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), fliplr([1:colms]), rrgValsStim(6), 'spline');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'spline');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-1),2);
for i = 1:(rows-1)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
figure()
for k = 2:1:(size(cellStatsAll,2)-1)
cellStatsAll{4,k} = ones((rows-1),2);
for i = 1:(rows-1)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-1),2);
for i = 2:(rows-2)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
cellStatsAll{4,1}= []
cellStatsAll{4,2}= []
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-1),2);
for i = 2:(rows-2)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}(i,:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
rows-2
2:rows-2
rows
1:rows
rows-1
i
1:rows-3
cellStatsAll{2,1}{6,7}
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), fliplr([1:colms]), rrgValsStim(6), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
rrgValsStim(6)
rrgValsStim
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), fliplr([1:colms]), rrgValsStim(7), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}, cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), fliplr([1:colms]), rrgValsStim(7), 'pchip');
cellStatsAll{2,1}{6,7}
cellStatsAll{3,k}((i+1),:)
fliplr([1:colms]),
interp1(cellStatsAll{3,k}((i+1),:), fliplr([1:colms]), rrgValsStim(7), 'pchip')
interp1(fliplr([1:13])
fliplr([1:13])
rrgValsStim
interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip')
cellStatsAll{2,1}{6,7}
cellStatsAll{4,k}(:,2)
cellStatsAll{4,k}(2:(end-1),2)
cellStatsAll{4,k}(:,2)
cellStatsAll{2,1}{6,7}(2:end-1,1)
cellStatsAll{2,1}{6,7}
cellStatsAll{2,1}{6,7}(2)
cellStatsAll{2,1}{6,7}(2:end-1)
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), fliplr([1:colms]), rrgValsStim(7), 'pchip');
cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,2), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
rrgValsStim
cellStatsAll{3,k}
figure
plot(fliplr(rrgValsStim), cellStatsAll{3,k}, '-o')
figure
plot(fliplr(rrgValsStim), cellStatsAll{3,1}, '-o')
legend(fig2Legend, 'NumColumns',2)
hold on
yline(rrgValsStim(7))
fliplr([1:colms])
rrgValsStim
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim), rrgValsStim(7), 'pchip');
%                 cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim, rrgValsStim(7), 'pchip');
%                 cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
% Figure for shift in UY
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim, rrgValsStim(7), 'spline');
%                 cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
figure
plot((rrgValsStim), cellStatsAll{3,k}, '-o')
figure
plot((rrgValsStim), cellStatsAll{3,1}, '-o')
hold on
yline(rrgValsStim(7))
yline(rrgValsStim(7), '-')
yline(rrgValsStim(7), '.')
axis square
xlim([0 1])
axis square
plot(0:0.1:1, 0:0.1:1, '-')
plot(0:0.1:1, 0:0.1:1, '-k')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
ylabel('Reconstruction $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
ylabel('Reconstruction $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('Stimulus $\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel(['Stimulus $\frac{r}{r+g}$'], 'Interpreter', 'latex', 'FontSize', 20);
xlim([0 1])
axis square
plot(0:0.1:1, 0:0.1:1, '.--')
plot(0:0.1:1, 0:0.1:1, '--k')
axis square
Figure for the stimuli color progression
figure()
plot(fliplr([1:colms]), rrgValsStim, '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
title('Change in Stimulus Appearance', 'FontSize', 16);
set(gcf, 'Position', [119   321   661   518]);
%         Figure for the stimuli color progression
figure()
plot(fliplr([1:colms]), rrgValsStim, '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
title('Change in Stimulus Appearance', 'FontSize', 16);
set(gcf, 'Position', [119   321   661   518]);
axis square
figure
plot(fliplr(rrgValsStim), cellStatsAll{3,1}, '-o')
figure
plot(fliplr(rrgValsStim), cellStatsAll{3,1}{2:10, :}, '-o')
plot(fliplr(rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
legend(fig2Legend, 'NumColumns',2)
plot(fliplr(rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
set(gcf, 'Position', [119   321   661   518]);
axis square
xlim([0 1])
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
plot(fliplr(rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
xlim([0 1])
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
plot(fliplr(rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
plot((rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
figure
plot((rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
figure
plot((rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
figure
plot((rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
title(['Reconstruction Quantification: 3.5 arcmin'], 'FontSize', 26)
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
fig1Legend
fig2Legend
legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')
figure
plot((rrgValsStim), cellStatsAll{3,2}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
title(['Reconstruction Quantification: 3.5 arcmin'], 'FontSize', 26)
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim, rrgValsStim(7), 'spline');
%                 cellStatsAll{4,k}(i,2) = interp1(fliplr([1:13]), rrgValsStim, cellStatsAll{4,k}(i,1),  'pchip');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
fig1Legend = {'2.5 arcmin', '10.5 arcmin'};
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim, rrgValsStim(7), 'spline');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('%L Cones', 'FontSize', 14);
title(['Stimulus Appearance for UY Recon'], 'FontSize', 16)
end
axis square
legend(fig1Legend, 'Location', 'northwest')
set(gcf, 'Position', [119   321   661   518]);
fig1Legend = {'2.5 arcmin', '10.5 arcmin'};
figure()
for k = 1:1:size(cellStatsAll,2)
cellStatsAll{4,k} = ones((rows-3),2);
for i = 1:(rows-3)
cellStatsAll{4,k}(i,1) = interp1(cellStatsAll{3,k}((i+1),:), rrgValsStim, rrgValsStim(7), 'spline');
end
plot(cellStatsAll{2,1}{6,7}(2:end-1), cellStatsAll{4,k}(:,1), '-o'); hold on;
xlim([0.1 0.9]); ylim([0 1]);
set(gca, 'XDir', 'reverse')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('%L Cones', 'FontSize', 24);
title(['Stimulus Needed for UY across Proportions'], 'FontSize', 26)
end
axis square
legend(fig1Legend, 'Location', 'northwest')
set(gcf, 'Position', [119   321   661   518]);
box off
figure
plot((rrgValsStim), cellStatsAll{3,1}(2:10, :), '-o'); hold on
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
title(['Reconstruction Quantification: 3.5 arcmin'], 'FontSize', 26)
xlim([0 1])
set(gcf, 'Position', [119   321   661   518]);
axis square
plot(0:0.1:1, 0:0.1:1, '--k')
legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')
figure()
plot(fliplr([1:colms]), rrgValsStim, '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Stimulus Number', 'FontSize', 14);
title('Change in Stimulus Appearance', 'FontSize', 16);
set(gcf, 'Position', [119   321   661   518]);
axis square
figure()
plot(fliplr([1:colms]), rrgValsStim, '-o')
ylabel('$\frac{r}{r+g}$', 'Interpreter', 'latex', 'FontSize', 44);
xlabel('Stimulus Number', 'FontSize', 24);
title('Change in Stimulus Appearance', 'FontSize', 26);
set(gcf, 'Position', [119   321   661   518]);
axis square
box off
box off
close all