%% Base variables for the rrf sequence
rrf = struct;
rrf.rerunImages = true;
rrf.slidePlots = true;
rrf.statPlots = true;
rrf.trueDisplayName = 'mono';
rrf.viewingDisplayName = 'conventional';
rrf.stimDispScale = 3;
rrf.reconDispScale = 1;

%% Estabish parameters for the Stat Plots
if (rrf.slidePlots)
    rows = 2;
    colms = 13;
    mosaicSpread =  fliplr([0:2:20] / 20); %fliplr([0:8] / 8);
    sizes = [3.5 ]; % [10.5 2.5 3.5 4.5 6.5]; 
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
