%% Metamer/Color Percept figs, as .tif no labels

% Establish directory names and use dir function to create the loop through
% sub directories. Manually replace outputMainName here. 
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
originalOutputMain = '/AO7_AO7_0.00_0.00_58_30_conventional_conventional_dichromTestsPeanutButter';
slideList = 0;
dirStruct = dir([aoReconDir originalOutputMain]);


zoomRegions = false;


for i=4:length(dirStruct)
    outputSubName = dirStruct(i).name;
    load([fullfile(aoReconDir, originalOutputMain, outputSubName), '/xRunOutput.mat'])
    
    outputMainName = [originalOutputMain, '_figures'];
    outputDir = fullfile(aoReconDir,outputMainName,outputSubName);
    if (~exist(outputDir,'dir'))
        mkdir(outputDir);
    end
    

    if(zoomRegions)
        % Establish the zoom regions
        zoomReg = find(forwardExcitationsToStimulusUse >= 150);
        ztest = forwardConeMosaic.Mosaic.coneRFpositionsDegs(zoomReg, :);
        buffer = 0.015; 
        xmin = min(ztest(:,1)) - buffer; xmax = max(ztest(:,1)) + buffer;
        ymin = min(ztest(:,2)) - buffer; ymax = max(ztest(:,2)) + buffer;
    
        % Specific to the stimulus and reconstruction case since they're weird,
        % needed to set manually for each level for centering. 
        % [span of recon] + (slide * span of stim)
        % bounds = [-0.44 0.36] + (slide * 0.374); % 2.4
        % bounds = [-1.05 0.95] + (slide * 0.374); % 9.6
        % bounds = [] + (slide * ); % 16.8 
        % bounds = [] + (slide * ); % 24
    end


    % Visualize stimulus, remove labeling, save
    visualizeScene(stimulusScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
    set(gca, 'xtick', []); set(gca, 'ytick', []); set(gca, 'xlabel', [])
    set(gca, 'ylabel', []); set(gca, 'title', []);
    saveas(gcf,fullfile(outputDir,'stimulus.tif'),'tif');
    if(zoomRegions)
        xlim(bounds); ylim(bounds);
        saveas(gcf,fullfile(outputDir,'stimulusZoom.tif'),'tif');
    end


    % Visualize forward cone mosaic, remove labeling, save
    forwardConeMosaic.visualizeMosaic();
    set(gca, 'xtick', []); set(gca, 'ytick', []); set(gca, 'xlabel', [])
    set(gca, 'ylabel', []); set(gca, 'title', [])
    saveas(gcf,fullfile(outputDir,'forwardMosaic.tif'),'tif');
    if(zoomRegions)
        xlim([xmin xmax]); ylim([ymin ymax]);
        saveas(gcf,fullfile(outputDir,'forwardMosaicZoom.tif'),'tif');
    end


    % Visualize cone excitations, remove labeling, save
    figureHandle = figure(); axesHandle = [];
    forwardConeMosaic.Mosaic.visualize(...
    'figureHandle', figureHandle, ...
    'axesHandle', axesHandle, ...
    'activation', reshape(forwardExcitationsToStimulusUse,1,1,length(forwardExcitationsToStimulusUse)), ...
    'activationRange', [0 max(forwardExcitationsToStimulusUse)], ...
    'plotTitle',  titleStr);
    set(gca, 'xtick', []); set(gca, 'ytick', []); set(gca, 'xlabel', []);
    set(gca, 'ylabel', []); set(gca, 'title', [])
    saveas(gcf,fullfile(outputDir,'forwardMosaicExcitations.tif'),'tif');
    if(zoomRegions)
        xlim([xmin xmax]); ylim([ymin ymax]);
        saveas(gcf,fullfile(outputDir,'forwardMosaicExcitationsZoom.tif'),'tif');
    end


    % Visualize recon cone mosaic, remove labeling, save
    reconConeMosaic.visualizeMosaic();
    set(gca, 'xtick', []); set(gca, 'ytick', []); set(gca, 'xlabel', [])
    set(gca, 'ylabel', []); set(gca, 'title', []);
    saveas(gcf,fullfile(outputDir,'reconMosaic.tif'),'tif');
    if(zoomRegions)
        xlim([xmin xmax]); ylim([ymin ymax]);
        saveas(gcf,fullfile(outputDir,'reconMosaicZoom.tif'),'tif');
    end


    % Visualize reconstruction, remove labeling, save 
    visualizeScene(recon1Scene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true, 'noTitle', true);
    set(gca, 'xtick', []); set(gca, 'ytick', []); set(gca, 'xlabel', []);
    set(gca, 'ylabel', []); set(gca, 'title', []);
    saveas(gcf,fullfile(outputDir,'recon.tif'),'tif');
    if(zoomRegions)
        xlim(bounds); ylim(bounds);
        saveas(gcf,fullfile(outputDir,'reconZoom.tif'),'tif');
    end

    
    close all
end

