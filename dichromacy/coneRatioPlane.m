% Optimality analysis of L/M/S cone ratio across the 2D plane

% LQZ, Jan 27, 2022

%% Set up retinal mosaic
display = displayCreate('CRT12BitDisplay');
imageSize = [128, 128, 3];

%% Set up cone mosaic
% Create a 1.0 deg foveal mosaic with a pupil size of 2.5 mm
% Note that to make the tutorial run faster, reduce the size of the mosaic
retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display, 'pupilSize', 2.5, 'integrationTime', 0.1);

%% Cone ratio analysis
% generate mtxs
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0];

% figure();
for i = 1:length(s)
    for j = 1:length(l)
        lR = l(j) * (1 - s(i));
        mR = (1 - l(j)) * (1 - s(i));
        % scatter(lR, mR, '*r'); hold on;

        retina.setConeRatio(lR, mR);
        renderMtx = retina.forwardRender(imageSize, false, true, false);

        fileName = sprintf('mtx_%d_%d.mat', i, j);
        save(fullfile('.', 'allMtx', fileName), 'renderMtx', '-v7.3');

    end
end

%% After running the reconstruction
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0];

lrs = zeros(1, length(s) * length(l));
mrs = zeros(1, length(s) * length(l));

figure();
for i = 1:length(s)
    for j = 1:length(l)        
        idx = (i - 1) * length(l) + j;
        lrs(idx) = l(j) * (1 - s(i));
        mrs(idx) = (1 - l(j)) * (1 - s(i));        
    end
end

scatter(lrs, mrs, 'r*');

%% Image reconstruction error
load('input_cone_ratio.mat');
load('cone_ratio_plane.mat');
inputSet = double(input_cone_ratio);

error = zeros(1, length(s) * length(l));

for i = 1:length(s)
    for j = 1:length(l)        
        idx = (i - 1) * length(l) + j;
        reconSet = squeeze(all_recon(idx, :, :, :, :));

        rss = 0;
        for imgIdx = 1:size(inputSet, 1)
            gt = gammaCorrection(reshape(inputSet(imgIdx, :, :, :), imageSize), display);
            rc = reshape(reconSet(imgIdx, :, :, :), imageSize);

            rss = rss + norm(rc(:) - gt(:));
        end

        error(idx) = rss;
    end
end

%% Plot slices / Extend data
l_extend = [0, 0.01, 0.05:0.05:0.95, 0.99, 1.0];
error_extend = zeros(1, length(s) * length(l_extend));
for idx = 1 : length(s)
    id_start = (idx - 1) * length(l) + 1;
    sub_error = error(id_start : (id_start + length(l) - 1));

    % figure();
    % plot(l, sub_error);
    % plot(l_extend, interp1(l, sub_error, l_extend));

    sub_error = interp1(l, sub_error, l_extend);
    id_start = (idx - 1) * length(l_extend) + 1;
    error_extend(id_start : (id_start + length(l_extend) - 1)) = sub_error;
end

%% Extended data
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = l_extend;

lrs = zeros(1, length(s) * length(l));
mrs = zeros(1, length(s) * length(l));

figure();
for i = 1:length(s)
    for j = 1:length(l)        
        idx = (i - 1) * length(l) + j;
        lrs(idx) = l(j) * (1 - s(i));
        mrs(idx) = (1 - l(j)) * (1 - s(i));        
    end
end

scatter(lrs, mrs, 'r*');
error = error_extend;

%% Error plot
plotIdx = error < 1e3;

lv = 0 : 0.02 : 1;
mv = 0 : 0.02 : 1;
[L, M] = meshgrid(lv, mv);

Z = griddata(lrs(plotIdx), mrs(plotIdx), error(plotIdx), L, M, 'natural');

figure(); subplot(1, 2, 1);
surf(L, M, Z);

subplot(1, 2, 2);
contour(L, M, Z); box off;

%% contour plot
figure();

levels = [519, 520:4:535, 538, 550:50:700];
contour(L, M, Z, levels, 'ShowText', 'off'); 
hold on; box off; axis equal;

scatter(0.45, 0.45, '*r');
set(gca,'TickDir','out');
