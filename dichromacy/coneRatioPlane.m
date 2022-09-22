% Optimality analysis of L/M/S cone ratio across the 2D plane

% LQZ, Jan 27, 2022
% Feb 3, Minor Changes

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

%% Human data
% Hofer et al., 2005
load('lm_ratio_hofer2015.mat');

%% Error surface, sparse prior
% Image reconstruction error

load('input_cone_ratio.mat');
load('../denoiser_2022/design/coneRatioRecon.mat')
inputSet = double(input_cone_ratio);

error = zeros(1, length(s) * length(l));

for i = 1:length(s)
    for j = 1:length(l)
        idx = (i - 1) * length(l) + j;

        % sparse prior set
        reconSet = allOutput{i, j};

        rss = 0;
        for imgIdx = 1:size(inputSet, 1)
            gt = gammaCorrection(reshape(inputSet(imgIdx, :, :, :), imageSize), display);
            rc = reshape(reconSet(imgIdx, :, :, :), imageSize);

            rss = rss + norm(rc(:) - gt(:));
        end

        error(idx) = rss;
    end
end

sparse_error = error;

% Plot slices / Extend data
error = sparse_error;

l_extend = [0, 0.01, 0.05:0.05:0.95, 0.99, 1.0];
error_extend = zeros(1, length(s) * length(l_extend));
for idx = 1 : length(s)
    id_start = (idx - 1) * length(l) + 1;
    sub_error = error(id_start : (id_start + length(l) - 1));

    % plot slices through the data
     figure();
     plot(l, sub_error);
    % plot(l_extend, interp1(l, sub_error, l_extend));

    sub_error = interp1(l, sub_error, l_extend);
    id_start = (idx - 1) * length(l_extend) + 1;
    error_extend(id_start : (id_start + length(l_extend) - 1)) = sub_error;
end

% Extended data
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = l_extend;

lrs = zeros(1, length(s) * length(l));
mrs = zeros(1, length(s) * length(l));

% figure();
for i = 1:length(s)
    for j = 1:length(l)
        idx = (i - 1) * length(l) + j;
        lrs(idx) = l(j) * (1 - s(i));
        mrs(idx) = (1 - l(j)) * (1 - s(i));
    end
end

% scatter(lrs, mrs, 'r*');
error = error_extend;

%% Error plot
% sparse prior
plotIdx = error < 2e3;

lv = 0 : 0.02 : 1;
mv = 0 : 0.02 : 1;
[L, M] = meshgrid(lv, mv);

Z = griddata(lrs(plotIdx), mrs(plotIdx), error(plotIdx), L, M, 'natural');

figure(); subplot(1, 2, 1);
surf(L, M, Z);

subplot(1, 2, 2);
contour(L, M, Z); box off;

% contour plot (sparse)
figure(10); subplot(1, 2, 1);
levels = [835, 840:10:890, 900:50:1600];
contour(L, M, Z, levels, 'ShowText', 'off');
hold on; box off; axis equal;

set(gca,'TickDir','out');
scatter(l_ratio, m_ratio, 50, '*r');

%% VSS Plot
figure();
levels = [0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56] + 1;
levels = min(Z(:)) * levels;
[~, c] = contour(L, M, Z, levels, 'ShowText', 'off');
c.LineWidth = 1.0;
hold on; box off; axis equal;

set(gca,'TickDir','out');
scatter(l_ratio, m_ratio, 50, '*r');

%% turn to log units
figure();
ratioZ = (Z - min(Z(:))) / min(Z(:));
logZ = log2(ratioZ);
levels = -7.5 : 0.8 : -0.5;
[~, c] = contour(L, M, logZ, levels);
c.LineWidth = 1.5;

levelLabel = (2 .^ levels) * 100;
levelLabel = arrayfun(@(x) sprintf('%.2f%%', x), levelLabel, 'UniformOutput', false);
colorbar('Ticks', levels, ...
         'TickLabels', levelLabel)

hold on; box off; axis equal;
set(gca,'TickDir','out');
scatter(l_ratio, m_ratio, 100, '*r');

xticks(0 : 0.2 : 1);
yticks(0 : 0.2 : 1);

%% Error surface, denoiser prior
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0];

lrs = zeros(1, length(s) * length(l));
mrs = zeros(1, length(s) * length(l));

for i = 1:length(s)
    for j = 1:length(l)
        idx = (i - 1) * length(l) + j;
        lrs(idx) = l(j) * (1 - s(i));
        mrs(idx) = (1 - l(j)) * (1 - s(i));
    end
end

% load input and output images
load('input_cone_ratio.mat');
load('../denoiser_2022/design/denoiseAverage.mat');
inputSet = double(input_cone_ratio);

error = zeros(1, length(s) * length(l));

for i = 1:length(s)
    for j = 1:length(l)
        idx = (i - 1) * length(l) + j;

        % denoiser set
        reconSet = squeeze(allRecon(idx, :, :, :, :));

        rss = 0;
        for imgIdx = 1:size(inputSet, 1)
            gt = gammaCorrection(reshape(inputSet(imgIdx, :, :, :), imageSize), display);
            rc = reshape(reconSet(imgIdx, :, :, :), imageSize);

            rss = rss + norm(rc(:) - gt(:));
        end

        error(idx) = rss;
    end
end

% Plot slices / Extend data
l_extend = [0, 0.01, 0.05:0.05:0.95, 0.99, 1.0];
error_extend = zeros(1, length(s) * length(l_extend));
for idx = 1 : length(s)
    id_start = (idx - 1) * length(l) + 1;
    sub_error = error(id_start : (id_start + length(l) - 1));

    % plot slices through the data
    % figure();
    % plot(l, sub_error);
    % plot(l_extend, interp1(l, sub_error, l_extend));

    sub_error = interp1(l, sub_error, l_extend);
    id_start = (idx - 1) * length(l_extend) + 1;
    error_extend(id_start : (id_start + length(l_extend) - 1)) = sub_error;
end

% Extended data
s = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
l = l_extend;

lrs = zeros(1, length(s) * length(l));
mrs = zeros(1, length(s) * length(l));

% figure();
for i = 1:length(s)
    for j = 1:length(l)
        idx = (i - 1) * length(l) + j;
        lrs(idx) = l(j) * (1 - s(i));
        mrs(idx) = (1 - l(j)) * (1 - s(i));
    end
end

% scatter(lrs, mrs, 'r*');
error = error_extend;

%% Error plot
% denoiser
plotIdx = error < 800;

lv = 0 : 0.025 : 1;
mv = 0 : 0.025 : 1;
[L, M] = meshgrid(lv, mv);

Z = griddata(lrs(plotIdx), mrs(plotIdx), error(plotIdx), L, M, 'natural');

figure(); subplot(1, 2, 1);
surf(L, M, Z);

subplot(1, 2, 2);
contour(L, M, Z); box off;

% contour plot (denoiser)
figure(10); subplot(1, 2, 2);

levels = [464:4:504, 550, 600];
contour(L, M, Z, levels, 'ShowText', 'off');
hold on; box off; axis equal;

set(gca,'TickDir','out');
scatter(l_ratio, m_ratio, '*r');

%% Plot, VSS 2022
figure();
ratioZ = (Z - min(Z(:))) / min(Z(:));
logZ = log2(ratioZ);
levels = -5.025 : 0.25 : -3;
[~, c] = contour(L, M, logZ, levels);
c.LineWidth = 1.5;

levelLabel = (2 .^ levels) * 100;
levelLabel = arrayfun(@(x) sprintf('%.1f%%', x), levelLabel, 'UniformOutput', false);
colorbar('Ticks', levels, ...
         'TickLabels', levelLabel)

hold on; box off; axis equal;
set(gca,'TickDir','out');
scatter(l_ratio, m_ratio, 100, '*r');

xticks(0 : 0.2 : 1);
yticks(0 : 0.2 : 1);
