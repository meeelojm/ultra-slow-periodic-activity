%% ====================== POPULATION FROM SAVED DENSITIES ======================
baseFolder = '\\forskning.it.ntnu.no\ntnu\mh-kin\yaksi5\anna\Data\Processed 2P data\mGluR multimodal\cell_detect_data\';
outRoot    = fullfile(baseFolder, '_perfish_stability_outputs');

perFishFiles = dir(fullfile(outRoot, '*_periodStability.mat'));
assert(~isempty(perFishFiles), 'No per-fish periodStability files found in %s', outRoot);

% accumulators (sum + count so we can handle any NaNs robustly)
acc_d1_sum   = 0;  acc_d1_cnt   = 0;
acc_d2_sum   = {0,0,0};  acc_d2_cnt = {0,0,0};

% (optional) also keep a stack to compute medians if you’d like
stack_d1   = [];               % 50x50xN
stack_d2   = {[],[],[]};       % each 50x50xN

grid_size_expected = [50 50];
n_used = 0;

for k = 1:numel(perFishFiles)
    PF = load(fullfile(perFishFiles(k).folder, perFishFiles(k).name));
    if ~isfield(PF,'density1') || ~isfield(PF,'density2_all')
        warning('Skipping %s (missing density fields).', perFishFiles(k).name);
        continue
    end
    d1 = PF.density1;
    d2 = PF.density2_all;   % cell {P2,P3,P4}

    % basic checks/reshape if needed
    if ~isequal(size(d1), grid_size_expected)
        warning('Skipping %s (density1 not %dx%d).', perFishFiles(k).name, grid_size_expected);
        continue
    end
    ok = true;
    for p = 1:3
        if ~iscell(d2) || numel(d2) < p || ~isequal(size(d2{p}), grid_size_expected)
            warning('Skipping %s (density2_all{%d} not %dx%d).', perFishFiles(k).name, p, grid_size_expected);
            ok = false; break
        end
    end
    if ~ok, continue; end

    % accumulate mean (handle potential NaNs)
    acc_d1_sum = acc_d1_sum + nansum(d1, 3);   % (2D, nansum is fine)
    acc_d1_cnt = acc_d1_cnt + double(isfinite(d1));

    for p = 1:3
        acc_d2_sum{p} = acc_d2_sum{p} + nansum(d2{p}, 3);
        acc_d2_cnt{p} = acc_d2_cnt{p} + double(isfinite(d2{p}));
    end

    % build stacks for median if desired
    stack_d1 = cat(3, stack_d1, d1);
    for p = 1:3, stack_d2{p} = cat(3, stack_d2{p}, d2{p}); end

    n_used = n_used + 1;
end

assert(n_used > 0, 'No usable per-fish density files.');
fprintf('Included %d fish (from %d files).\n', n_used, numel(perFishFiles));

% population MEAN maps (pixel-wise, ignoring NaNs)
pop_density1_mean = acc_d1_sum ./ max(1, acc_d1_cnt);
pop_density2_mean = cell(1,3);
for p = 1:3
    pop_density2_mean{p} = acc_d2_sum{p} ./ max(1, acc_d2_cnt{p});
end
pop_density_diff_mean = {pop_density2_mean{1}-pop_density1_mean, ...
                         pop_density2_mean{2}-pop_density1_mean, ...
                         pop_density2_mean{3}-pop_density1_mean};

% (optional) population MEDIAN maps
pop_density1_median = median(stack_d1, 3, 'omitnan');
pop_density2_median = cellfun(@(S) median(S,3,'omitnan'), stack_d2, 'UniformOutput', false);
pop_density_diff_median = {pop_density2_median{1}-pop_density1_median, ...
                           pop_density2_median{2}-pop_density1_median, ...
                           pop_density2_median{3}-pop_density1_median};

% save
save(fullfile(outRoot, 'population_density_from_saved.mat'), ...
    'pop_density1_mean','pop_density2_mean','pop_density_diff_mean', ...
    'pop_density1_median','pop_density2_median','pop_density_diff_median', ...
    'n_used','-v7.3');

%% ====================== quick plots (mean) ======================
figure('Color','w','Position',[50 50 1400 600]);
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(1); imagesc(pop_density1_mean); axis image off; camroll(90); title('Pop Mean P1 (0–2m)'); colorbar;
caxis([0 1]);
for p = 1:3
    nexttile(1+p); imagesc(pop_density2_mean{p}); axis image off; camroll(90); title(sprintf('Pop Mean P%d', p+1)); colorbar;
    caxis([0 1]);
end
for p = 1:3
    nexttile(4+p); imagesc(pop_density_diff_mean{p}); axis image off; camroll(90); title(sprintf('Mean: P%d - P1', p+1)); colorbar;
    %caxis([0 ]);
end
colormap(coolwarm);

%% ====================== MSE / SSIM (population-wise) ======================
G = load(fullfile(outRoot, 'group_periodStability_summary.mat'));
MSE = [G.group.mse12,  G.group.mse13,  G.group.mse14];
SSM = [G.group.ssim12, G.group.ssim13, G.group.ssim14];
labels = {'P1vsP2','P1vsP3','P1vsP4'};

figure('Color','w','Position',[100 100 900 380]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
boxchart(repelem(1:3, size(MSE,1))', MSE(:), 'BoxFaceAlpha',0.3); hold on;
scatter(repelem(1:3, size(MSE,1)), MSE(:), 12, 'filled', 'MarkerFaceAlpha',0.6);
xlim([0.5 3.5]); xticks(1:3); xticklabels(labels); ylabel('MSE'); title('Population MSE'); grid on;

nexttile;
boxchart(repelem(1:3, size(SSM,1))', SSM(:), 'BoxFaceAlpha',0.3); hold on;
scatter(repelem(1:3, size(SSM,1)), SSM(:), 12, 'filled', 'MarkerFaceAlpha',0.6);
xlim([0.5 3.5]); xticks(1:3); xticklabels(labels); ylabel('SSIM'); title('Population SSIM'); ylim([0 1]); grid on;

% quick printouts
for j = 1:3
    mj = MSE(:,j); mj = mj(~isnan(mj));
    sj = SSM(:,j); sj = sj(~isnan(sj));
    fprintf('MSE %s: n=%d, mean=%.4g, median=%.4g\n', labels{j}, numel(mj), mean(mj), median(mj));
    fprintf('SSIM %s: n=%d, mean=%.3f, median=%.3f\n', labels{j}, numel(sj), mean(sj), median(sj));
end
