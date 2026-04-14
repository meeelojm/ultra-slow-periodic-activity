%% 1. Find All Fish Folders with Data
baseFolder = 'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';
allFolders = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
for i = 1:length(allFolders)
    currentFolder = allFolders{i};
    if isempty(currentFolder), continue; end
    resultsFile = fullfile(currentFolder, 'dffs_repact_respcells.mat');
    metadataFile = fullfile(currentFolder, 'metadata_multimodal.mat');
    if exist(resultsFile, 'file') == 2 && exist(metadataFile, 'file') == 2
        validFolders{end+1} = currentFolder;
    end
end
disp(['N fish: ', num2str(length(validFolders))]);

%% 2. Optional: Get fps from .ini
fpss = zeros(1, length(validFolders));
for i = 1:length(validFolders)
    dataFolder = validFolders{i};
    iniFiles = dir(fullfile(dataFolder, '*.ini'));
    if ~isempty(iniFiles)
        iniContent = fileread(fullfile(dataFolder, iniFiles(1).name));
        pattern = 'volume\.rate\.\(in\.Hz\)\s*=\s*([0-9.]+)';
        tokens = regexp(iniContent, pattern, 'tokens');
        if ~isempty(tokens)
            fpss(i) = str2double(tokens{1}{1});
        else
            fpss(i) = 14.64;
        end
    else
        fpss(i) = 14.64;
    end
end

%% 3. Compute A-P and L-R PC organization for each fish
real_corr_AP = nan(1, length(validFolders));
real_corr_LR = nan(1, length(validFolders));

for iFish = 1:length(validFolders)
    dataFolder = validFolders{iFish};
    load(fullfile(dataFolder, "dffs_repact_respcells.mat"),'fps','dff_new','position'); % loads 'trace', 'position', etc.
    fps = fpss(iFish);
    dff = dff_new; % [neurons x time]
    pos = position;
    
    % Frame index mapping for 4 min window
    original_start_frame = 450; original_fps = 2.5;
    start_time_sec = original_start_frame / original_fps;
    start_frame = round(start_time_sec * fps);
    end_frame = start_frame + round(600 * fps);

    framesPerPer   = round(120 * fps);
    num_lags       = framesPerPer;    % 2-min lags
    dff_g  = dff(:, start_frame:end_frame);
    pos_g  = pos(:, 1:2);  % [X Y]
    [coeff, ~, ~] = pca(dff_g', 'NumComponents', 2); 
    angles = atan2(coeff(:,2), coeff(:,1));   
    [~, sortIdx] = sort(angles);

    % Ensure order + sorted axes (position must match dff_g subset)
    N   = numel(sortIdx);
    ord = (1:N)';
    
    sorted_y = position(sortIdx, 2);   % A-P (Y)
    sorted_x = position(sortIdx, 1);   % L-R (X)
    
    % --- 1. Real (observed) metrics ---
    real_corr_AP(iFish)  = corr(ord, sorted_y, 'type','Spearman', 'rows','complete');
    p                = polyfit(ord, sorted_y, 1);
    real_slope_AP(iFish) = p(1);
    
    real_corr_LR(iFish)  = corr(ord, sorted_x, 'type','Spearman', 'rows','complete');
    q                = polyfit(ord, sorted_x, 1);
    real_slope_LR(iFish) = q(1);
    
    % --- 2. Null Distribution (shuffle mapping between order and axis) ---
    if ~exist('n_perm','var') || isempty(n_perm), n_perm = 1000; end
    null_corrs = nan(n_perm,1); null_slopes = nan(n_perm,1);
    
    % A-P (Y)
    for pidx = 1:n_perm
        shuffled_y       = sorted_y(randperm(N));
        null_corrs(pidx) = corr(ord, shuffled_y, 'type','Spearman', 'rows','complete');
        sp               = polyfit(ord, shuffled_y, 1);
        null_slopes(pidx)= sp(1);
    end
    null_corr_AP(iFish,:)  = null_corrs;
    null_slope_AP(iFish,:) = null_slopes;
    
    % L-R (X)
    for pidx = 1:n_perm
        shuffled_x       = sorted_x(randperm(N));
        null_corrs(pidx) = corr(ord, shuffled_x, 'type','Spearman', 'rows','complete');
        sp               = polyfit(ord, shuffled_x, 1);
        null_slopes(pidx)= sp(1);
    end
    null_corr_LR(iFish,:)  = null_corrs;
    null_slope_LR(iFish,:) = null_slopes;
    
    % --- Optional: two-sided p-values from nulls (per fish) ---
    p_corr_AP(iFish)   = mean(abs(null_corr_AP(iFish,:))   >= abs(real_corr_AP(iFish)));
    p_corr_LR(iFish)   = mean(abs(null_corr_LR(iFish,:))   >= abs(real_corr_LR(iFish)));
    p_slope_AP(iFish)  = mean(abs(null_slope_AP(iFish,:))  >= abs(real_slope_AP(iFish)));
    p_slope_LR(iFish)  = mean(abs(null_slope_LR(iFish,:))  >= abs(real_slope_LR(iFish)));

end

%% 4. Plot paired comparison: A-P vs L-R
figure('Color','w','Position',[100 100 980 350]);

% (1) A-P vs null
subplot(1,3,1); hold on; box off; grid off;
v1 = abs(real_corr_AP(:));
v2 = abs(mean(null_corr_AP, 2));   % mean across permutations per fish
valid = isfinite(v1) & isfinite(v2);
A = [v1(valid) v2(valid)];
N = size(A,1);

groups = repelem(1:2, N)';          % 1..1,2..2
values = A(:);

boxchart(groups, values, 'BoxFaceAlpha',0.2);
for i = 1:N
    plot(1:2, A(i,:), '-', 'Color',[.7 .7 .7]);
end
scatter(groups, values, 40, 'filled', 'MarkerFaceAlpha',0.7);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'A-P','null'});
ylabel('|Spearman \rho|'); title('A-P vs null');

% (2) L-R vs null
subplot(1,3,2); hold on; box off; grid off;
v1 = abs(real_corr_LR(:));
v2 = abs(mean(null_corr_LR, 2));
valid = isfinite(v1) & isfinite(v2);
A = [v1(valid) v2(valid)];
N = size(A,1);

groups = repelem(1:2, N)';
values = A(:);

boxchart(groups, values, 'BoxFaceAlpha',0.2);
for i = 1:N
    plot(1:2, A(i,:), '-', 'Color',[.7 .7 .7]);
end
scatter(groups, values, 40, 'filled', 'MarkerFaceAlpha',0.7);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'L-R','null'});
ylabel('|Spearman \rho|'); title('L-R vs null');

% (3) A-P vs L-R
subplot(1,3,3); hold on; box off; grid off;
v1 = abs(real_corr_AP(:));
v2 = abs(real_corr_LR(:));
valid = isfinite(v1) & isfinite(v2);
A = [v1(valid) v2(valid)];
N = size(A,1);

groups = repelem(1:2, N)';
values = A(:);

boxchart(groups, values, 'BoxFaceAlpha',0.2);
for i = 1:N
    plot(1:2, A(i,:), '-', 'Color',[.7 .7 .7]);
end
scatter(groups, values, 40, 'filled', 'MarkerFaceAlpha',0.7);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'A-P','L-R'});
ylabel('|Spearman \rho|'); title('A-P vs L-R');


%% ===== Stats for the three panels =====
% Uses |Spearman rho|, paired within entry (fish/hemi)

% ---------- (1) A-P vs null ----------
v_real = abs(real_corr_AP(:));
v_null = abs(mean(null_corr_AP, 2, 'omitnan'));   % or median(...,2,'omitnan')
valid  = isfinite(v_real) & isfinite(v_null);
A_AP   = [v_real(valid), v_null(valid)];
N_AP   = size(A_AP,1);

[p_APnull,~,st_APnull] = signrank(A_AP(:,1), A_AP(:,2), 'method','approximate');
r_APnull = NaN; if isfield(st_APnull,'zval'), r_APnull = st_APnull.zval/sqrt(N_AP); end

med_AP_real = median(A_AP(:,1),'omitnan');  iqr_AP_real = iqr(A_AP(:,1));
med_AP_null = median(A_AP(:,2),'omitnan');  iqr_AP_null = iqr(A_AP(:,2));
dmed_AP     = med_AP_real - med_AP_null;

fprintf('\n[Stats] A-P vs null (|rho|): N=%d, p=%.4g, signedrank=%d, r=%.3f\n', ...
    N_AP, p_APnull, st_APnull.signedrank, r_APnull);
fprintf('        med real=%.3f (IQR %.3f–%.3f) | med null=%.3f (IQR %.3f–%.3f) | Δmed=%.3f\n', ...
    med_AP_real, prctile(A_AP(:,1),25), prctile(A_AP(:,1),75), ...
    med_AP_null, prctile(A_AP(:,2),25), prctile(A_AP(:,2),75), dmed_AP);

% ---------- (2) L-R vs null ----------
v_real = abs(real_corr_LR(:));
v_null = abs(mean(null_corr_LR, 2, 'omitnan'));
valid  = isfinite(v_real) & isfinite(v_null);
A_LR   = [v_real(valid), v_null(valid)];
N_LR   = size(A_LR,1);

[p_LRnull,~,st_LRnull] = signrank(A_LR(:,1), A_LR(:,2), 'method','approximate');
r_LRnull = NaN; if isfield(st_LRnull,'zval'), r_LRnull = st_LRnull.zval/sqrt(N_LR); end

med_LR_real = median(A_LR(:,1),'omitnan');  iqr_LR_real = iqr(A_LR(:,1));
med_LR_null = median(A_LR(:,2),'omitnan');  iqr_LR_null = iqr(A_LR(:,2));
dmed_LR     = med_LR_real - med_LR_null;

fprintf('\n[Stats] L-R vs null (|rho|): N=%d, p=%.4g, signedrank=%d, r=%.3f\n', ...
    N_LR, p_LRnull, st_LRnull.signedrank, r_LRnull);
fprintf('        med real=%.3f (IQR %.3f–%.3f) | med null=%.3f (IQR %.3f–%.3f) | Δmed=%.3f\n', ...
    med_LR_real, prctile(A_LR(:,1),25), prctile(A_LR(:,1),75), ...
    med_LR_null, prctile(A_LR(:,2),25), prctile(A_LR(:,2),75), dmed_LR);

% ---------- (3) A-P vs L-R ----------
v_AP = abs(real_corr_AP(:));
v_LR = abs(real_corr_LR(:));
valid = isfinite(v_AP) & isfinite(v_LR);
A_AL  = [v_AP(valid), v_LR(valid)];
N_AL  = size(A_AL,1);

[p_APvLR,~,st_APvLR] = signrank(A_AL(:,1), A_AL(:,2), 'method','approximate');
r_APvLR = NaN; if isfield(st_APvLR,'zval'), r_APvLR = st_APvLR.zval/sqrt(N_AL); end

med_AP = median(A_AL(:,1),'omitnan');  iqr_AP = iqr(A_AL(:,1));
med_LR = median(A_AL(:,2),'omitnan');  iqr_LR = iqr(A_AL(:,2));
dmed_AL = med_AP - med_LR;

fprintf('\n[Stats] A-P vs L-R (|rho|): N=%d, p=%.4g, signedrank=%d, r=%.3f\n', ...
    N_AL, p_APvLR, st_APvLR.signedrank, r_APvLR);
fprintf('        med A-P=%.3f (IQR %.3f–%.3f) | med L-R=%.3f (IQR %.3f–%.3f) | Δmed=%.3f\n\n', ...
    med_AP, prctile(A_AL(:,1),25), prctile(A_AL(:,1),75), ...
    med_LR, prctile(A_AL(:,2),25), prctile(A_AL(:,2),75), dmed_AL);

