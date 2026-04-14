%% ===================== CCN2A: Corr vs Distance (ALL CELLS) =====================
% Requires: corr_vs_distance_AO.m, shadedErrorBar.m

%% ---------- CONFIG ----------
baseFolder     = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\5\Emiliano\ccn2a\light_tap_free_tail';
origStartFrame = 200;          % where your analysis window starts in the raw movie
origFps        = 2.5;          % fps used for origStartFrame
secPerPeriod   = 120;          % concatenate 4 x 120 s (P1..P4) = 8 min
please_split   = 1;            % split hemispheres
side2split     = 2;            % split by Y
nRand          = 200;          % position-shuffle control per fish
rng(0);

% colors
blue = [0 0.4470 0.7410];

% optional saving
outRoot = fullfile(baseFolder, '_corr_vs_distance_outputs');
if ~exist(outRoot,'dir'), mkdir(outRoot); end

%% ---------- FIND FISH FOLDERS WITH DATA ----------
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
resultsFiles = {};
metaFiles    = {};
fishLabel    = {};
for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    rf = fullfile(allFolders{i}, 'dffs_repact_respcells.mat');
    mf = fullfile(allFolders{i}, 'metadata_multimodal.mat');
    if exist(rf,'file')==2 && exist(mf,'file')==2
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
        resultsFiles{end+1} = rf;            %#ok<SAGROW>
        metaFiles{end+1}    = mf;            %#ok<SAGROW>
        [~, fishLabel{end+1}] = fileparts(allFolders{i}); %#ok<SAGROW>
    end
end
fprintf('N fish folders: %d\n', numel(validFolders));

%% ---------- FPS PER FISH ----------
fps_by_fish = nan(1, numel(validFolders));
for i = 1:numel(validFolders)
    fps = try_get_fps_from_metadata(metaFiles{i});
    if ~isfinite(fps)
        iniFiles = dir(fullfile(validFolders{i}, '*.ini'));
        if ~isempty(iniFiles)
            fps = try_get_fps_from_ini(fullfile(iniFiles(1).folder, iniFiles(1).name));
        end
    end
    if ~isfinite(fps), fps = 14.64; end     % sensible default for this dataset
    fps_by_fish(i) = fps;
end
%% ---------- ACCUMULATORS & FIXED DISTANCE GRID ----------
% distance axis (same for every fish)
binStep_um   = 5;                      % resolution (µm)
maxRange_um  = 100;                    % x-axis limit (µm)
F_binEdges   = 0:binStep_um:maxRange_um;
F_binCenters = F_binEdges(1:end-1) + diff(F_binEdges)/2;
nBins        = numel(F_binCenters);

% per-fish stores (rows = fish, cols = distance bins)
F_curve_perFish      = zeros(0, nBins);   % all-cells corr vs distance
F_ctrl_mean_perFish  = zeros(0, nBins);   % shuffle mean
F_ctrl_sem_perFish   = zeros(0, nBins);   % shuffle SEM
F_pairs_perFish      = zeros(0, nBins);   % # of pairs contributing per bin
F_keptFish           = strings(0,1);      % fish IDs
F_nCells             = [];                % # cells per fish

%% ---------- PER-FISH LOOP (robust, no brace arithmetic) ----------
nBins = numel(F_binCenters);
for f = 1:numel(resultsFiles)
    fid = fishLabel{f};
    fps = fps_by_fish(f);
    fprintf('[%d/%d] %s (fps=%.3g)\n', f, numel(resultsFiles), fid, fps);

    S = load(resultsFiles{f});

    % ---- dFF [cells x frames] ----
    if isfield(S,'dff_new')
        dff = double(S.dff_new);
    else
        warning('  %s: dff_new not found → skip.', fid);
        continue;
    end
    if size(dff,1) > size(dff,2), dff = dff.'; end
    Ncells = size(dff,1);

    % ---- positions [cells x 2/3], convert to µm (÷10) ----
    pos = [];
    cand = {'positions_nodoubles','position','xyz','XYZ','coords','centroids'};
    for c = 1:numel(cand)
        if isfield(S, cand{c}) && ~isempty(S.(cand{c}))
            pos = double(S.(cand{c}));
            break;
        end
    end
    if isempty(pos)
        warning('  %s: positions not found → skip.', fid);
        continue;
    end
    if size(pos,1) ~= Ncells
        if size(pos,2) == Ncells
            pos = pos.';
        else
            warning('  %s: positions size mismatch → skip.', fid);
            continue;
        end
    end
    if size(pos,2) < 2
        warning('  %s: positions need >=2 columns → skip.', fid);
        continue;
    end
    if size(pos,2) > 3, pos = pos(:,1:3); end
    pos = pos ./ 10;  % px -> µm

    % ---- analysis window (concat P1..P4) ----
    start_time_sec = origStartFrame / origFps;
    start_frame    = max(1, round(start_time_sec * fps));
    framesPerPer   = round(secPerPeriod * fps);
    t0 = start_frame + (0:3)*framesPerPer;
    t1 = start_frame + (1:4)*framesPerPer;
    t1 = min(t1, size(dff,2));
    anal_frames = t0(1):t1(end);
    anal_frames = anal_frames(anal_frames>=1 & anal_frames<=size(dff,2));
    if numel(anal_frames) < 10
        warning('  %s: too few frames → skip.', fid);
        continue;
    end
    dffA = dff(:, anal_frames);

    % ---- fixed distance bins for every fish ----
    steps_used = F_binEdges;

    % ---- corr vs distance (split hemispheres by Y) ----
    [cd_out, ~, ~, pc_out] = corr_vs_distance_AO( ...
        dffA, pos, please_split, side2split, [], steps_used);

    % ---- Coerce outputs to [nHemis x nBins] with padding/truncation ----
    if iscell(cd_out)
        nH = numel(cd_out);
    else
        if isvector(cd_out)
            nH = 1;
        else
            nH = size(cd_out,1);
        end
    end

    cd_mat = nan(nH, nBins);     % correlations per hemisphere × bin
    pc_mat = zeros(nH, nBins);   % pair counts per hemisphere × bin

    for h = 1:nH
        % correlations
        if iscell(cd_out)
            v = cd_out{h};
        else
            if nH == 1
                v = cd_out(:).';
            else
                v = cd_out(h,:);
            end
        end
        v = v(:).';
        k = min(numel(v), nBins);
        cd_mat(h,1:k) = v(1:k);     % pad trailing with NaN

        % pair counts
        if isempty(pc_out)
            p = [];
        elseif iscell(pc_out)
            p = pc_out{h};
        else
            if nH == 1
                p = pc_out(:).';
            else
                p = pc_out(h,:);
            end
        end
        p = p(:).';
        k2 = min(numel(p), nBins);
        if k2 > 0
            pc_mat(h,1:k2) = p(1:k2);   % pad trailing with 0
        end
    end

    % aggregate hemispheres
    curve      = mean(cd_mat, 1, 'omitnan');   % 1×nBins
    pairs_this = sum(pc_mat, 1, 'omitnan');    % 1×nBins

    % ---- randomized control: position shuffles (keep activity, shuffle geometry) ----
    ctrl_cur = nan(nRand, nBins);
    for r = 1:nRand
        ridx = randperm(Ncells);
        [cdr, ~, ~, ~] = corr_vs_distance_AO( ...
            dffA, pos(ridx,:), please_split, side2split, [], steps_used);

        % coerce shuffled correlation to 1×nBins
        if iscell(cdr)
            if isempty(cdr)
                vr = [];
            else
                % if function returns both hemispheres, average here first
                tmp = nan(numel(cdr), nBins);
                for hh = 1:numel(cdr)
                    vv = cdr{hh}(:).';
                    kk = min(numel(vv), nBins);
                    tmp(hh,1:kk) = vv(1:kk);
                end
                vr = mean(tmp,1,'omitnan');
            end
        else
            if isvector(cdr)
                vr = cdr(:).';
            else
                vr = mean(cdr,1,'omitnan'); % average hemispheres if rows
            end
        end
        vr = vr(:).';
        kr = min(numel(vr), nBins);
        rowr = nan(1, nBins);
        if kr > 0
            rowr(1:kr) = vr(1:kr);
        end
        ctrl_cur(r,:) = rowr;
    end
    ctrl_mean = mean(ctrl_cur, 1, 'omitnan');
    ctrl_sem  = std(ctrl_cur, 0, 1, 'omitnan') ./ sqrt(sum(all(isfinite(ctrl_cur),2)));

    % ---- append one row for this fish ----
    F_curve_perFish      = [F_curve_perFish;      curve];
    F_ctrl_mean_perFish  = [F_ctrl_mean_perFish;  ctrl_mean];
    F_ctrl_sem_perFish   = [F_ctrl_sem_perFish;   ctrl_sem];
    F_pairs_perFish      = [F_pairs_perFish;      pairs_this];
    F_keptFish           = [F_keptFish;           string(fid)];
    F_nCells             = [F_nCells;             Ncells];
end

nFish = size(F_curve_perFish,1);
nBins = size(F_curve_perFish,2);

% For each fish, locate the last column where BOTH series are finite
lastFinite = zeros(nFish,1);
for k = 1:nFish
    ok = isfinite(F_curve_perFish(k,:)) & isfinite(F_ctrl_mean_perFish(k,:));
    if any(ok)
        lastFinite(k) = find(ok, 1, 'last');
    else
        lastFinite(k) = 0;  % this fish contributed nothing usable
    end
end

% Largest uniform prefix where *every* fish has data
L_common = min(lastFinite(lastFinite>0));
if isempty(L_common) || L_common < 2
    error('No common prefix across fish; nothing to pool/plot.');
end

% Optional: also restrict to ≤100 µm if desired
max_um = 100;
if any(F_binCenters <= max_um)
    L_um = find(F_binCenters <= max_um, 1, 'last');
    L_use = min(L_common, L_um);
else
    L_use = L_common;
end

% Trim all arrays to the uniform prefix
F_curve_perFish      = F_curve_perFish(:,      1:L_use);
F_ctrl_mean_perFish  = F_ctrl_mean_perFish(:,  1:L_use);
F_ctrl_sem_perFish   = F_ctrl_sem_perFish(:,   1:L_use);
F_pairs_perFish      = F_pairs_perFish(:,      1:L_use);
F_binCenters         = F_binCenters(                1:L_use);
F_binEdges           = F_binEdges(                 1:(L_use+1));  % keep edges aligned

% Quick diagnostics
fprintf('Common usable bins: %d (max distance ≈ %.1f µm)\n', L_use, F_binCenters(end));


%% ---------- PLOTS ----------
if isempty(F_curve_perFish)
    error('No fish completed corr–vs–distance.');
end

minPairsPerBin = 200;   % total across fish; tune
minFishPerBin  = 8;     % at least this many fish contributing

pairs_total  = nansum(F_pairs_perFish,1);                 % 1 x nbins
fish_counts  = sum(isfinite(F_curve_perFish),1);          % fish contributing per bin

w_curve = F_pairs_perFish ./ max(eps, nansum(F_pairs_perFish,1));   % rows sum to 1 per bin
w_ctrl  = w_curve;  % same weights (pair counts identical for shuffles)

pool_curve = nansum(w_curve .* F_curve_perFish,1);        % weighted mean
pool_ctrl  = nansum(w_ctrl  .* F_ctrl_mean_perFish,1);    % weighted mean

% Weighted SEM (delta method): sqrt( sum(w^2 * var) )
curve_var  = (F_curve_perFish - pool_curve).^2;
ctrl_var   = (F_ctrl_mean_perFish - pool_ctrl).^2;
pool_sem   = sqrt(nansum(w_curve.^2 .* curve_var,1));
pool_ctrl_sem = sqrt(nansum(w_ctrl.^2  .* ctrl_var,1));

% Mask poorly supported bins
mask = pairs_total >= minPairsPerBin & fish_counts >= minFishPerBin;
pool_curve(~mask)    = NaN;
pool_ctrl(~mask)     = NaN;
pool_sem(~mask)      = NaN;
pool_ctrl_sem(~mask) = NaN;


% Per-fish panels
figure('Color','w','Position',[60 60 1200 900]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for k = 1:size(F_curve_perFish,1)
    nexttile; hold on;
    Hc = shadedErrorBar(F_binCenters, F_ctrl_mean_perFish(k,:), F_ctrl_sem_perFish(k,:), ...
                        'lineProps', {'--k','LineWidth',1.8});
    Hb = plot(F_binCenters, F_curve_perFish(k,:), '-', 'LineWidth', 2, 'Color', blue);
    grid on; box off;
    xlabel('Distance (\mum)'); ylabel('Mean correlation');
    title(sprintf('Fish %s (N=%d cells)', F_keptFish(k), F_nCells(k)));
    legend([Hc.mainLine, Hb], {'Position-shuffle ctrl (mean \pm SEM)','All cells'}, ...
           'Location','northeast');
end
sgtitle('CCN2A — Per-fish correlation vs distance (all cells)');

% Pooled (equal weight per fish)
pool_curve    = mean(F_curve_perFish, 1, 'omitnan');
pool_ctrl     = mean(F_ctrl_mean_perFish, 1, 'omitnan');
pool_sem      = std(F_curve_perFish, 0, 1, 'omitnan') ./ sqrt(size(F_curve_perFish,1));
pool_ctrl_sem = std(F_ctrl_mean_perFish, 0, 1, 'omitnan') ./ sqrt(size(F_ctrl_mean_perFish,1));

figure('Color','w','Position',[120 120 900 500]); hold on;
Hc = shadedErrorBar(F_binCenters, pool_ctrl,  pool_ctrl_sem,  'lineProps', {'--k','LineWidth',2});
Hb = shadedErrorBar(F_binCenters, pool_curve, pool_sem,       'lineProps', {'-','Color',blue,'LineWidth',2});
grid on; box off;
xlabel('Distance (\mum)'); ylabel('Mean correlation');
title(sprintf('CCN2A — Pooled across fish (N=%d)', size(F_curve_perFish,1)));
legend([Hc.mainLine, Hb.mainLine], {'Position-shuffle ctrl (mean \pm SEM)','All cells (mean \pm SEM)'}, ...
       'Location','northeast');
xlim([0 max(F_binEdges)]);
xlabel('Distance (\mum)');
%%



%% ---------- Console summary ----------
auc_diff       = trapz(F_binCenters, F_curve_perFish - F_ctrl_mean_perFish, 2);
nf_bins        = 1:min(5, numel(F_binCenters));   % near field
nearfield_diff = mean(F_curve_perFish(:,nf_bins) - F_ctrl_mean_perFish(:,nf_bins), 2, 'omitnan');

fprintf('\n=== CCN2A Corr–vs–Distance (all cells) ===\n');
fprintf('Fish kept: %s\n', strjoin(cellstr(F_keptFish), ', '));
fprintf('Cells per fish: %s\n', mat2str(F_nCells.'));
fprintf('AUC(all - ctrl): median = %.4f  [IQR %.4f–%.4f]\n', ...
    median(auc_diff,'omitnan'), prctile(auc_diff,25), prctile(auc_diff,75));
fprintf('Near-field Δcorr (bins 1–%d): median = %.4f  [IQR %.4f–%.4f]\n', ...
    nf_bins(end), median(nearfield_diff,'omitnan'), prctile(nearfield_diff,25), prctile(nearfield_diff,75));
try
    [pAUC,~,stAUC] = signrank(auc_diff, 0);
    [pNF,~,stNF]   = signrank(nearfield_diff, 0);
    fprintf('Signrank AUC vs 0:        p=%.3g (z=%.2f)\n', pAUC, stAUC.zval);
    fprintf('Signrank Near-field vs 0: p=%.3g (z=%.2f)\n', pNF, stNF.zval);
catch
    fprintf('Signrank unavailable (Statistics Toolbox not found).\n');
end

%% ===================== HELPERS =====================
function fps = try_get_fps_from_metadata(mfile)
    fps = NaN;
    try
        S = load(mfile);
        fps = sniff_fps_in_struct(S);
    catch
        fps = NaN;
    end
end

function fps = sniff_fps_in_struct(S)
    fps = NaN;
    try
        fn = fieldnames(S);
        for i = 1:numel(fn)
            f = fn{i};
            val = S.(f);
            if isstruct(val)
                fps = sniff_fps_in_struct(val);
                if ~isnan(fps), return; end
            elseif isnumeric(val) && isscalar(val)
                fl = lower(f);
                if contains(fl,'fps') || strcmp(fl,'fs') || contains(fl,'framerate') || ...
                   contains(fl,'frame_rate') || contains(fl,'volume_rate') || ...
                   contains(fl,'volumerate') || contains(fl,'volume_rate_hz')
                    fps = double(val); return;
                end
            end
        end
    catch
        fps = NaN;
    end
end

function fps = try_get_fps_from_ini(iniPath)
    fps = NaN;
    try
        txt = fileread(iniPath);
        tok = regexp(txt, 'volume\.rate\.\(in\.Hz\)\s*=\s*([0-9.]+)', 'tokens');
        if ~isempty(tok), fps = str2double(tok{1}{1}); return; end
        tok = regexp(txt, 'frame[_ ]?rate\s*=\s*([0-9.]+)', 'tokens');
        if ~isempty(tok), fps = str2double(tok{1}{1}); return; end
    catch
        fps = NaN;
    end
end
