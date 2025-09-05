%% ============================================================
%  FIGURE 4 — Corr vs Distance for 0.05–0.07 Hz band (per fish)
%  Uses shadedErrorBar for SEM ribbons
%  Requires: corr_vs_distance_AO.m, shadedErrorBar.m
%  ============================================================

%% ====================== CONFIG ======================
baseFolder     = '\\forskning.it.ntnu.no\ntnu\mh-kin\yaksi5\anna\Data\Processed 2P data\mGluR multimodal\cell_detect_data\';
folderPattern  = '^[A-Za-z0-9]{4}_Wt$';     % or 'Het'

origStartFrame = 200;  % canonical start in raw movie domain
origFps        = 2.5;  % canonical fps used for origStartFrame

% frequency analysis params
windowSize     = 256;                 % frames (Welch window)
maxPeakCount   = 15;                  %#ok<NASGU> % kept for context
minPeakHeight  = 2;                   %#ok<NASGU>
frequencyRange = [0.01 0.10];         % Hz (search band for dominant frequency)

% low-frequency-of-interest (selection band for Figure 4)
lfs  = 0.05;                          % Hz
lffs = 0.07;                          % Hz

% randomized control
nRand = 200;                          % shuffles per fish
rng(0);                               % reproducible shuffles

% KDE grid (not used here; kept for continuity)
kdeGridN = 50;                        %#ok<NASGU>

% output folder (optional saving)
outRoot = fullfile(baseFolder, '_perfish_stability_outputs');
if ~exist(outRoot, 'dir'), mkdir(outRoot); end

%% ====================== FIND FISH FOLDERS ======================
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    [~, nm] = fileparts(allFolders{i});
    if ~isempty(regexp(nm, folderPattern, 'once'))
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
    end
end
fprintf('Found %d matching folders.\n', numel(validFolders));

resultsFiles = {};
fishLabel    = {};
for i = 1:numel(validFolders)
    D = dir(fullfile(validFolders{i}, '*_Results_dff.mat'));
    if ~isempty(D)
        resultsFiles{end+1} = fullfile(validFolders{i}, D(1).name); %#ok<SAGROW>
        [~, fishLabel{end+1}] = fileparts(validFolders{i});         %#ok<SAGROW>
    else
        warning('No *_Results_dff.mat in %s', validFolders{i});
    end
end

nF = numel(resultsFiles);
fprintf('Will analyze %d fish.\n', nF);

%% ====================== GROUP STORAGE (light) ======================
group = struct( ...
    'id',        strings(nF,1), ...
    'fps',       nan(nF,1), ...
    'mse12',     nan(nF,1), ...  % (unused here; kept for continuity)
    'mse13',     nan(nF,1), ...
    'mse14',     nan(nF,1), ...
    'ssim12',    nan(nF,1), ...
    'ssim13',    nan(nF,1), ...
    'ssim14',    nan(nF,1) ...
    );

%% ====================== OUTPUT ACCUMULATORS (Figure 4) ======================
F4_keptFish            = strings(0,1);
F4_nBandCells          = [];

F4_band_perFish        = [];    % rows = fish, cols = distance bins
F4_ctrl_mean_perFish   = [];    % same shape (mean across shuffles)
F4_ctrl_sem_perFish    = [];    % same shape (SEM across shuffles)
F4_binCenters          = [];    % set from first fish

%% ====================== PER-FISH ANALYSIS ======================
for f = 1:nF
    fprintf('\n[%d/%d] %s\n', f, nF, fishLabel{f});
    S = load(resultsFiles{f});

    % ---- DFF extraction ([cells x frames]) ----
    if isfield(S,'results') && isfield(S.results,'DV_DFFmovwindow')
        dff = double(S.results.DV_DFFmovwindow);
    else
        error('Expected S.results.DV_DFFmovwindow in %s', resultsFiles{f});
    end
    if size(dff,1) > size(dff,2)
        dff = dff.'; % ensure [cells x frames]
    end

    % ---- fps (fixed) ----
    fps = 2.5;   % your fixed fps
    group.id(f)  = string(fishLabel{f});
    group.fps(f) = fps;

    % ---- compute canonical analysis start in *this* fps ----
    start_time_sec = origStartFrame / origFps;         % seconds
    start_frame    = max(1, round(start_time_sec * fps));

    % ---- 4 fixed periods (each 2 min = 120 s) ----
    secPerPeriod = 120;
    framesPerPer = round(secPerPeriod * fps);

    P = struct( ...
        't0', start_frame + [0, 1, 2, 3]*framesPerPer, ...
        't1', start_frame + [1, 2, 3, 4]*framesPerPer  ...
    );
    % clamp to available frames
    P.t0 = min(P.t0, size(dff,2));
    P.t1 = min(P.t1, size(dff,2));

    % ---- positions (best-effort) ----
    positions = [];
    if isfield(S,'results')
        cand = {'positions_nodoubles','pos','xyz','XYZ','coords','centroids'};
        for c = 1:numel(cand)
            if isfield(S.results, cand{c})
                positions = S.results.(cand{c});
                break;
            end
        end
    end
    if isempty(positions)
        % try top-level fallbacks
        cand2 = {'positions_nodoubles','pos','xyz','coords'};
        for c = 1:numel(cand2)
            if isfield(S, cand2{c})
                positions = S.(cand2{c});
                break;
            end
        end
    end
    if ~isempty(positions) && size(positions,1) ~= size(dff,1)
        % try transpose
        if size(positions,2) == size(dff,1)
            positions = positions.';
        else
            warning('positions size mismatch; corr-vs-distance will be skipped for %s', fishLabel{f});
            positions = [];
        end
    end

    if isempty(positions)
        fprintf('  -> No positions found. Skipping fish %s.\n', fishLabel{f});
        continue
    end

    % ======== choose analysis window (concatenate P1–P4) ========
    anal_frames = P.t0(1):P.t1(4);
    anal_frames = anal_frames(anal_frames >= 1 & anal_frames <= size(dff,2));
    dffA = dff(:, anal_frames);

    % ======== dominant frequency per neuron (Hz) on analysis window ========
    T            = size(dffA,2);
    wlen         = min(windowSize, T);
    noverlap     = floor(wlen/2);
    nfft         = 2^nextpow2(max(wlen, 8*wlen));
    fmin         = frequencyRange(1);
    fmax         = frequencyRange(2);

    domF = nan(size(dffA,1),1);
    % precompute frequency grid (from first cell)
    [~, Fgrid] = pwelch(dffA(1,:)-mean(dffA(1,:)), hanning(wlen), noverlap, nfft, fps);
    bandMaskGrid = Fgrid >= fmin & Fgrid <= fmax;

    for iCell = 1:size(dffA,1)
        x = dffA(iCell,:) - mean(dffA(iCell,:));
        [Pxx, F] = pwelch(x, hanning(wlen), noverlap, nfft, fps);
        if numel(F) ~= numel(Fgrid) || any(F ~= Fgrid)
            bandMask = F >= fmin & F <= fmax;
            if any(bandMask)
                [~, imax] = max(Pxx(bandMask));
                fband = F(bandMask);
                domF(iCell) = fband(imax);
            end
        else
            if any(bandMaskGrid)
                [~, imax] = max(Pxx(bandMaskGrid));
                fband = Fgrid(bandMaskGrid);
                domF(iCell) = fband(imax);
            end
        end
    end

    % ======== select band cells (0.05–0.07 Hz) ========
    inBand  = isfinite(domF) & domF >= lfs & domF <= lffs;
    bandIdx = find(inBand);
    nBand   = numel(bandIdx);

    if nBand < 4
        fprintf('  -> Too few band cells (n=%d). Skipping fish %s.\n', nBand, fishLabel{f});
        continue
    end

    % ======== distance bins based on band cells ========
    posB    = positions(bandIdx, 1:2);
    pd      = squareform(pdist(posB));
    maxDist = max(pd(:), [], 'omitnan');
    if ~isfinite(maxDist) || maxDist < 50, maxDist = 100; end
    steps_used = floor(linspace(0, maxDist, 51));         % 50 bins
    binCenters = steps_used(2:end) - diff(steps_used)/2;

    % ======== corr vs distance (your function), split hemispheres on Y ========
    please_split = 1;
    side2split   = 2;    % split by Y

    [cd_band, ~, ~, ~] = corr_vs_distance_AO(dffA(bandIdx,:), positions(bandIdx,:), ...
                                             please_split, side2split, [], steps_used);
    band_curve = mean(cd_band, 1, 'omitnan');   % average hemispheres

    % ======== randomized control (same n, many shuffles; same steps) ========
    Ncells   = size(dffA,1);
    ctrl_cur = nan(nRand, numel(band_curve));
    for r = 1:nRand
        ridx = randperm(Ncells, nBand);
        [cd_r, ~, ~, ~] = corr_vs_distance_AO(dffA(ridx,:), positions(bandIdx,:), ...
                                              please_split, side2split, [], steps_used);
        ctrl_cur(r,:) = mean(cd_r, 1, 'omitnan');
    end

    % ======== stash per-fish results (align bins if needed) ========
    if isempty(F4_binCenters)
        F4_binCenters = binCenters;
    elseif numel(binCenters) ~= numel(F4_binCenters) || any(binCenters ~= F4_binCenters)
        band_curve = interp1(binCenters, band_curve, F4_binCenters, 'linear', 'extrap');
        ctrl_cur   = interp1(binCenters, ctrl_cur.', F4_binCenters, 'linear', 'extrap').';
    end

    F4_band_perFish(end+1,:)       = band_curve;
    F4_ctrl_mean_perFish(end+1,:)  = mean(ctrl_cur, 1, 'omitnan');
    F4_ctrl_sem_perFish(end+1,:)   = std(ctrl_cur, 0, 1, 'omitnan') ./ sqrt(sum(all(isfinite(ctrl_cur),2)));
    F4_keptFish(end+1,1)           = group.id(f);
    F4_nBandCells(end+1,1)         = nBand;
end

%% ====================== Figure 4 — PLOTS (with shadedErrorBar) ======================
if isempty(F4_band_perFish)
    error('Figure 4: no fish had >=4 band cells in 0.05–0.07 Hz.');
end

% ---- Per-fish panels ----
figure('Color','w','Position',[60 60 1200 900]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for k = 1:size(F4_band_perFish,1)
    nexttile; hold on;
    % control (mean ± SEM) in gray, dashed black line
    Hc = shadedErrorBar(F4_binCenters, F4_ctrl_mean_perFish(k,:), F4_ctrl_sem_perFish(k,:), ...
                        'lineProps', {'--k','LineWidth',1.8});
    % band (mean only; use SEM=0 to avoid extra ribbon OR reuse ctrl SEM if desired)
    % If you also want band SEM across hemispheres, you'd need hemisphere rows separately.
    Hb = plot(F4_binCenters, F4_band_perFish(k,:), 'r-', 'LineWidth', 2);

    grid on; box off;
    xlabel('Distance (pixels)'); ylabel('Mean correlation');
    title(sprintf('Fish %s (n_{band}=%d)', F4_keptFish(k), F4_nBandCells(k)));
    % Build legend using mainLine handle from shadedErrorBar
    legend([Hc.mainLine, Hb], {'Random ctrl (mean±SEM)','Band 0.05–0.07 Hz'}, 'Location','northeast');
end
sgtitle('Figure 4A — Per fish: correlation vs distance (band vs randomized)');

% ---- Pooled across fish (equal weight per fish) ----
pool_band      = mean(F4_band_perFish, 1, 'omitnan');
pool_ctrl      = mean(F4_ctrl_mean_perFish, 1, 'omitnan');
pool_band_sem  = std(F4_band_perFish, 0, 1, 'omitnan') ./ sqrt(size(F4_band_perFish,1));
pool_ctrl_sem  = std(F4_ctrl_mean_perFish, 0, 1, 'omitnan') ./ sqrt(size(F4_ctrl_mean_perFish,1));

figure('Color','w','Position',[120 120 900 500]); hold on;
% control pooled (mean ± SEM)
Hc = shadedErrorBar(F4_binCenters, pool_ctrl, pool_ctrl_sem, ...
                    'lineProps', {'--k','LineWidth',2});
% band pooled (mean ± SEM) in red
Hb = shadedErrorBar(F4_binCenters, pool_band, pool_band_sem, ...
                    'lineProps', {'-r','LineWidth',2});

grid on; box off;
xlabel('Distance (pixels)'); ylabel('Mean correlation');
title(sprintf('Figure 4B — Pooled across fish (N=%d)', size(F4_band_perFish,1)));
legend([Hc.mainLine, Hb.mainLine], {'Random ctrl (mean±SEM)', 'Band 0.05–0.07 Hz (mean±SEM)'}, ...
       'Location','northeast');

%% ====================== Quick console summary ======================
auc_diff       = trapz(F4_binCenters, F4_band_perFish - F4_ctrl_mean_perFish, 2);
nf_bins        = 1:min(5, numel(F4_binCenters));              % near-field: first bins
nearfield_diff = mean(F4_band_perFish(:,nf_bins) - F4_ctrl_mean_perFish(:,nf_bins), 2, 'omitnan');

fprintf('\n=== Figure 4 (0.05–0.07 Hz band) ===\n');
fprintf('Fish kept: %s\n', strjoin(cellstr(F4_keptFish), ', '));
fprintf('Band cells per fish: %s\n', mat2str(F4_nBandCells.'));
fprintf('AUC(band - ctrl): median = %.4f  [IQR %.4f–%.4f]\n', ...
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

%% ====================== (Optional) SAVE FIGURES ======================
% saveas(gcf, fullfile(outRoot, 'Figure4B_pooled_corr_vs_distance_shadedError.png'));
