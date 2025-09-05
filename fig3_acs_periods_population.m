%% ================== POPULATION ACS + POWER AGGREGATION ==================
baseFolder = '\\forskning.it.ntnu.no\ntnu\mh-kin\yaksi5\anna\Data\Processed 2P data\mGluR multimodal\cell_detect_data\';
outRoot    = fullfile(baseFolder, '_perfish_stability_outputs');
folderPattern = '^[A-Za-z0-9]{4}_Wt$';  % keep whatever you use
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
for ii = 1:numel(allFolders)
    if isempty(allFolders{ii}), continue; end
    [~, nm] = fileparts(allFolders{ii});
    if ~isempty(regexp(nm, folderPattern, 'once'))
        validFolders{end+1} = allFolders{ii}; %#ok<SAGROW>
    end
end

% canonical start (in raw movie coordinates)
origStartFrame = 200; origFps = 2.5;
secPerPeriod   = 120;              % 2 min
lfs = 0.04; lffs = 0.20;           % low-frequency band

% locate per-fish outputs
PFF = dir(fullfile(outRoot, '*_periodStability.mat'));
assert(~isempty(PFF), 'No per-fish periodStability files in %s', outRoot);

% prep accumulators
ACS_stack = {[],[],[],[]};          % {P1,P2,P3,P4} rows=cells (sorted by P1), cols=lags
pow_vals_by_period  = {[],[],[],[]};% y-values per period (peak power in [lfs lffs])
fish_idx_by_period  = {[],[],[],[]};% x-index (fish) for scatter/boxchart
fish_names          = strings(0,1);



fprintf('Aggregating ACS + Pxx across %d fish…\n', numel(PFF));
for f = 1:numel(PFF)
    % ---- load per-fish summary ----
    PF = load(fullfile(PFF(f).folder, PFF(f).name));
    [~, fid] = fileparts(PFF(f).name);
    fid = erase(fid, '_periodStability');
    fish_names(end+1) = string(fid);

    % sanity
    assert(isfield(PF,'low_idx') && numel(PF.low_idx)>=4, '%s missing low_idx{1..4}', fid);
    assert(isfield(PF,'isdf_list') && numel(PF.isdf_list)>=1, '%s missing isdf_list', fid);
    assert(isfield(PF,'Pxxs') && numel(PF.Pxxs)>=4 && isfield(PF,'Fss'), '%s missing Pxxs/Fss', fid);
    if isfield(PF,'fps') && ~isempty(PF.fps) && ~isnan(PF.fps), fps = PF.fps; else, fps = 2.5; end

    % ---- load raw DFF to re-compute ACS per period (acs wasn’t saved) ----
    resFile = findRes(fid, validFolders);
    assert(~isempty(resFile), 'Cannot find original *_Results_dff.mat for %s', fid);
    S = load(resFile);
    if isfield(S,'results') && isfield(S.results,'DV_DFFmovwindow')
        dff = double(S.results.DV_DFFmovwindow);
    else
        error('Expected S.results.DV_DFFmovwindow in %s', resFile);
    end
    if size(dff,1) > size(dff,2), dff = dff.'; end
    nCells = size(dff,1);
    baseOrder = PF.isdf_list{1}(:);                   % global within-fish order
    if numel(baseOrder) ~= nCells
        % fallback: if isdf_list was stored as a subset, create identity
        warning('%s: isdf_list{1} length != nCells. Using 1:nCells order.', fid);
        baseOrder = (1:nCells).';
    end

    % period windows (frames in this fps)
    start_time_sec = origStartFrame / origFps;
    start_frame    = max(1, round(start_time_sec * fps));
    framesPerPer   = round(secPerPeriod * fps);
    num_lags       = framesPerPer;    % 2-min lags
    t0 = start_frame + (0:3)*framesPerPer;
    t1 = start_frame + (1:4)*framesPerPer;

    % ---- compute ACS per period, pick low-band cells in P1-order ----
    for p = 1:4
        s0 = t0(p); s1 = min(size(dff,2), t1(p));
        dff_win = dff(:, s0:s1);

        % returns acs and other outputs (we only need acs here)
        [acs,~,~,~,~,~,~,~,~,~,~,~] = calcsort_autocorr_freq_analysis_acs( ...
            dff_win, num_lags, fps, 0.35, 128, 120);

        % normalize orientation to [cells x lags]
        if size(acs,2) == nCells
            % acs is [lags x cells] -> transpose
            acs_c_by_lag = acs(:, baseOrder)';  % reorder columns by baseOrder first
        else
            % acs is [cells x lags]
            acs_c_by_lag = acs(baseOrder, :);
        end

        mask = PF.low_idx{p};
        if islogical(mask) && numel(mask)==nCells
            sel = baseOrder(mask(baseOrder));        % keep P1 order but filter
        else
            warning('%s: low_idx{%d} invalid; using all cells.', fid, p);
            sel = baseOrder;
        end

        acs_sel = acs_c_by_lag(ismember(baseOrder, sel), :);   % rows in P1 order
        ACS_stack{p} = [ACS_stack{p}; acs_sel];                %#ok<AGROW>

        % ---- collect power for this period (band-peak per selected cell) ----
        Fp  = PF.Fss{p}(:);           % column vector of freqs
        Pxx = PF.Pxxs{p};             % expected [nFreq x nCells]
        
        % 1) Fix orientation if needed → make Pxx = [nFreq x nCells]
        if size(Pxx,2) ~= nCells && size(Pxx,1) == nCells
            Pxx = Pxx.';  % transpose
        end
        
        % 2) Align Fp and Pxx row counts
        nFreqPxx = size(Pxx,1);
        nFreqFs  = numel(Fp);
        nf = min(nFreqPxx, nFreqFs);
        if nFreqPxx ~= nFreqFs
            warning('%s P%d: Fss(%d) vs Pxx rows(%d) mismatch → aligning to %d.', ...
                fid, p, nFreqFs, nFreqPxx, nf);
            Pxx = Pxx(1:nf, :);
            Fp  = Fp(1:nf);
        end
        
        % 3) Band mask (logical vector length == size(Pxx,1))
        band = (Fp >= lfs) & (Fp <= lffs);
        
        % 4) Ensure 'sel' are valid column indices into Pxx
        sel = baseOrder(mask(baseOrder));                     % as before
        sel = sel(sel >= 1 & sel <= size(Pxx,2));             % clamp to cols
        sel = unique(sel, 'stable');                          % avoid dupes
        
        % 5) Compute band-limited peak power for selected cells
        if any(band) && ~isempty(sel)
            pow_vec = max(Pxx(band, sel), [], 1);             % 1 x nSel
            pow_vals_by_period{p} = [pow_vals_by_period{p}, pow_vec]; %#ok<AGROW>
            fish_idx_by_period{p} = [fish_idx_by_period{p}, f * ones(1,numel(pow_vec))]; %#ok<AGROW>
        else
            warning('%s P%d: empty band or no selected cells — skipping power.', fid, p);
        end

    end
end

%% ================== PLOT: ACS heatmaps (same per-fish ordering) ==================
% compute common color scale across all periods (robust to outliers)
% allvals = cell2mat(cellfun(@(A) A(:), ACS_stack, 'UniformOutput', false));
% if isempty(allvals), error('ACS stacks are empty.'); end
% v = prctile(abs(allvals), 99);
% cl = [-v v];

figure('Color','w','Position',[50 50 1600 800]);
tl = tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
titles = {'ACS 0–2 min (P1)','ACS 2–4 min (P2)','ACS 4–6 min (P3)','ACS 6–8 min (P4)'};
for p = 1:4
    nexttile;
    imagesc(ACS_stack{p}); axis xy;  % rows=cells (sorted), cols=lags
    title(titles{p});
    colorbar; caxis([0 0.5]);
    title('Autocorrelations sorted by dominant frequencies');
    xtick = get(gca, 'XTick');  % Get current x-tick values in frames
    xticks_in_seconds = xtick / fps;
    xtickangle(90);
    set(gca, 'XTickLabel', xticks_in_seconds);
    xlabel('Time (seconds)');
    ylabel('Neurons');
    c.Label.String = 'Autocorrelation value';  % Add label to the colorbar
end
colormap(turbo);

% ===== pooled (all neurons) band power per period =====
ptitles = {'P1 (0–2 m)','P2 (2–4 m)','P3 (4–6 m)','P4 (6–8 m)'};

% gather all values first to set a common y-axis
allY = [];
Ns   = zeros(1,4);
for p = 1:4
    y = pow_vals_by_period{p};
    if isempty(y), y = []; else, y = y(:); end
    y = y(isfinite(y));
    Ns(p) = numel(y);
    allY  = [allY; y];
end
if isempty(allY), error('No power values collected.'); end
yl = [prctile(allY, 1) prctile(allY, 99)];   % robust y-limits

% --- Figure A: four panels (one per period), pooled neurons ---
figure('Color','w','Position',[80 80 1200 700]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for p = 1:4
    nexttile; hold on;
    y = pow_vals_by_period{p};
    if isempty(y), axis off; title([ptitles{p} ' (no data)']); continue; end
    y = y(:); y = y(isfinite(y));
    n = numel(y);

    % box at x=1 with all neurons
    boxchart(ones(n,1), y, 'BoxFaceAlpha',0.25, 'WhiskerLineStyle','-', 'MarkerStyle','none');

    % jittered scatter: one dot per neuron (all fish pooled)
    xjit = 1 + 0.18*(rand(n,1) - 0.5);
    scatter(xjit, y, 6, 'filled', 'MarkerFaceAlpha', 0.25);  % small & semi-transparent

    % mean ± SEM overlay
    mu = mean(y); se = std(y)/sqrt(max(1,n));
    errorbar(1.25, mu, se, 'k', 'LineWidth', 1.3, 'CapSize', 0);

    xlim([0.6 1.4]); xticks(1); xticklabels({'All neurons'});
    ylim(yl); grid on;
    ylabel('Peak Pxx in band');
    title(sprintf('%s — N=%d cells', ptitles{p}, n));
end

% --- Figure B: one axis, four periods side-by-side, pooled neurons ---
figure('Color','w','Position',[100 100 900 450]); hold on;
xAll = []; yAll = [];
for p = 1:4
    y = pow_vals_by_period{p};
    if isempty(y), continue; end
    y = y(:); y = y(isfinite(y));
    xAll = [xAll; p*ones(numel(y),1)];
    yAll = [yAll; y];
end
% boxes per period
boxchart(xAll, yAll, 'BoxFaceAlpha',0.25, 'WhiskerLineStyle','-', 'MarkerStyle','none');
% scatter all neurons per period
xjit = xAll + 0.18*(rand(size(xAll)) - 0.5);
scatter(xjit, yAll, 6, 'filled', 'MarkerFaceAlpha', 0.25);

% means ± SEM per period
for p = 1:4
    y = yAll(xAll==p);
    if isempty(y), continue; end
    mu = mean(y); se = std(y)/sqrt(numel(y));
    errorbar(p+0.25, mu, se, 'k', 'LineWidth', 1.3, 'CapSize', 0);
end

xlim([0.5 4.5]); xticks(1:4); xticklabels(ptitles);
ylim(yl); grid on;
ylabel('Peak Pxx in band');
title('Pooled per-neuron band power across periods');

% quick console summary
means = zeros(1,4); sems = zeros(1,4);
for p = 1:4
    y = pow_vals_by_period{p}; y = y(:); y = y(isfinite(y));
    means(p) = mean(y); sems(p) = std(y)/sqrt(max(1,numel(y)));
end
fprintf('Pooled per-neuron peak power (mean ± SEM):\n');
fprintf('  P1: %.3g ± %.3g (N=%d)\n', means(1),sems(1), Ns(1));
fprintf('  P2: %.3g ± %.3g (N=%d)\n', means(2),sems(2), Ns(2));
fprintf('  P3: %.3g ± %.3g (N=%d)\n', means(3),sems(3), Ns(3));
fprintf('  P4: %.3g ± %.3g (N=%d)\n', means(4),sems(4), Ns(4));

%% ================== SAVE AGGREGATED OUTPUT ==================
popOut = struct();
popOut.ACS_stack = ACS_stack;                 % {P1..P4}, each [cells x lags]
popOut.pow_vals_by_period = pow_vals_by_period;
popOut.fish_idx_by_period = fish_idx_by_period;
popOut.fish_names = fish_names;
popOut.band = [lfs lffs];
save(fullfile(outRoot, 'population_ACS_PXX_aggregated.mat'), '-struct', 'popOut', '-v7.3');

fprintf('Done. Saved aggregated ACS & PXX to %s\n', fullfile(outRoot, 'population_ACS_PXX_aggregated.mat'));
