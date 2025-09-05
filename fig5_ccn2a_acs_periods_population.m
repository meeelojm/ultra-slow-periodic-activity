%% ===================== CONFIG (edit) =====================
% CCN2A sources (try all; first that exists wins)
baseRoots = { ...
  'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\high2low', ...
  'Z:\temp\Emilian temp\ccn2a\high2low', ...
  'Z:\Emilian temp\ccn2a\high2low'};

% if you want to restrict, list fish IDs; otherwise leave {} to auto-discover
fishList  = {};                          % e.g., {'f010','f011','f012'}
condList  = {'high','mid','low'};        % we search both cond and cond+"act"
condAliases = @(c){c, [c 'act']};

% where *_periodStability.mat live (we search recursively under this root)
outRoot   = 'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\high2low';

% canonical start (in RAW movie coordinates)
origStartFrame   = 200; 
origFps          = 2.5;                  % the fps of the *raw* movie that origStartFrame refers to
secPerPeriod     = 120;                  % 2 min per period window
lfs = 0.04;  lffs = 0.20;               % band for peak-power extraction

% target plotting fps for ACS x-axis & column alignment across fish
targetFpsForACS  = 2.5;                  % choose a nominal fps (2.5 keeps your original)
targetLags       = round(secPerPeriod * targetFpsForACS);

rng(13);                                 % for stable jitter in scatters

%% ============ 0) Build an index of Suite2p result folders ============
fprintf('Indexing CCN2A Suite2p folders…\n');
[idx, validFolders] = build_ccn2a_index(baseRoots, fishList, condList, condAliases);
fprintf('  Found %d Suite2p folders across %d fish.\n', numel(validFolders), numel(unique({idx.fish})));

%% ============ 1) Locate per-fish *_periodStability.mat files ============
PFF = dir(fullfile(outRoot, '**', '*_periodStability.mat'));
assert(~isempty(PFF), 'No *_periodStability.mat found under %s', outRoot);

% prep accumulators
ACS_stack           = {[],[],[],[]};    % {P1..P4} rows=cells (sorted by P1), cols=lags
pow_vals_by_period  = {[],[],[],[]};
fish_idx_by_period  = {[],[],[],[]};
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
    assert(isfield(PF,'isdf_list') && ~isempty(PF.isdf_list), '%s missing isdf_list', fid);
    assert(isfield(PF,'Pxxs') && numel(PF.Pxxs)>=4 && isfield(PF,'Fss'), '%s missing Pxxs/Fss', fid);

    % fps for THIS fish (from PF if stored, else try to sniff from results folder, else 2.5)
    if isfield(PF,'fps') && ~isempty(PF.fps) && isfinite(PF.fps)
        fps = PF.fps;
    else
        [resFile, fps_sniff] = findRes_ccn2a(fid, idx);
        fps = iff(isfinite(fps_sniff), fps_sniff, 2.5);
    end

    % ---- load raw DFF to re-compute ACS per period (acs wasn’t saved) ----
    if ~exist('resFile','var') || isempty(resFile)
        [resFile, ~] = findRes_ccn2a(fid, idx);
    end
    assert(~isempty(resFile) && isfile(resFile), 'Cannot find results.mat for %s', fid);

    S = load(resFile);
    if isfield(S,'results') && isfield(S.results,'DV_DFFmovwindow')
        dff = double(S.results.DV_DFFmovwindow);
    else
        error('Expected S.results.DV_DFFmovwindow in %s', resFile);
    end
    if size(dff,1) > size(dff,2), dff = dff.'; end
    nCells = size(dff,1);

    % global within-fish order
    baseOrder = PF.isdf_list{1}(:);
    if numel(baseOrder) ~= nCells
        warning('%s: isdf_list{1} length (%d) != nCells (%d). Using 1:nCells order.', ...
            fid, numel(baseOrder), nCells);
        baseOrder = (1:nCells).';
    end

    % period windows (frames in this fps)
    start_time_sec = origStartFrame / origFps;
    start_frame    = max(1, round(start_time_sec * fps));
    framesPerPer   = round(secPerPeriod * fps);
    num_lags       = framesPerPer;
    t0 = start_frame + (0:3)*framesPerPer;
    t1 = start_frame + (1:4)*framesPerPer;

    % ---- compute ACS per period, pick low-band cells in P1-order ----
    for p = 1:4
        s0 = t0(p); s1 = min(size(dff,2), t1(p));
        if s1 <= s0, warning('%s P%d: empty window; skipping.', fid, p); continue; end
        dff_win = dff(:, s0:s1);

        % returns acs and other outputs (we only need acs here)
        [acs,~,~,~,~,~,~,~,~,~,~,~] = calcsort_autocorr_freq_analysis_acs( ...
            dff_win, num_lags, fps, 0.35, 128, 120);

        % normalize orientation to [cells x lags]
        if size(acs,2) == nCells
            acs_c_by_lag = acs(:, baseOrder)';   % [lags x cells] -> [cells x lags], reorder
        else
            acs_c_by_lag = acs(baseOrder, :);    % already [cells x lags]
        end

        % --- align columns to common targetLags (so stacks concatenate cleanly)
        if size(acs_c_by_lag,2) ~= targetLags
            x_src = linspace(0, secPerPeriod, size(acs_c_by_lag,2));
            x_tar = linspace(0, secPerPeriod, targetLags);
            acs_c_by_lag = interp1(x_src, acs_c_by_lag.', x_tar, 'linear', 'extrap').';
        end

        % selection mask from PF (low-band cells for that period)
        mask = PF.low_idx{p};
        if islogical(mask) && numel(mask)==nCells
            sel = baseOrder(mask(baseOrder));  % keep P1 order but filter
        else
            warning('%s: low_idx{%d} invalid; using all cells.', fid, p);
            sel = baseOrder;
        end

        acs_sel = acs_c_by_lag(ismember(baseOrder, sel), :);   % rows in P1 order
        ACS_stack{p} = [ACS_stack{p}; acs_sel];                %#ok<AGROW>

        % ---- collect band-peak power for this period (per selected cell) ----
        Fp  = PF.Fss{p}(:);           % [nFreq x 1] or [1 x nFreq], make column
        Pxx = PF.Pxxs{p};             % expected [nFreq x nCells]

        % 1) Fix orientation → Pxx = [nFreq x nCells]
        if size(Pxx,2) ~= nCells && size(Pxx,1) == nCells, Pxx = Pxx.'; end

        % 2) Align F and Pxx row counts
        nFreqPxx = size(Pxx,1);
        nFreqFs  = numel(Fp);
        nf = min(nFreqPxx, nFreqFs);
        if nFreqPxx ~= nFreqFs
            warning('%s P%d: Fss(%d) vs Pxx rows(%d) mismatch → aligning to %d.', ...
                fid, p, nFreqFs, nFreqPxx, nf);
            Pxx = Pxx(1:nf, :);
            Fp  = Fp(1:nf);
        end

        % 3) Band mask (logical, length==size(Pxx,1))
        band = (Fp >= lfs) & (Fp <= lffs);

        % 4) Ensure 'sel' are valid col indices into Pxx
        sel = baseOrder(mask(baseOrder));
        sel = sel(sel >= 1 & sel <= size(Pxx,2));
        sel = unique(sel, 'stable');

        % 5) Compute band-limited peak power
        if any(band) && ~isempty(sel)
            pow_vec = max(Pxx(band, sel), [], 1);    % 1 x nSel
            pow_vals_by_period{p} = [pow_vals_by_period{p}, pow_vec]; %#ok<AGROW>
            fish_idx_by_period{p} = [fish_idx_by_period{p}, f * ones(1,numel(pow_vec))]; %#ok<AGROW>
        else
            warning('%s P%d: empty band or no selected cells — skipping power.', fid, p);
        end
    end
end

%% ================== PLOT: ACS heatmaps (aligned, pooled) ==================
figure('Color','w','Position',[50 50 1600 800]);
tl = tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
ptitles = {'ACS 0–2 min (P1)','ACS 2–4 min (P2)','ACS 4–6 min (P3)','ACS 6–8 min (P4)'};
for p = 1:4
    nexttile;
    imagesc(ACS_stack{p}); axis xy;
    title(ptitles{p});
    c = colorbar; caxis([0 0.5]);
    c.Label.String = 'Autocorrelation';
    % x in seconds using the common target length
    nt = size(ACS_stack{p},2);
    xtick = get(gca,'XTick');
    xtick(xtick<1 | xtick>nt) = [];
    xticklabels( (xtick-1) * (secPerPeriod/(nt-1)) );
    xtickangle(90);
    xlabel('Time (s)');
    ylabel('Neurons (sorted by P1)');
end
colormap(turbo);

%% ===== pooled (all neurons) band power per period =====
ptitles = {'P1 (0–2 m)','P2 (2–4 m)','P3 (4–6 m)','P4 (6–8 m)'};

% gather all values first to set a common y-axis
allY = []; Ns = zeros(1,4);
for p = 1:4
    y = pow_vals_by_period{p};
    if isempty(y), y = []; else, y = y(:); end
    y = y(isfinite(y)); Ns(p) = numel(y);
    allY  = [allY; y];
end
if isempty(allY), error('No power values collected.'); end
yl = [prctile(allY, 1) prctile(allY, 99)];

% Figure A: 4 panels, pooled neurons
figure('Color','w','Position',[80 80 1200 700]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for p = 1:4
    nexttile; hold on;
    y = pow_vals_by_period{p}; if isempty(y), axis off; title([ptitles{p} ' (no data)']); continue; end
    y = y(:); y = y(isfinite(y)); n = numel(y);
    boxchart(ones(n,1), y, 'BoxFaceAlpha',0.25, 'WhiskerLineStyle','-', 'MarkerStyle','none');
    xjit = 1 + 0.18*(rand(n,1) - 0.5);
    scatter(xjit, y, 6, 'filled', 'MarkerFaceAlpha', 0.25);
    mu = mean(y); se = std(y)/sqrt(max(1,n));
    errorbar(1.25, mu, se, 'k', 'LineWidth', 1.3, 'CapSize', 0);
    xlim([0.6 1.4]); xticks(1); xticklabels({'All neurons'});
    ylim(yl); grid on; ylabel('Peak Pxx in band');
    title(sprintf('%s — N=%d cells', ptitles{p}, n));
end

% Figure B: one axis, four periods side-by-side
figure('Color','w','Position',[100 100 900 450]); hold on;
xAll = []; yAll = [];
for p = 1:4
    y = pow_vals_by_period{p}; if isempty(y), continue; end
    y = y(:); y = y(isfinite(y));
    xAll = [xAll; p*ones(numel(y),1)];
    yAll = [yAll; y];
end
boxchart(xAll, yAll, 'BoxFaceAlpha',0.25, 'WhiskerLineStyle','-', 'MarkerStyle','none');
xjit = xAll + 0.18*(rand(size(xAll)) - 0.5);
scatter(xjit, yAll, 6, 'filled', 'MarkerFaceAlpha', 0.25);
for p = 1:4
    y = yAll(xAll==p); if isempty(y), continue; end
    mu = mean(y); se = std(y)/sqrt(numel(y));
    errorbar(p+0.25, mu, se, 'k', 'LineWidth', 1.3, 'CapSize', 0);
end
xlim([0.5 4.5]); xticks(1:4); xticklabels(ptitles);
ylim(yl); grid on; ylabel('Peak Pxx in band');
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
popOut.ACS_stack = ACS_stack;                 % {P1..P4}, each [cells x targetLags]
popOut.pow_vals_by_period = pow_vals_by_period;
popOut.fish_idx_by_period = fish_idx_by_period;
popOut.fish_names = fish_names;
popOut.band = [lfs lffs];
save(fullfile(outRoot, 'population_ACS_PXX_aggregated.mat'), '-struct', 'popOut', '-v7.3');
fprintf('Done. Saved aggregated ACS & PXX to %s\n', fullfile(outRoot, 'population_ACS_PXX_aggregated.mat'));

%% ================== HELPERS ==================

function [idx, validFolders] = build_ccn2a_index(baseRoots, fishList, condList, condAliases)
% Return an index with (fish, cond, suite2pPath, resultsPath, fps_guess)
    idx = struct('fish',{},'cond',{},'suite2p',{},'results',{},'fps',NaN);
    validFolders = {};
    % expand roots that exist
    roots = baseRoots(cellfun(@(p) exist(p,'dir')==7, baseRoots));
    if isempty(roots)
        error('None of the baseRoots exist on this machine.');
    end
    % discover fish by folder if not provided
    if isempty(fishList)
        fishSet = string.empty(0,1);
        for r = 1:numel(roots)
            D = dir(fullfile(roots{r}, 'f*'));
            fishSet = [fishSet; string({D([D.isdir]).name})'];
        end
        fishList = unique(cellstr(fishSet));
    end
    k = 0;
    for r = 1:numel(roots)
        for fi = 1:numel(fishList)
            fID = fishList{fi};
            for ci = 1:numel(condList)
                for alias = condAliases(condList{ci})
                    c = alias{1};
                    suiteDir = fullfile(roots{r}, fID, [c filesep 'suite2p']);
                    if exist(suiteDir,'dir') ~= 7
                        % also try without the explicit subfolder nesting
                        suiteDir = fullfile(roots{r}, fID, c, 'suite2p');
                    end
                    if exist(suiteDir,'dir') == 7
                        resMat = fullfile(suiteDir, 'results.mat');
                        if isfile(resMat)
                            k = k+1;
                            idx(k).fish    = fID;
                            idx(k).cond    = c;
                            idx(k).suite2p = suiteDir;
                            idx(k).results = resMat;
                            idx(k).fps     = try_get_fps_from_ini_or_default(suiteDir,  NaN);
                            validFolders{end+1} = suiteDir; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
    validFolders = unique(validFolders);
end

function [resFile, fps] = findRes_ccn2a(fid, idx)
% Pick a reasonable results.mat for a fish ID
    resFile = '';
    fps = NaN;
    if isempty(idx), return; end
    hits = find(strcmpi({idx.fish}, fid));
    if isempty(hits)
        % sometimes fid includes extra tokens; try substring search
        hits = find(contains(lower({idx.fish}), lower(fid)));
    end
    if isempty(hits), return; end
    % prefer 'low'/'lowact' by default; else first hit
    lowish = hits(contains(lower({idx(hits).cond}), 'low'));
    pick = iff(~isempty(lowish), lowish(1), hits(1));
    resFile = idx(pick).results;
    fps     = idx(pick).fps;
end

function v = iff(cond, a, b), if cond, v=a; else, v=b; end

function fps = try_get_fps_from_ini_or_default(folder, defaultVal)
    fps = defaultVal;
    try
        ii = dir(fullfile(folder, '*.ini'));
        if isempty(ii), ii = dir(fullfile(fileparts(folder), '*.ini')); end
        if ~isempty(ii)
            txt = fileread(fullfile(ii(1).folder, ii(1).name));
            tok = regexp(txt, 'volume\.rate\.\(in\.Hz\)\s*=\s*([0-9.]+)', 'tokens');
            if ~isempty(tok), fps = str2double(tok{1}{1}); return; end
            tok = regexp(txt, 'frame[_ ]?rate\s*=\s*([0-9.]+)', 'tokens');
            if ~isempty(tok), fps = str2double(tok{1}{1}); return; end
        end
    catch
        % keep default
    end
end
