%% ====================== CONFIG ======================
baseFolder     = '\\forskning.it.ntnu.no\ntnu\mh-kin\yaksi5\anna\Data\Processed 2P data\mGluR multimodal\cell_detect_data\';
folderPattern  = '^[A-Za-z0-9]{4}_Wt$';     % or 'Het'

origStartFrame = 200;   % canonical start defined in raw movie domain
origFps        = 2.5;

% frequency analysis params
windowSize     = 256;
maxPeakCount   = 15;
minPeakHeight  = 2;
frequencyRange = [0.01 0.10];

% low-frequency-of-interest for densities/stability
lfs  = 0.05;
lffs = 0.07;

% KDE grid (in pixel coords of 'positions')
kdeGridN = 50;

% output folders
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
    'mse12',     nan(nF,1), ...  % P1 vs P2
    'mse13',     nan(nF,1), ...  % P1 vs P3
    'mse14',     nan(nF,1), ...  % P1 vs P4
    'ssim12',    nan(nF,1), ...
    'ssim13',    nan(nF,1), ...
    'ssim14',    nan(nF,1) ...
    );

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

    % ---- fps ----
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
            warning('positions size mismatch; density maps will be NaN for %s', fishLabel{f});
            positions = [];
        end
    end

    % ---- container for period-wise outputs ----
    nPeriods = 4;   % P1 (0–2), P2 (2–4), P3 (4–6), P4 (6–8)
    domF        = cell(nPeriods,1);
    Pxxs        = cell(nPeriods,1);
    Fss         = cell(nPeriods,1);
    isdf_list   = cell(nPeriods,1);
    low_idx     = cell(nPeriods,1);  % logical low-freq mask per period

    % ---- run frequency analysis per period ----
    num_lags = framesPerPer;  % = 2 min window for autocorr lags

    for p = 1:nPeriods
        s0 = P.t0(p); s1 = P.t1(p);
        if s1 <= s0, error('Period %d window invalid for %s', p, fishLabel{f}); end
        dff_win = dff(:, s0:s1);

        [acs,~,~,~,~,dominant_frequencies,~,Pxx,Fss_p,~,~,isdf] = ...
            calcsort_autocorr_freq_analysis_acs( ...
                dff_win, num_lags, fps, 0.35, 128, 120);

        domF{p}      = dominant_frequencies; %#ok<NASGU>
        Pxxs{p}      = Pxx;
        Fss{p}       = Fss_p;
        isdf_list{p} = isdf;

        lf = (dominant_frequencies >= lfs) & (dominant_frequencies <= lffs);
        low_idx{p} = lf;
    end

    % ================= DENSITY MAPS =================
    % density1 = low-freq map for Period 1 (0–2)
    % density2_{k} = low-freq map for Period k (k=2,3,4)
    % density_diff_{k} = density2_{k} - density1
    density1        = NaN;
    density2_list   = cell(3,1);   % P2,P3,P4
    density_diff    = cell(3,1);
    kde_grid        = struct();

    if ~isempty(positions)
        % grid
        xlims = [min(positions(:,1)) max(positions(:,1))];
        ylims = [min(positions(:,2)) max(positions(:,2))];
        [xgrid, ygrid] = meshgrid(linspace(xlims(1), xlims(2), kdeGridN), ...
                                  linspace(ylims(1), ylims(2), kdeGridN));
        kde_grid.xgrid = xgrid; kde_grid.ygrid = ygrid;

        % period 1
        idx1 = find(low_idx{1});
        if ~isempty(idx1)
            d1 = ksdensity([positions(idx1,1), positions(idx1,2)], [xgrid(:), ygrid(:)]);
            d1 = d1 ./ max(d1);
            density1 = reshape(d1, size(xgrid));
        else
            density1 = zeros(size(xgrid));
        end

        % later periods
        for k = 2:4
            idxk = find(low_idx{k});
            if ~isempty(idxk)
                d2 = ksdensity([positions(idxk,1), positions(idxk,2)], [xgrid(:), ygrid(:)]);
                d2 = d2 ./ max(d2);
                d2 = reshape(d2, size(xgrid));
            else
                d2 = zeros(size(xgrid));
            end
            density2_list{k-1} = d2;
            density_diff{k-1}  = d2 - density1;
        end
    else
        % no positions available
        density1 = NaN; density2_list(:) = {NaN}; density_diff(:) = {NaN};
    end

    % ================= STABILITY METRICS =================
    % Compare Pxx of common low-freq neurons between P1 and each later period
    mse_vals  = nan(1,3);
    ssim_vals = nan(1,3);

    for k = 2:4
        % neuron sets that are low-freq in BOTH periods
        lf1 = low_idx{1};
        lfk = low_idx{k};
        common_idx = find(lf1 & lfk);

        if numel(common_idx) < 2
            % too few to compute meaningful SSIM; MSE will still work with 1 col but keep consistent
            fprintf('  Period %d: too few common low-freq neurons (%d)\n', k, numel(common_idx));
            mse_vals(k-1)  = NaN;
            ssim_vals(k-1) = NaN;
            continue;
        end

        A = Pxxs{1}(:, common_idx);
        B = Pxxs{k}(:, common_idx);

        % shape check
        minN = min(size(A,2), size(B,2));
        A = A(:,1:minN); B = B(:,1:minN);

        % MSE
        mse_vals(k-1) = mean((A(:) - B(:)).^2, 'omitnan');

        % SSIM (fallback to NaN if function missing)
        if exist('ssim','file') == 2
            try
                % Normalize to [0,1] for numeric stability
                AA = A - min(A(:)); AA = AA ./ max(eps, max(AA(:)));
                BB = B - min(B(:)); BB = BB ./ max(eps, max(BB(:)));
                ssim_vals(k-1) = ssim(AA, BB);
            catch
                ssim_vals(k-1) = NaN;
            end
        else
            ssim_vals(k-1) = NaN;
        end
    end

    % store in group summary
    group.mse12(f)  = mse_vals(1);
    group.mse13(f)  = mse_vals(2);
    group.mse14(f)  = mse_vals(3);
    group.ssim12(f) = ssim_vals(1);
    group.ssim13(f) = ssim_vals(2);
    group.ssim14(f) = ssim_vals(3);

    % ================= SAVE PER-FISH =================
    perFishOut = struct();
    perFishOut.fish_id      = fishLabel{f};
    perFishOut.fps          = fps;
    perFishOut.periods_sec  = [0 2; 2 4; 4 6; 6 8];  % relative to canonical start
    perFishOut.Fss          = Fss;                   % 1x4
    perFishOut.Pxxs         = Pxxs;                  % 1x4 (each [nFreq x nCommonCells])
    perFishOut.low_idx      = low_idx;               % 1x4 logical
    perFishOut.isdf_list    = isdf_list;             % 1x4
    perFishOut.density1     = density1;
    perFishOut.density2_all = density2_list;         % {P2,P3,P4}
    perFishOut.density_diff = density_diff;          % {P2-P1, P3-P1, P4-P1}
    perFishOut.kde_grid     = kde_grid;
    perFishOut.mse_vals     = struct('P1vP2',mse_vals(1), 'P1vP3',mse_vals(2), 'P1vP4',mse_vals(3));
    perFishOut.ssim_vals    = struct('P1vP2',ssim_vals(1),'P1vP3',ssim_vals(2),'P1vP4',ssim_vals(3));
    perFishOut.freq_band    = [lfs lffs];

    save(fullfile(outRoot, sprintf('%s_periodStability.mat', fishLabel{f})), '-struct', 'perFishOut', '-v7.3');

    % quick log
    fprintf('  MSE P1vsP2/3/4 = [%.4g  %.4g  %.4g]\n', mse_vals);
    fprintf('  SSIM P1vsP2/3/4= [%.3f  %.3f  %.3f]\n', ssim_vals);
end

%% ====================== SAVE GROUP SUMMARY ======================
save(fullfile(outRoot, 'group_periodStability_summary.mat'), 'group', '-v7.3');

fprintf('\nDone. Per-fish files in: %s\n', outRoot);
