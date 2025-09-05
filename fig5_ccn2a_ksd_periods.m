%% ===================== CCN2A: Corr vs Distance + Period Spectra =====================
% Requires on path: corr_vs_distance_AO.m, calcsort_autocorr_freq_analysis_acs.m,
% shadedErrorBar.m (optional), try_get_fps_from_metadata.m, try_get_fps_from_ini.m

%% ---------- CONFIG ----------
baseFolder     = 'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';
origStartFrame = 200;      % where analysis window starts in the *raw movie*
origFps        = 2.5;      % fps used for origStartFrame
secPerPeriod   = 120;      % P1..P4: 4 x 120 s = 8 min

please_split   = 1;        % corr-vs-distance options
side2split     = 2;        % split by Y
nRand          = 200;      % position-shuffle control per fish
rng(0);

% low-frequency band for “stable” cells (match what you used elsewhere)
lfs       = 0.02;          % Hz
lffs      = 0.08;          % Hz
kdeGridN  = 60;            % KDE grid resolution

% output root
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
    fps_by_fish(i) = resolve_fps(resultsFiles{i}, metaFiles{i});
end


%% ---------- ACCUMULATORS (group level) ----------
F_keptFish          = strings(0,1);
F_nCells            = [];
binStep_um   = 5;
maxRange_um  = 300;
F_binEdges   = 0:binStep_um:maxRange_um;
F_binCenters = F_binEdges(1:end-1) + diff(F_binEdges)/2;
F_curve_perFish     = [];
F_ctrl_mean_perFish = [];
F_ctrl_sem_perFish  = [];

group = struct('mse12',nan(1,numel(resultsFiles)), 'mse13',nan(1,numel(resultsFiles)), ...
               'mse14',nan(1,numel(resultsFiles)), 'ssim12',nan(1,numel(resultsFiles)), ...
               'ssim13',nan(1,numel(resultsFiles)), 'ssim14',nan(1,numel(resultsFiles)));

%% ---------- PER-FISH LOOP ----------
for f = 1:numel(resultsFiles)
    fid = fishLabel{f};
    fps = fps_by_fish(f);
    fprintf('[%d/%d] %s (fps=%.3g)\n', f, numel(resultsFiles), fid, fps);

    S = load(resultsFiles{f});   % Suite2p outputs (we assume dff_new + position)

    % ---- dFF [cells x frames] ----
    if isfield(S,'dff_new')
        dff = double(S.dff_new);
    else
        warning('  %s: dff_new not found → skip.', fid); continue;
    end
    if size(dff,1) > size(dff,2), dff = dff.'; end  % cells x frames
    Ncells = size(dff,1);

    % ---- positions [cells x 2/3] ----
    pos = [];
    if isfield(S,'position')
        cand = {'positions_nodoubles','position','xyz','XYZ','coords','centroids'};
        for c = 1:numel(cand)
            if isfield(S, cand{c}) && ~isempty(S.(cand{c}))
                pos = double(S.(cand{c}));
                % your original script divides by 10 — keep that to be consistent
                pos = pos./10;  
                break;
            end
        end
    end
    if isempty(pos), warning('  %s: positions not found → skip.', fid); continue; end
    if size(pos,2) < 2, warning('  %s: positions need >=2 cols → skip.', fid); continue; end
    if size(pos,2) > 3, pos = pos(:,1:3); end
    if size(pos,1) ~= Ncells
        if size(pos,2) == Ncells, pos = pos.'; else, warning('  %s: positions size mismatch → skip.', fid); continue; end
    end

    %% ---------- analysis window (P1..P4) ----------
    start_time_sec = origStartFrame / origFps;
    start_frame    = max(1, round(start_time_sec * fps));
    framesPerPer   = round(secPerPeriod * fps);

    P.t0 = start_frame + (0:3)*framesPerPer;
    P.t1 = start_frame + (1:4)*framesPerPer;
    P.t1 = min(P.t1, size(dff,2));

    analysis_frames = P.t0(1):P.t1(end);
    analysis_frames = analysis_frames(analysis_frames>=1 & analysis_frames<=size(dff,2));
    if numel(analysis_frames) < 10
        warning('  %s: too few frames → skip.', fid); continue;
    end
    dffA = dff(:, analysis_frames);


    %% ---------- PERIOD-WISE SPECTRAL / KDE STABILITY ----------
    nPeriods = 4;
    domF        = cell(nPeriods,1);
    Pxxs        = cell(nPeriods,1);
    Fss         = cell(nPeriods,1);
    isdf_list   = cell(nPeriods,1);
    low_idx     = cell(nPeriods,1);

    num_lags = framesPerPer;        % autocorr lags cover one 2-min window

    for p = 1:nPeriods
        s0 = P.t0(p); s1 = P.t1(p);
        if s1 <= s0, error('Period %d window invalid for %s', p, fid); end
        dff_win = dff(:, s0:s1);

        % calcsort_* returns per-cell spectra and dominant frequencies
        [acs,~,~,~,~,dominant_frequencies,~,Pxx,Fss_p,~,~,isdf] = ...
            calcsort_autocorr_freq_analysis_acs( ...
                dff_win, num_lags, fps, 0.35, 128, 120);
        close
        domF{p}      = dominant_frequencies; %#ok<NASGU>
        Pxxs{p}      = Pxx;       % [nFreq x nCellsKept]
        Fss{p}       = Fss_p;     % frequency axis
        isdf_list{p} = isdf;

        lf = (dominant_frequencies >= lfs) & (dominant_frequencies <= lffs);
        low_idx{p} = lf;
    end

    % KDE density maps (2D) of low-freq cell positions
    positions = pos(:,1:2);
    density1      = NaN;
    density2_list = cell(3,1);
    density_diff  = cell(3,1);
    kde_grid      = struct();

    if ~isempty(positions)
        xlims = [min(positions(:,1)) max(positions(:,1))];
        ylims = [min(positions(:,2)) max(positions(:,2))];
        [xgrid, ygrid] = meshgrid(linspace(xlims(1), xlims(2), kdeGridN), ...
                                  linspace(ylims(1), ylims(2), kdeGridN));
        kde_grid.xgrid = xgrid; kde_grid.ygrid = ygrid;

        % P1
        idx1 = find(low_idx{1});
        if ~isempty(idx1)
            d1 = ksdensity([positions(idx1,1), positions(idx1,2)], [xgrid(:), ygrid(:)]);
            d1 = d1 ./ max(d1);
            density1 = reshape(d1, size(xgrid));
        else
            density1 = zeros(size(xgrid));
        end

        % P2–P4
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
        density1 = NaN; density2_list(:) = {NaN}; density_diff(:) = {NaN};
    end

    % MSE/SSIM between P1 and later periods, on common low-freq cells
    mse_vals  = nan(1,3);
    ssim_vals = nan(1,3);
    for k = 2:4
        lf1 = low_idx{1}; lfk = low_idx{k};
        common_idx = find(lf1 & lfk);
        if numel(common_idx) < 2
            fprintf('  %s Period %d: too few common low-freq cells (%d)\n', fid, k, numel(common_idx));
            mse_vals(k-1)  = NaN; ssim_vals(k-1) = NaN; continue;
        end
        A = Pxxs{1}(:, common_idx);
        B = Pxxs{k}(:, common_idx);
        minN = min(size(A,2), size(B,2));
        A = A(:,1:minN); B = B(:,1:minN);

        mse_vals(k-1) = mean((A(:) - B(:)).^2, 'omitnan');

        if exist('ssim','file') == 2
            try
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

    % store group metrics
    group.mse12(f)  = mse_vals(1);
    group.mse13(f)  = mse_vals(2);
    group.mse14(f)  = mse_vals(3);
    group.ssim12(f) = ssim_vals(1);
    group.ssim13(f) = ssim_vals(2);
    group.ssim14(f) = ssim_vals(3);

    % ---------- SAVE PER-FISH ----------
    perFishOut = struct();
    perFishOut.fish_id      = fid;
    perFishOut.fps          = fps;
    perFishOut.periods_sec  = [0 2; 2 4; 4 6; 6 8];
    perFishOut.Fss          = Fss;
    perFishOut.Pxxs         = Pxxs;
    perFishOut.low_idx      = low_idx;
    perFishOut.isdf_list    = isdf_list;
    perFishOut.density1     = density1;
    perFishOut.density2_all = density2_list;
    perFishOut.density_diff = density_diff;
    perFishOut.kde_grid     = kde_grid;
    perFishOut.mse_vals     = struct('P1vP2',mse_vals(1),'P1vP3',mse_vals(2),'P1vP4',mse_vals(3));
    perFishOut.ssim_vals    = struct('P1vP2',ssim_vals(1),'P1vP3',ssim_vals(2),'P1vP4',ssim_vals(3));
    perFishOut.freq_band    = [lfs lffs];

    save(fullfile(outRoot, sprintf('%s_periodStability.mat', fid)), '-struct', 'perFishOut', '-v7.3');
    fprintf('  %s  MSE P1vsP2/3/4 = [%.4g  %.4g  %.4g]\n', fid, mse_vals);
    fprintf('  %s SSIM P1vsP2/3/4 = [%.3f  %.3f  %.3f]\n', fid, ssim_vals);
end

%% ---------- SAVE GROUP SUMMARY ----------
save(fullfile(outRoot, 'group_periodStability_summary.mat'), ...
     'F_keptFish','F_nCells','F_binCenters','F_curve_perFish','F_ctrl_mean_perFish','F_ctrl_sem_perFish','group','-v7.3');

fprintf('\nDone. Outputs in: %s\n', outRoot);
