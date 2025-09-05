%% ===================== CCN2A: Corr vs Distance (ALL CELLS) =====================
% Requires: corr_vs_distance_AO.m, shadedErrorBar.m

%% ---------- CONFIG ----------
baseFolder     = 'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';
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

%% ---------- ACCUMULATORS ----------
F_keptFish          = strings(0,1);
F_nCells            = [];
F_binCenters        = [];                 % set by first fish
F_curve_perFish     = [];                 % rows = fish, cols = distance bins
F_ctrl_mean_perFish = [];
F_ctrl_sem_perFish  = [];
% ----- fixed distance grid for all fish -----
binStep_um   = 5;                 % choose your resolution (e.g., 5 µm)
maxRange_um  = 300;               % target x-axis limit
F_binEdges   = 0:binStep_um:maxRange_um;
F_binCenters = F_binEdges(1:end-1) + diff(F_binEdges)/2;

%% ---------- PER-FISH LOOP ----------
for f = 1:numel(resultsFiles)
    fid = fishLabel{f};
    fps = fps_by_fish(f);
    fprintf('[%d/%d] %s (fps=%.3g)\n', f, numel(resultsFiles), fid, fps);

    S = load(resultsFiles{f});

    % ---- dFF [cells x frames] ----
    if isfield(S,'dff_new') 
        dff = double(S.dff_new);
    else
        warning('  %s: results.DV_DFFmovwindow not found → skip.', fid); continue;
    end
    if size(dff,1) > size(dff,2), dff = dff.'; end
    Ncells = size(dff,1);

    % ---- positions [cells x 2/3] ----
    pos = [];
    if isfield(S,'position')
        cand = {'positions_nodoubles','position','xyz','XYZ','coords','centroids'};
        for c = 1:numel(cand)
            if isfield(S, cand{c}) && ~isempty(S.(cand{c}))
                pos = double(S.(cand{c}));
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

    % ---- analysis window (P1..P4) ----
    start_time_sec = origStartFrame / origFps;
    start_frame    = max(1, round(start_time_sec * fps));
    framesPerPer   = round(secPerPeriod * fps);
    t0 = start_frame + (0:3)*framesPerPer;
    t1 = start_frame + (1:4)*framesPerPer;
    t1 = min(t1, size(dff,2));
    anal_frames = t0(1):t1(end);
    anal_frames = anal_frames(anal_frames>=1 & anal_frames<=size(dff,2));
    if numel(anal_frames) < 10
        warning('  %s: too few frames → skip.', fid); continue;
    end
    dffA = dff(:, anal_frames);

    % ---- distance bins from ALL cells ----
    % ---- distance bins: use the fixed grid for every fish ----
    steps_used = F_binEdges;

    % ---- corr vs distance (hemispheres split by Y) ----
    [cd_all, ~, ~, ~] = corr_vs_distance_AO(dffA, pos, ...
                                            please_split, side2split, [], steps_used);
    curve = mean(cd_all, 1, 'omitnan');                 % average hemispheres

    % ---- randomized control: position shuffles (keep activity, shuffle geometry) ----
    ctrl_cur = nan(nRand, numel(curve));
    for r = 1:nRand
        ridx = randperm(Ncells);
        [cd_r, ~, ~, ~] = corr_vs_distance_AO(dffA, pos(ridx,:), ...
                                              please_split, side2split, [], steps_used);
        ctrl_cur(r,:) = mean(cd_r, 1, 'omitnan');
    end



    % ---- stash ----
    F_curve_perFish(end+1,:)     = curve;                             %#ok<AGROW>
    F_ctrl_mean_perFish(end+1,:) = mean(ctrl_cur, 1, 'omitnan');      %#ok<AGROW>
    F_ctrl_sem_perFish(end+1,:)  = std(ctrl_cur, 0, 1, 'omitnan') ./ sqrt(sum(all(isfinite(ctrl_cur),2))); %#ok<AGROW>
    F_keptFish(end+1,1)          = string(fid);                       %#ok<AGROW>
    F_nCells(end+1,1)            = Ncells;                            %#ok<AGROW>
end

%% ---------- PLOTS ----------
if isempty(F_curve_perFish)
    error('No fish completed corr–vs–distance.');
end

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
