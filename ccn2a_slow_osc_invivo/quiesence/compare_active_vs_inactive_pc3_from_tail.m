%% =========================================================================
% compare_active_vs_inactive_pc3_from_tail.m
%
% Tail-defined least-active vs most-active windows:
% compare neural trajectory metrics in PC1-3 space
%
% Metrics:
%   - mean radius
%   - mean speed
%   - mean angular velocity
%   - cumulative rotations
%
% One script = one analysis / one figure set
% =========================================================================

clear; clc; close all;

%% =========================
% 0) USER SETTINGS
% =========================

rootFolder   = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low';
tailFileName = 'tail_quick.mat';
dffFileName  = 'dffs_repact_respcells.mat';

% Tail-defined windows
winMin  = 20;   % minutes
stepMin = 5;    % minutes

% NaN handling
maxNanFracPerNeuron = 0.20;

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_tail_pc3_results');

set(0,'DefaultFigureColor',[1 1 1]);

if saveOutputs && ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% =========================
% 1) FIND MATCHING RECORDINGS
% =========================

fprintf('\nSearching for recordings under:\n%s\n\n', rootFolder);

tailFiles = dir(fullfile(rootFolder, '**', tailFileName));
dffFiles  = dir(fullfile(rootFolder, '**', dffFileName));

if isempty(tailFiles), error('No %s files found.', tailFileName); end
if isempty(dffFiles),  error('No %s files found.', dffFileName);  end

dffMap = containers.Map;
for i = 1:numel(dffFiles)
    dffMap(dffFiles(i).folder) = fullfile(dffFiles(i).folder, dffFiles(i).name);
end

recordings = struct('folder', {}, 'tailPath', {}, 'dffPath', {});
for i = 1:numel(tailFiles)
    thisFolder = tailFiles(i).folder;
    if isKey(dffMap, thisFolder)
        recordings(end+1).folder = thisFolder; %#ok<SAGROW>
        recordings(end).tailPath = fullfile(tailFiles(i).folder, tailFiles(i).name);
        recordings(end).dffPath  = dffMap(thisFolder);
    end
end

if isempty(recordings)
    error('No folders contain both %s and %s.', tailFileName, dffFileName);
end

fprintf('Found %d matched recordings.\n', numel(recordings));

%% =========================
% 2) PREALLOCATE RESULTS
% =========================

nR = numel(recordings);

R = struct();
R.recordingLabel       = cell(nR,1);
R.folder               = cell(nR,1);

R.nNeurons_total       = nan(nR,1);
R.nNeurons_used_low    = nan(nR,1);
R.nNeurons_used_high   = nan(nR,1);
R.nNeurons_sharedPCA   = nan(nR,1);

R.fps                  = nan(nR,1);
R.nFrames              = nan(nR,1);
R.durationMin          = nan(nR,1);

R.low_startFrame       = nan(nR,1);
R.low_endFrame         = nan(nR,1);
R.high_startFrame      = nan(nR,1);
R.high_endFrame        = nan(nR,1);

R.low_boutCount        = nan(nR,1);
R.high_boutCount       = nan(nR,1);
R.low_boutRate_perMin  = nan(nR,1);
R.high_boutRate_perMin = nan(nR,1);

R.low_pc_radius_mean   = nan(nR,1);
R.high_pc_radius_mean  = nan(nR,1);

R.low_pc_speed_mean    = nan(nR,1);
R.high_pc_speed_mean   = nan(nR,1);

R.low_pc_angvel_mean   = nan(nR,1);
R.high_pc_angvel_mean  = nan(nR,1);

R.low_pc_rotations     = nan(nR,1);
R.high_pc_rotations    = nan(nR,1);

R.pc1_explained        = nan(nR,1);
R.pc2_explained        = nan(nR,1);
R.pc3_explained        = nan(nR,1);

allLowPCmetrics        = cell(nR,1);
allHighPCmetrics       = cell(nR,1);

%% =========================
% 3) PROCESS EACH RECORDING
% =========================

for r = 1:nR

    fprintf('\n============================================================\n');
    fprintf('Recording %d / %d\n', r, nR);
    fprintf('%s\n', recordings(r).folder);

    R.folder{r} = recordings(r).folder;
    R.recordingLabel{r} = make_recording_label(recordings(r).folder);

    % ---------------------------------------------------------------------
    % 3.1 Load DFF and FPS
    % ---------------------------------------------------------------------
    S_dff = load(recordings(r).dffPath);
    [dff, fps] = extract_dff_and_fps(S_dff);

    if isempty(dff)
        warning('Could not find DFF matrix in %s. Skipping.', recordings(r).dffPath);
        continue;
    end
    if isempty(fps) || ~isscalar(fps) || isnan(fps) || fps <= 0
        warning('Could not find valid fps in %s. Skipping.', recordings(r).dffPath);
        continue;
    end

    dff = double(dff);
    if size(dff,1) > size(dff,2)
        dff = dff.';
    end

    nNeurons = size(dff,1);
    nFrames  = size(dff,2);

    R.nNeurons_total(r) = nNeurons;
    R.fps(r)            = fps;
    R.nFrames(r)        = nFrames;
    R.durationMin(r)    = nFrames / fps / 60;

    fprintf('  DFF: %d neurons x %d frames | fps = %.4f | duration = %.2f min\n', ...
        nNeurons, nFrames, fps, R.durationMin(r));

    % ---------------------------------------------------------------------
    % 3.2 Load tail and extract bout onsets
    % ---------------------------------------------------------------------
    S_tail = load(recordings(r).tailPath);
    [boutOnsetsFrames, ~] = extract_bout_onsets_from_tail(S_tail, nFrames, fps);

    if isempty(boutOnsetsFrames)
        warning('No bout onsets extracted in %s. Skipping.', recordings(r).tailPath);
        continue;
    end

    boutOnsetsFrames = unique(round(boutOnsetsFrames(:)));
    boutOnsetsFrames = boutOnsetsFrames(isfinite(boutOnsetsFrames));
    boutOnsetsFrames = boutOnsetsFrames(boutOnsetsFrames >= 1 & boutOnsetsFrames <= nFrames);

    fprintf('  Extracted %d bout onsets.\n', numel(boutOnsetsFrames));

    % ---------------------------------------------------------------------
    % 3.3 Rank windows by tail activity
    % ---------------------------------------------------------------------
    winFrames  = max(1, round(winMin  * 60 * fps));
    stepFrames = max(1, round(stepMin * 60 * fps));

    if nFrames < winFrames
        warning('Recording shorter than one analysis window. Skipping.');
        continue;
    end

    [winTable, lowIdx, highIdx] = rank_windows_by_bout_activity( ...
        boutOnsetsFrames, nFrames, winFrames, stepFrames, fps);

    if isempty(winTable)
        warning('Could not build sliding windows for %s. Skipping.', recordings(r).folder);
        continue;
    end

    low_startFrame  = winTable.startFrame(lowIdx);
    low_endFrame    = winTable.endFrame(lowIdx);
    high_startFrame = winTable.startFrame(highIdx);
    high_endFrame   = winTable.endFrame(highIdx);

    R.low_startFrame(r)  = low_startFrame;
    R.low_endFrame(r)    = low_endFrame;
    R.high_startFrame(r) = high_startFrame;
    R.high_endFrame(r)   = high_endFrame;

    R.low_boutCount(r)        = winTable.boutCount(lowIdx);
    R.high_boutCount(r)       = winTable.boutCount(highIdx);
    R.low_boutRate_perMin(r)  = winTable.boutRatePerMin(lowIdx);
    R.high_boutRate_perMin(r) = winTable.boutRatePerMin(highIdx);

    fprintf('  Least active window: frames %d-%d | bouts = %d | %.3f bouts/min\n', ...
        low_startFrame, low_endFrame, R.low_boutCount(r), R.low_boutRate_perMin(r));
    fprintf('  Most active window : frames %d-%d | bouts = %d | %.3f bouts/min\n', ...
        high_startFrame, high_endFrame, R.high_boutCount(r), R.high_boutRate_perMin(r));

    % ---------------------------------------------------------------------
    % 3.4 Select neurons for shared PCA basis
    % ---------------------------------------------------------------------
    dff_low  = dff(:, low_startFrame:low_endFrame);
    dff_high = dff(:, high_startFrame:high_endFrame);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;

    dff_low_use  = dff_low(keepLow, :);
    dff_high_use = dff_high(keepHigh, :);

    dff_low_use  = dff_low_use(~all(isnan(dff_low_use),2), :);
    dff_high_use = dff_high_use(~all(isnan(dff_high_use),2), :);

    R.nNeurons_used_low(r)  = size(dff_low_use,1);
    R.nNeurons_used_high(r) = size(dff_high_use,1);

    keepBoth = keepLow & keepHigh & ~all(isnan(dff),2);
    dff_pca  = dff(keepBoth, :);
    R.nNeurons_sharedPCA(r) = size(dff_pca,1);

    if size(dff_pca,1) < 3
        warning('Not enough shared neurons for PCA in %s. Skipping.', recordings(r).folder);
        continue;
    end

    % ---------------------------------------------------------------------
    % 3.5 PCA on whole recording, project low/high windows, compute metrics
    % ---------------------------------------------------------------------
    dff_pca_filled = fill_nans_rowwise(dff_pca);

    muN = mean(dff_pca_filled, 2, 'omitnan');
    sdN = std(dff_pca_filled, 0, 2, 'omitnan');
    sdN(sdN==0 | ~isfinite(sdN)) = 1;
    dff_pca_z = (dff_pca_filled - muN) ./ sdN;

    Xfull = dff_pca_z.';  % frames x neurons

    try
        [coeff, ~, ~, ~, explained, pcaMu] = pca(Xfull, 'NumComponents', 3);
    catch ME
        warning('PCA failed in %s\n%s', recordings(r).folder, ME.message);
        continue;
    end

    if size(coeff,2) < 3
        warning('PCA returned fewer than 3 PCs in %s.', recordings(r).folder);
        continue;
    end

    Xlow  = dff_pca_z(:, low_startFrame:low_endFrame).';
    Xhigh = dff_pca_z(:, high_startFrame:high_endFrame).';

    score_low  = (Xlow  - pcaMu) * coeff(:,1:3);
    score_high = (Xhigh - pcaMu) * coeff(:,1:3);

    pcLow  = compute_pc3_metrics(score_low, fps);
    pcHigh = compute_pc3_metrics(score_high, fps);

    R.low_pc_radius_mean(r)  = pcLow.radius_mean;
    R.high_pc_radius_mean(r) = pcHigh.radius_mean;

    R.low_pc_speed_mean(r)   = pcLow.speed_mean;
    R.high_pc_speed_mean(r)  = pcHigh.speed_mean;

    R.low_pc_angvel_mean(r)  = pcLow.angvel_mean;
    R.high_pc_angvel_mean(r) = pcHigh.angvel_mean;

    R.low_pc_rotations(r)    = pcLow.rotations;
    R.high_pc_rotations(r)   = pcHigh.rotations;

    R.pc1_explained(r) = explained(1);
    R.pc2_explained(r) = explained(2);
    R.pc3_explained(r) = explained(3);

    allLowPCmetrics{r}  = pcLow;
    allHighPCmetrics{r} = pcHigh;

    fprintf('  PCA explained variance PC1-3: %.2f / %.2f / %.2f %%\n', ...
        explained(1), explained(2), explained(3));
    fprintf('  PC low : radius=%.4f, speed=%.4f, angvel=%.4f, rotations=%.4f\n', ...
        pcLow.radius_mean, pcLow.speed_mean, pcLow.angvel_mean, pcLow.rotations);
    fprintf('  PC high: radius=%.4f, speed=%.4f, angvel=%.4f, rotations=%.4f\n', ...
        pcHigh.radius_mean, pcHigh.speed_mean, pcHigh.angvel_mean, pcHigh.rotations);
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used_low, R.nNeurons_used_high, R.nNeurons_sharedPCA, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.low_startFrame, R.low_endFrame, R.high_startFrame, R.high_endFrame, ...
    R.low_boutCount, R.high_boutCount, R.low_boutRate_perMin, R.high_boutRate_perMin, ...
    R.low_pc_radius_mean, R.high_pc_radius_mean, ...
    R.low_pc_speed_mean, R.high_pc_speed_mean, ...
    R.low_pc_angvel_mean, R.high_pc_angvel_mean, ...
    R.low_pc_rotations, R.high_pc_rotations, ...
    R.pc1_explained, R.pc2_explained, R.pc3_explained, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used_low','nNeurons_used_high','nNeurons_sharedPCA', ...
    'fps','nFrames','durationMin', ...
    'low_startFrame','low_endFrame','high_startFrame','high_endFrame', ...
    'low_boutCount','high_boutCount','low_boutRate_perMin','high_boutRate_perMin', ...
    'low_pc_radius_mean','high_pc_radius_mean', ...
    'low_pc_speed_mean','high_pc_speed_mean', ...
    'low_pc_angvel_mean','high_pc_angvel_mean', ...
    'low_pc_rotations','high_pc_rotations', ...
    'pc1_explained','pc2_explained','pc3_explained'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_pc3_active_vs_inactive.csv'));
    save(fullfile(outFolder, 'summary_pc3_active_vs_inactive.mat'), ...
        'R', 'T', 'allLowPCmetrics', 'allHighPCmetrics');
end

%% =========================
% 5) PAIRED STATISTICS
% =========================

statsSummary = struct();
statsSummary.pcRadius   = run_paired_signrank(R.low_pc_radius_mean, R.high_pc_radius_mean);
statsSummary.pcSpeed    = run_paired_signrank(R.low_pc_speed_mean, R.high_pc_speed_mean);
statsSummary.pcAngVel   = run_paired_signrank(R.low_pc_angvel_mean, R.high_pc_angvel_mean);
statsSummary.pcRotation = run_paired_signrank(R.low_pc_rotations, R.high_pc_rotations);

disp('================ PAIRED SIGNRANK TESTS ================')
fprintf('PC radius:            n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.pcRadius.n, statsSummary.pcRadius.p, statsSummary.pcRadius.deltaMedian);
fprintf('PC speed:             n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.pcSpeed.n, statsSummary.pcSpeed.p, statsSummary.pcSpeed.deltaMedian);
fprintf('PC angular velocity:  n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.pcAngVel.n, statsSummary.pcAngVel.p, statsSummary.pcAngVel.deltaMedian);
fprintf('PC rotations:         n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.pcRotation.n, statsSummary.pcRotation.p, statsSummary.pcRotation.deltaMedian);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_pc3.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURE
% =========================

validRec = isfinite(R.low_pc_radius_mean) & isfinite(R.high_pc_radius_mean) & ...
           isfinite(R.low_pc_speed_mean) & isfinite(R.high_pc_speed_mean) & ...
           isfinite(R.low_pc_angvel_mean) & isfinite(R.high_pc_angvel_mean) & ...
           isfinite(R.low_pc_rotations) & isfinite(R.high_pc_rotations);

if any(validRec)

    xPair = [1 2];
    idxC  = find(validRec(:))';

    figC = figure('Name','Active vs inactive PC3 trajectory comparison', ...
        'Position',[100 100 1400 700]);
    tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact');

    % Radius
    nexttile; hold on;
    yLow = R.low_pc_radius_mean(validRec);
    yHigh = R.high_pc_radius_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxC
        plot(xPair, [R.low_pc_radius_mean(i), R.high_pc_radius_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean radius in PC1-3');
    title(sprintf('Radius\nn=%d, p=%.3g %s', ...
        statsSummary.pcRadius.n, statsSummary.pcRadius.p, p_to_stars(statsSummary.pcRadius.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.pcRadius.p);

    % Speed
    nexttile; hold on;
    yLow = R.low_pc_speed_mean(validRec);
    yHigh = R.high_pc_speed_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxC
        plot(xPair, [R.low_pc_speed_mean(i), R.high_pc_speed_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean speed in PC1-3');
    title(sprintf('Speed\nn=%d, p=%.3g %s', ...
        statsSummary.pcSpeed.n, statsSummary.pcSpeed.p, p_to_stars(statsSummary.pcSpeed.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.pcSpeed.p);

    % Angular velocity
    nexttile; hold on;
    yLow = R.low_pc_angvel_mean(validRec);
    yHigh = R.high_pc_angvel_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxC
        plot(xPair, [R.low_pc_angvel_mean(i), R.high_pc_angvel_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean angular velocity');
    title(sprintf('Angular velocity\nn=%d, p=%.3g %s', ...
        statsSummary.pcAngVel.n, statsSummary.pcAngVel.p, p_to_stars(statsSummary.pcAngVel.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.pcAngVel.p);

    % Rotations
    nexttile; hold on;
    yLow = R.low_pc_rotations(validRec);
    yHigh = R.high_pc_rotations(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxC
        plot(xPair, [R.low_pc_rotations(i), R.high_pc_rotations(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Cumulative rotations');
    title(sprintf('Rotations\nn=%d, p=%.3g %s', ...
        statsSummary.pcRotation.n, statsSummary.pcRotation.p, p_to_stars(statsSummary.pcRotation.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.pcRotation.p);

    sgtitle('PC1-3 trajectory comparison: tail-defined inactive vs active windows');

    if saveOutputs
        saveas(figC, fullfile(outFolder, 'summary_pc3_active_vs_inactive.png'));
    end
else
    warning('No valid recordings available for PC trajectory summary figure.');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% =========================
% LOCAL FUNCTIONS
% =========================

function [dff, fps] = extract_dff_and_fps(S)
% Try to robustly extract dff matrix and fps from loaded dff mat file

    dff = [];
    fps = [];

    vars = fieldnames(S);

    % -------- find dff
    preferredDffNames = { ...
        'dff', 'dffs', 'dFF', 'DFF', ...
        'dff_respcells', 'dffs_respcells', ...
        'dff_repact_respcells', 'dffs_repact_respcells', ...
        'F7', 'F', 'allDff'};

    for i = 1:numel(preferredDffNames)
        nm = preferredDffNames{i};
        if isfield(S, nm) && isnumeric(S.(nm)) && ndims(S.(nm)) == 2
            dff = S.(nm);
            break;
        end
    end

    % fallback: choose biggest 2D numeric matrix
    if isempty(dff)
        bestSize = 0;
        for i = 1:numel(vars)
            v = S.(vars{i});
            if isnumeric(v) && ismatrix(v)
                sz = numel(v);
                if sz > bestSize
                    bestSize = sz;
                    dff = v;
                end
            end
        end
    end

    % -------- find fps
    preferredFpsNames = {'fps','Fs','fs','frameRate','framerate','samplingRate'};
    for i = 1:numel(preferredFpsNames)
        nm = preferredFpsNames{i};
        if isfield(S, nm) && isnumeric(S.(nm)) && isscalar(S.(nm))
            fps = double(S.(nm));
            break;
        end
    end

    % nested struct search for fps
    if isempty(fps)
        for i = 1:numel(vars)
            v = S.(vars{i});
            if isstruct(v)
                fns = fieldnames(v);
                for j = 1:numel(fns)
                    if any(strcmpi(fns{j}, preferredFpsNames))
                        val = v.(fns{j});
                        if isnumeric(val) && isscalar(val)
                            fps = double(val);
                            return;
                        end
                    end
                end
            end
        end
    end
end

function [boutOnsetsFrames, tailActivity] = extract_bout_onsets_from_tail(S_tail, nFrames, fps)
% Robustly extract bout onsets from tail_quick.mat
%
% Tries, in order:
%   1) direct bout onset variable in frames
%   2) direct bout onset variable in seconds -> converts to frames
%   3) direct binary bout train
%   4) thresholded tail activity peaks from continuous tail trace

    boutOnsetsFrames = [];
    tailActivity = [];

    vars = fieldnames(S_tail);

    % -------- try obvious direct frame variables
    candidateFrameNames = { ...
        'boutOnsets', 'bout_onsets', 'boutStarts', 'bout_starts', ...
        'boutOnsetFrames', 'boutStartFrames', 'bout_frames', ...
        'onsets', 'bout_train_onsets'};

    for i = 1:numel(candidateFrameNames)
        nm = candidateFrameNames{i};
        if isfield(S_tail, nm) && isnumeric(S_tail.(nm))
            x = S_tail.(nm);
            if isvector(x)
                x = x(:);
                if max(x) <= nFrames*1.2
                    boutOnsetsFrames = round(x);
                    return;
                end
            end
        end
    end

    % -------- try obvious direct second variables
    candidateSecNames = { ...
        'boutOnsetsSec', 'bout_onsets_sec', 'boutTimes', 'bout_times', ...
        'onsetTimes', 'onsets_sec'};

    for i = 1:numel(candidateSecNames)
        nm = candidateSecNames{i};
        if isfield(S_tail, nm) && isnumeric(S_tail.(nm))
            x = S_tail.(nm);
            if isvector(x)
                x = x(:);
                boutOnsetsFrames = round(x * fps) + 1;
                return;
            end
        end
    end

    % -------- try binary bout train
    candidateTrainNames = { ...
        'bout_train', 'boutTrain', 'bouts', 'bout_binary', ...
        'boutVec', 'boutVector', 'tail_bout_train'};

    for i = 1:numel(candidateTrainNames)
        nm = candidateTrainNames{i};
        if isfield(S_tail, nm) && isnumeric(S_tail.(nm))
            bt = S_tail.(nm);
            if isvector(bt)
                bt = bt(:);
                if numel(bt) ~= nFrames
                    bt = try_resample_vector(bt, nFrames);
                end
                bt = bt > 0;
                boutOnsetsFrames = find(diff([0; bt]) == 1);
                return;
            end
        end
    end

    % -------- try continuous tail activity variables
    candidateTraceNames = { ...
        'tail', 'tail_trace', 'tailTrace', 'tail_angle', 'tailAngle', ...
        'smoothedTail', 'tai', 'trace', 'tail_mov', 'tailMovement', ...
        'tail_power', 'tail_activity'};

    for i = 1:numel(candidateTraceNames)
        nm = candidateTraceNames{i};
        if isfield(S_tail, nm) && isnumeric(S_tail.(nm)) && isvector(S_tail.(nm))
            tailActivity = double(S_tail.(nm)(:));
            break;
        end
    end

    % fallback: choose numeric vector with length closest to nFrames
    if isempty(tailActivity)
        bestDist = inf;
        bestVec = [];
        for i = 1:numel(vars)
            v = S_tail.(vars{i});
            if isnumeric(v) && isvector(v) && numel(v) > 100
                d = abs(numel(v) - nFrames);
                if d < bestDist
                    bestDist = d;
                    bestVec = double(v(:));
                end
            end
        end
        tailActivity = bestVec;
    end

    if isempty(tailActivity)
        return;
    end

    % match length to imaging frames
    if numel(tailActivity) ~= nFrames
        tailActivity = try_resample_vector(tailActivity, nFrames);
    end

    % smooth lightly and detect bouts from envelope / peaks
    tailActivity = fillmissing(tailActivity, 'linear', 'EndValues', 'nearest');
    tailActivity = smoothdata(tailActivity, 'gaussian', 5);

    % robust threshold
    thr = median(tailActivity, 'omitnan') + 2 * mad(tailActivity, 1);
    active = tailActivity > thr;

    % remove tiny fragments
    minBoutDurFrames = max(2, round(0.10 * fps));  % 100 ms minimum
    active = remove_short_true_runs(active, minBoutDurFrames);

    boutOnsetsFrames = find(diff([0; active(:)]) == 1);
end

function [winTable, lowIdx, highIdx] = rank_windows_by_bout_activity(boutOnsetsFrames, nFrames, winFrames, stepFrames, fps)

    startFrames = 1:stepFrames:(nFrames - winFrames + 1);
    if isempty(startFrames)
        winTable = table();
        lowIdx = [];
        highIdx = [];
        return;
    end

    nW = numel(startFrames);
    endFrames    = startFrames + winFrames - 1;
    centerFrames = round((startFrames + endFrames)/2);

    boutCount = nan(nW,1);
    for w = 1:nW
        boutCount(w) = sum(boutOnsetsFrames >= startFrames(w) & boutOnsetsFrames <= endFrames(w));
    end

    boutRatePerMin = boutCount ./ (winFrames / fps / 60);

    winTable = table(startFrames(:), endFrames(:), centerFrames(:), boutCount(:), boutRatePerMin(:), ...
        'VariableNames', {'startFrame','endFrame','centerFrame','boutCount','boutRatePerMin'});

    % choose least active and most active windows
    % tie-breaker: earliest one among ties
    [~, lowIdx]  = min(winTable.boutCount);
    [~, highIdx] = max(winTable.boutCount);
end

function out = run_paired_signrank(xLow, xHigh)
    idx = isfinite(xLow) & isfinite(xHigh);
    xLow = xLow(idx);
    xHigh = xHigh(idx);

    out = struct('p',NaN,'h',NaN,'signedrank',NaN,'n',numel(xLow), ...
                 'deltaMedian',NaN);

    if numel(xLow) >= 2
        [p,h,stats] = signrank(xLow, xHigh);
        out.p = p;
        out.h = h;
        out.signedrank = stats.signedrank;
        out.deltaMedian = median(xHigh - xLow, 'omitnan');
    end
end

function s = p_to_stars(p)
    if ~isfinite(p)
        s = 'n.s.';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end

function add_sig_bar(ax, xPair, yVals, p)
    if isempty(yVals) || all(~isfinite(yVals))
        return;
    end

    yVals = yVals(isfinite(yVals));
    yMax = max(yVals);
    yMin = min(yVals);

    if yMax == yMin
        yRange = max(abs(yMax), 1);
    else
        yRange = yMax - yMin;
    end

    yBar = yMax + 0.08*yRange;
    yTxt = yMax + 0.12*yRange;

    axes(ax); hold(ax,'on');
    plot(ax, [xPair(1) xPair(1) xPair(2) xPair(2)], ...
        [yBar-0.015*yRange yBar yBar yBar-0.015*yRange], ...
        'k-', 'LineWidth', 1.2);

    text(mean(xPair), yTxt, p_to_stars(p), ...
        'Parent', ax, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold', ...
        'FontSize', 12);

    ylim(ax, [min(get(ax,'YLim')) yMax + 0.20*yRange]);
end

function label = make_recording_label(folderPath)
    parts = regexp(folderPath, '[\\/]', 'split');
    if numel(parts) >= 3
        label = strjoin(parts(end-2:end), '_');
    else
        label = folderPath;
    end
end

function x = remove_short_true_runs(x, minLen)
    x = logical(x(:));
    d = diff([false; x; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    lens   = ends - starts + 1;
    shortRuns = find(lens < minLen);
    for i = 1:numel(shortRuns)
        x(starts(shortRuns(i)):ends(shortRuns(i))) = false;
    end
end

function y = try_resample_vector(x, targetLen)
    x = x(:);
    if isempty(x) || targetLen < 2
        y = x;
        return;
    end
    xi = linspace(1, numel(x), targetLen);
    y = interp1(1:numel(x), x, xi, 'linear', 'extrap')';
end

function X = fill_nans_rowwise(X)
    for i = 1:size(X,1)
        xi = X(i,:);
        if any(isnan(xi))
            xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
            if any(isnan(xi))
                xi = fillmissing(xi, 'constant', 0);
            end
        end
        X(i,:) = xi;
    end
end

function out = compute_pc3_metrics(score3, fps)
% score3: frames x 3

    out = struct();

    if isempty(score3) || size(score3,1) < 3 || size(score3,2) < 3
        out.radius_mean = NaN;
        out.speed_mean  = NaN;
        out.angvel_mean = NaN;
        out.rotations   = NaN;
        out.radius_t    = [];
        out.speed_t     = [];
        out.angvel_t    = [];
        out.theta_step  = [];
        out.centroid    = [];
        return;
    end

    dt = 1 / fps;

    % Center trajectory in PC space
    ctr = mean(score3, 1, 'omitnan');
    R = score3 - ctr;   % frames x 3

    % Radius from centroid
    radius_t = sqrt(sum(R.^2, 2));

    % Velocity in PC space
    V = diff(score3, 1, 1) / dt;   % (T-1) x 3

    % Speed
    speed_t = sqrt(sum(V.^2, 2));

    % Angular velocity magnitude around centroid:
    % omega = |r x v| / |r|^2
    Rmid = R(1:end-1,:);   % align with V
    crossRV = cross(Rmid, V, 2);
    num = sqrt(sum(crossRV.^2, 2));
    den = sum(Rmid.^2, 2);
    den(den < eps) = NaN;
    angvel_t = num ./ den;

    % Rotations: cumulative angular displacement / (2*pi)
    r1 = R(1:end-1,:);
    r2 = R(2:end,:);
    n1 = sqrt(sum(r1.^2,2));
    n2 = sqrt(sum(r2.^2,2));
    denom = n1 .* n2;
    denom(denom < eps) = NaN;

    cosang = sum(r1 .* r2, 2) ./ denom;
    cosang = max(-1, min(1, cosang));
    theta_step = acos(cosang);   % radians, unsigned
    rotations = nansum(theta_step) / (2*pi);

    out.radius_mean = mean(radius_t, 'omitnan');
    out.speed_mean  = mean(speed_t, 'omitnan');
    out.angvel_mean = mean(angvel_t, 'omitnan');
    out.rotations   = rotations;

    out.radius_t    = radius_t;
    out.speed_t     = speed_t;
    out.angvel_t    = angvel_t;
    out.theta_step  = theta_step;
    out.centroid    = ctr;
end

function s = sanitize_filename(s)
    s = regexprep(s, '[^\w\-]', '_');
end