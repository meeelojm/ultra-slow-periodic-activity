%% =========================================================================
% compare_active_vs_inactive_interhem_xcorr_from_tail.m
%
% SELF-CONTAINED SCRIPT
%
% What it does:
%   1) Recursively finds recordings that contain:
%        - tail_quick.mat
%        - dffs_repact_respcells.mat
%   2) Loads neural dff, fps, and neuron positions
%   3) Extracts tail bout onsets
%   4) Uses sliding windows to rank periods from least active to most active
%      based on tail bout onsets
%   5) Selects:
%        - least active 20-min window
%        - most active 20-min window
%   6) Splits neurons into left vs right hemisphere using x position
%   7) Builds left and right population traces in each selected window
%   8) Computes interhemispheric cross-correlation metrics:
%        - zero-lag correlation
%        - peak cross-correlation
%        - lag at peak
%        - area under xcorr near zero lag
%   9) Builds summary table, paired signrank statistics, and summary figure
%
% IMPORTANT:
%   - This script quantifies interhemispheric coupling via cross-correlation,
%     not spectral coherence in the strict signal-processing sense.
%   - It assumes the first position dimension corresponds to left-right.
%   - If no usable neuron positions are found in a recording, that recording
%     is skipped for interhemispheric analysis.
% =========================================================================

clear; clc; close all;

%% =========================
% 0) USER SETTINGS
% =========================

% Root folder that contains fish/session subfolders
rootFolder = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low';

% File names to search for
tailFileName = 'tail_quick.mat';
dffFileName  = 'dffs_repact_respcells.mat';

% Sliding-window settings for finding active/inactive stretches
winMin  = 20;   % window length in minutes
stepMin = 5;    % overlap step in minutes

% Neuron inclusion criterion inside selected windows
maxNanFracPerNeuron = 0.20;   % remove neurons with >20% NaNs in selected window

% Cross-correlation settings
maxLagSec     = 120;   % maximum lag (seconds) for xcorr
aucHalfWidthSec = 10;  % integrate xcorr within +/- this many seconds around zero lag

% Whether to z-score each population trace before xcorr
zscorePopulationTraces = true;

% Whether to save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_tail_interhem_xcorr_results');

% Figure style
set(0, 'DefaultFigureColor', [1 1 1]);

if saveOutputs && ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%% =========================
% 1) FIND MATCHING RECORDINGS
% =========================

fprintf('\nSearching for recordings under:\n%s\n\n', rootFolder);

tailFiles = dir(fullfile(rootFolder, '**', tailFileName));
dffFiles  = dir(fullfile(rootFolder, '**', dffFileName));

if isempty(tailFiles)
    error('No %s files found under rootFolder.', tailFileName);
end
if isempty(dffFiles)
    error('No %s files found under rootFolder.', dffFileName);
end

% Match recordings by folder
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
    error('No folders were found that contain BOTH %s and %s.', tailFileName, dffFileName);
end

fprintf('Found %d matched recordings.\n', numel(recordings));

%% =========================
% 2) PREALLOCATE RESULTS
% =========================

nR = numel(recordings);

R = struct();
R.recordingLabel          = cell(nR,1);
R.folder                  = cell(nR,1);

R.nNeurons_total          = nan(nR,1);
R.nNeurons_used_low       = nan(nR,1);
R.nNeurons_used_high      = nan(nR,1);
R.nNeurons_used_shared    = nan(nR,1);

R.nLeft                   = nan(nR,1);
R.nRight                  = nan(nR,1);

R.fps                     = nan(nR,1);
R.nFrames                 = nan(nR,1);
R.durationMin             = nan(nR,1);

R.low_startFrame          = nan(nR,1);
R.low_endFrame            = nan(nR,1);
R.high_startFrame         = nan(nR,1);
R.high_endFrame           = nan(nR,1);

R.low_boutCount           = nan(nR,1);
R.high_boutCount          = nan(nR,1);
R.low_boutRate_perMin     = nan(nR,1);
R.high_boutRate_perMin    = nan(nR,1);

R.low_xcorr_zeroLag       = nan(nR,1);
R.high_xcorr_zeroLag      = nan(nR,1);

R.low_xcorr_peak          = nan(nR,1);
R.high_xcorr_peak         = nan(nR,1);

R.low_xcorr_peakLagSec    = nan(nR,1);
R.high_xcorr_peakLagSec   = nan(nR,1);

R.low_xcorr_aucNearZero   = nan(nR,1);
R.high_xcorr_aucNearZero  = nan(nR,1);

R.low_left_mean           = nan(nR,1);
R.low_right_mean          = nan(nR,1);
R.high_left_mean          = nan(nR,1);
R.high_right_mean         = nan(nR,1);

allLowXcorr  = cell(nR,1);
allHighXcorr = cell(nR,1);
allLagsSec   = cell(nR,1);

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
    % 3.1 Load DFF, FPS, and positions
        % ---------------------------------------------------------------------
    S_dff = load(recordings(r).dffPath);
    [dff, fps, pos] = extract_dff_fps_positions(S_dff);
    
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
        fprintf('Transposing dff because rows > columns.\n');
        dff = dff.';
    end
    
    nNeurons = size(dff,1);
    nFrames  = size(dff,2);
    recDurMin = nFrames / fps / 60;
    
    R.nNeurons_total(r) = nNeurons;
    R.fps(r)            = fps;
    R.nFrames(r)        = nFrames;
    R.durationMin(r)    = recDurMin;
    
    fprintf('  DFF: %d neurons x %d frames | fps = %.4f | duration = %.2f min\n', ...
        nNeurons, nFrames, fps, recDurMin);
    
    % If positions were not found in dff file, try metadata / Fall.mat
    if isempty(pos) || size(pos,1) ~= nNeurons || size(pos,2) < 1
        pos = [];
    
        if isfield(recordings(r), 'metaPath') && ~isempty(recordings(r).metaPath) && isfile(recordings(r).metaPath)
            S_meta = load(recordings(r).metaPath);
            pos = extract_positions_only(S_meta, nNeurons);
            if ~isempty(pos)
                fprintf('  Loaded positions from metadata_multimodal.mat\n');
            end
        end
    
        if (isempty(pos) || size(pos,1) ~= nNeurons) && ...
           isfield(recordings(r), 'fallPath') && ~isempty(recordings(r).fallPath) && isfile(recordings(r).fallPath)
            S_fall = load(recordings(r).fallPath);
            pos = extract_positions_only(S_fall, nNeurons);
            if ~isempty(pos)
                fprintf('  Loaded positions from Fall.mat\n');
            end
        end
    end
    
    fprintf('  Position matrix size: [%d x %d]\n', size(pos,1), size(pos,2));

    dff = double(dff);

    % Make sure neurons x frames
    if size(dff,1) > size(dff,2)
        fprintf('Transposing dff because rows > columns.\n');
        dff = dff.';
    end

    nNeurons = size(dff,1);
    nFrames  = size(dff,2);
    recDurMin = nFrames / fps / 60;

    R.nNeurons_total(r) = nNeurons;
    R.fps(r)            = fps;
    R.nFrames(r)        = nFrames;
    R.durationMin(r)    = recDurMin;

    fprintf('  DFF: %d neurons x %d frames | fps = %.4f | duration = %.2f min\n', ...
        nNeurons, nFrames, fps, recDurMin);

    if isempty(pos) || size(pos,1) ~= nNeurons || size(pos,2) < 1
        warning('Could not find usable neuron positions in %s. Skipping.', recordings(r).dffPath);
        continue;
    end

    % ---------------------------------------------------------------------
    % 3.2 Load tail data and extract bout onsets
    % ---------------------------------------------------------------------
    S_tail = load(recordings(r).tailPath);
    [boutOnsetsFrames, tailActivity] = extract_bout_onsets_from_tail(S_tail, nFrames, fps); %#ok<NASGU>

    if isempty(boutOnsetsFrames)
        warning('No bout onsets could be extracted for %s. Skipping.', recordings(r).tailPath);
        continue;
    end

    boutOnsetsFrames = boutOnsetsFrames(:);
    boutOnsetsFrames = boutOnsetsFrames(~isnan(boutOnsetsFrames));
    boutOnsetsFrames = round(boutOnsetsFrames);
    boutOnsetsFrames = unique(boutOnsetsFrames);
    boutOnsetsFrames = boutOnsetsFrames(boutOnsetsFrames >= 1 & boutOnsetsFrames <= nFrames);

    fprintf('  Extracted %d bout onsets.\n', numel(boutOnsetsFrames));

    % ---------------------------------------------------------------------
    % 3.3 Rank sliding windows by activity
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
    % 3.4 Select neurons with acceptable NaN fraction in both windows
    % ---------------------------------------------------------------------
    dff_low  = dff(:, low_startFrame:low_endFrame);
    dff_high = dff(:, high_startFrame:high_endFrame);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;

    R.nNeurons_used_low(r)  = sum(keepLow);
    R.nNeurons_used_high(r) = sum(keepHigh);

    keepShared = keepLow & keepHigh & isfinite(pos(:,1)) & ~all(isnan(dff),2);
    R.nNeurons_used_shared(r) = sum(keepShared);

    if sum(keepShared) < 4
        warning('Too few shared usable neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    % ---------------------------------------------------------------------
    % 3.5 Split left and right hemisphere by x position
    % ---------------------------------------------------------------------
    posShared = pos(keepShared, :);
    dffShared = dff(keepShared, :);
    
    % Use the left-right coordinate
    xcoord = posShared(:,2);
    
    % remove non-finite entries just in case
    validPos = isfinite(xcoord);
    if sum(validPos) < 4
        warning('Too few valid x positions for hemisphere split in %s. Skipping.', recordings(r).folder);
        continue;
    end
    
    % k-means on 1D x-coordinate
    [idxK, C] = kmeans(xcoord(validPos), 2, ...
        'Replicates', 20, ...
        'MaxIter', 1000);
    
    % assign smaller center = left, larger center = right
    [~, order] = sort(C(:), 'ascend');
    leftLabel  = order(1);
    rightLabel = order(2);
    
    isLeft  = false(size(xcoord));
    isRight = false(size(xcoord));
    
    tmpLeft  = idxK == leftLabel;
    tmpRight = idxK == rightLabel;
    
    isLeft(validPos)  = tmpLeft;
    isRight(validPos) = tmpRight;

    % discard exact midline neurons if any
    nLeft  = sum(isLeft);
    nRight = sum(isRight);

    R.nLeft(r)  = nLeft;
    R.nRight(r) = nRight;

    if nLeft < 2 || nRight < 2
        warning('Not enough left/right neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    fprintf('  Shared neurons for xcorr: left = %d, right = %d\n', nLeft, nRight);

    % ---------------------------------------------------------------------
    % 3.6 Build population traces in low/high windows
    % ---------------------------------------------------------------------
    leftTrace_low   = mean(dffShared(isLeft,  low_startFrame:low_endFrame),  1, 'omitnan');
    rightTrace_low  = mean(dffShared(isRight, low_startFrame:low_endFrame), 1, 'omitnan');

    leftTrace_high  = mean(dffShared(isLeft,  high_startFrame:high_endFrame),  1, 'omitnan');
    rightTrace_high = mean(dffShared(isRight, high_startFrame:high_endFrame), 1, 'omitnan');

    % fill any residual NaNs in population traces
    leftTrace_low   = fillmissing(leftTrace_low(:),   'linear', 'EndValues', 'nearest');
    rightTrace_low  = fillmissing(rightTrace_low(:),  'linear', 'EndValues', 'nearest');
    leftTrace_high  = fillmissing(leftTrace_high(:),  'linear', 'EndValues', 'nearest');
    rightTrace_high = fillmissing(rightTrace_high(:), 'linear', 'EndValues', 'nearest');

    if any(~isfinite(leftTrace_low)) || any(~isfinite(rightTrace_low)) || ...
       any(~isfinite(leftTrace_high)) || any(~isfinite(rightTrace_high))
        warning('Population traces contain non-finite values in %s. Skipping.', recordings(r).folder);
        continue;
    end

    R.low_left_mean(r)   = mean(leftTrace_low,  'omitnan');
    R.low_right_mean(r)  = mean(rightTrace_low, 'omitnan');
    R.high_left_mean(r)  = mean(leftTrace_high, 'omitnan');
    R.high_right_mean(r) = mean(rightTrace_high,'omitnan');

    % ---------------------------------------------------------------------
    % 3.7 Cross-correlation metrics
    % ---------------------------------------------------------------------
    maxLagFrames = round(maxLagSec * fps);
    aucHalfWidthFrames = round(aucHalfWidthSec * fps);

    lowMetrics = compute_interhem_xcorr_metrics( ...
        leftTrace_low, rightTrace_low, fps, maxLagFrames, aucHalfWidthFrames, zscorePopulationTraces);

    highMetrics = compute_interhem_xcorr_metrics( ...
        leftTrace_high, rightTrace_high, fps, maxLagFrames, aucHalfWidthFrames, zscorePopulationTraces);

    R.low_xcorr_zeroLag(r)      = lowMetrics.zeroLag;
    R.high_xcorr_zeroLag(r)     = highMetrics.zeroLag;

    R.low_xcorr_peak(r)         = lowMetrics.peakCorr;
    R.high_xcorr_peak(r)        = highMetrics.peakCorr;

    R.low_xcorr_peakLagSec(r)   = lowMetrics.peakLagSec;
    R.high_xcorr_peakLagSec(r)  = highMetrics.peakLagSec;

    R.low_xcorr_aucNearZero(r)  = lowMetrics.aucNearZero;
    R.high_xcorr_aucNearZero(r) = highMetrics.aucNearZero;

    allLowXcorr{r}  = lowMetrics.xc(:);
    allHighXcorr{r} = highMetrics.xc(:);
    allLagsSec{r}   = lowMetrics.lagsSec(:);

    fprintf('  Low  xcorr: zero=%.4f, peak=%.4f @ %.2f s, auc0=%.4f\n', ...
        lowMetrics.zeroLag, lowMetrics.peakCorr, lowMetrics.peakLagSec, lowMetrics.aucNearZero);
    fprintf('  High xcorr: zero=%.4f, peak=%.4f @ %.2f s, auc0=%.4f\n', ...
        highMetrics.zeroLag, highMetrics.peakCorr, highMetrics.peakLagSec, highMetrics.aucNearZero);
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used_low, R.nNeurons_used_high, R.nNeurons_used_shared, ...
    R.nLeft, R.nRight, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.low_startFrame, R.low_endFrame, R.high_startFrame, R.high_endFrame, ...
    R.low_boutCount, R.high_boutCount, R.low_boutRate_perMin, R.high_boutRate_perMin, ...
    R.low_xcorr_zeroLag, R.high_xcorr_zeroLag, ...
    R.low_xcorr_peak, R.high_xcorr_peak, ...
    R.low_xcorr_peakLagSec, R.high_xcorr_peakLagSec, ...
    R.low_xcorr_aucNearZero, R.high_xcorr_aucNearZero, ...
    R.low_left_mean, R.low_right_mean, R.high_left_mean, R.high_right_mean, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used_low','nNeurons_used_high','nNeurons_used_shared', ...
    'nLeft','nRight', ...
    'fps','nFrames','durationMin', ...
    'low_startFrame','low_endFrame','high_startFrame','high_endFrame', ...
    'low_boutCount','high_boutCount','low_boutRate_perMin','high_boutRate_perMin', ...
    'low_xcorr_zeroLag','high_xcorr_zeroLag', ...
    'low_xcorr_peak','high_xcorr_peak', ...
    'low_xcorr_peakLagSec','high_xcorr_peakLagSec', ...
    'low_xcorr_aucNearZero','high_xcorr_aucNearZero', ...
    'low_left_mean','low_right_mean','high_left_mean','high_right_mean'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_interhem_xcorr_active_vs_inactive.csv'));
    save(fullfile(outFolder, 'summary_interhem_xcorr_active_vs_inactive.mat'), ...
        'R', 'T', 'allLowXcorr', 'allHighXcorr', 'allLagsSec');
end

%% =========================
% 5) PAIRED STATISTICS
% =========================

statsSummary = struct();
statsSummary.zeroLag   = run_paired_signrank(R.low_xcorr_zeroLag,     R.high_xcorr_zeroLag);
statsSummary.peakCorr  = run_paired_signrank(R.low_xcorr_peak,        R.high_xcorr_peak);
statsSummary.peakLag   = run_paired_signrank(R.low_xcorr_peakLagSec,  R.high_xcorr_peakLagSec);
statsSummary.aucNear0  = run_paired_signrank(R.low_xcorr_aucNearZero, R.high_xcorr_aucNearZero);

disp('================ PAIRED SIGNRANK TESTS ================')
fprintf('Zero-lag xcorr:       n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.zeroLag.n, statsSummary.zeroLag.p, statsSummary.zeroLag.deltaMedian);
fprintf('Peak xcorr:           n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.peakCorr.n, statsSummary.peakCorr.p, statsSummary.peakCorr.deltaMedian);
fprintf('Peak lag (s):         n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.peakLag.n, statsSummary.peakLag.p, statsSummary.peakLag.deltaMedian);
fprintf('AUC near zero:        n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.aucNear0.n, statsSummary.aucNear0.p, statsSummary.aucNear0.deltaMedian);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_interhem_xcorr.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURE
% =========================

validRec = isfinite(R.low_xcorr_zeroLag) & isfinite(R.high_xcorr_zeroLag) & ...
           isfinite(R.low_xcorr_peak) & isfinite(R.high_xcorr_peak) & ...
           isfinite(R.low_xcorr_peakLagSec) & isfinite(R.high_xcorr_peakLagSec) & ...
           isfinite(R.low_xcorr_aucNearZero) & isfinite(R.high_xcorr_aucNearZero);

if any(validRec)

    xPair = [1 2];
    idxV  = find(validRec(:))';

    figX = figure('Name','Active vs inactive interhemispheric xcorr comparison', ...
        'Position',[100 100 1400 700]);
    tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact');

    % ---------------------------------------------------------------------
    % Zero-lag correlation
    % ---------------------------------------------------------------------
    nexttile; hold on;
    yLow = R.low_xcorr_zeroLag(validRec);
    yHigh = R.high_xcorr_zeroLag(validRec);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxV
        plot(xPair, [R.low_xcorr_zeroLag(i), R.high_xcorr_zeroLag(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Zero-lag cross-correlation');
    title(sprintf('Zero-lag\nn=%d, p=%.3g %s', ...
        statsSummary.zeroLag.n, statsSummary.zeroLag.p, p_to_stars(statsSummary.zeroLag.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.zeroLag.p);

    % ---------------------------------------------------------------------
    % Peak cross-correlation
    % ---------------------------------------------------------------------
    nexttile; hold on;
    yLow = R.low_xcorr_peak(validRec);
    yHigh = R.high_xcorr_peak(validRec);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxV
        plot(xPair, [R.low_xcorr_peak(i), R.high_xcorr_peak(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Peak cross-correlation');
    title(sprintf('Peak xcorr\nn=%d, p=%.3g %s', ...
        statsSummary.peakCorr.n, statsSummary.peakCorr.p, p_to_stars(statsSummary.peakCorr.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.peakCorr.p);

    % ---------------------------------------------------------------------
    % Peak lag
    % ---------------------------------------------------------------------
    nexttile; hold on;
    yLow = R.low_xcorr_peakLagSec(validRec);
    yHigh = R.high_xcorr_peakLagSec(validRec);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxV
        plot(xPair, [R.low_xcorr_peakLagSec(i), R.high_xcorr_peakLagSec(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Lag at peak xcorr (s)');
    title(sprintf('Peak lag\nn=%d, p=%.3g %s', ...
        statsSummary.peakLag.n, statsSummary.peakLag.p, p_to_stars(statsSummary.peakLag.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.peakLag.p);

    % ---------------------------------------------------------------------
    % Mean xcorr curves with SEM
    % ---------------------------------------------------------------------
    nexttile; hold on;

    validCurves = validRec & ~cellfun(@isempty, allLowXcorr) & ~cellfun(@isempty, allHighXcorr) & ~cellfun(@isempty, allLagsSec);

    if any(validCurves)
        refIdx = find(validCurves, 1, 'first');
        lagsSec = allLagsSec{refIdx}(:);

        % stack curves only if same length
        keepIdx = [];
        for i = find(validCurves(:))'
            if numel(allLowXcorr{i}) == numel(lagsSec) && numel(allHighXcorr{i}) == numel(lagsSec)
                keepIdx(end+1) = i; %#ok<SAGROW>
            end
        end

        lowMat = nan(numel(lagsSec), numel(keepIdx));
        highMat = nan(numel(lagsSec), numel(keepIdx));

        for k = 1:numel(keepIdx)
            i = keepIdx(k);
            lowMat(:,k)  = allLowXcorr{i}(:);
            highMat(:,k) = allHighXcorr{i}(:);
        end

        lowMean = mean(lowMat, 2, 'omitnan');
        highMean = mean(highMat, 2, 'omitnan');

        lowSEM = std(lowMat, 0, 2, 'omitnan') ./ sqrt(size(lowMat,2));
        highSEM = std(highMat, 0, 2, 'omitnan') ./ sqrt(size(highMat,2));

        fill([lagsSec; flipud(lagsSec)], ...
             [lowMean-lowSEM; flipud(lowMean+lowSEM)], ...
             'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([lagsSec; flipud(lagsSec)], ...
             [highMean-highSEM; flipud(highMean+highSEM)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        plot(lagsSec, lowMean,  'b', 'LineWidth', 1.8);
        plot(lagsSec, highMean, 'r', 'LineWidth', 1.8);

        xlabel('Lag (s)');
        ylabel('Cross-correlation');
        title('Mean interhem xcorr');
        legend({'Least active SEM','Most active SEM','Least active mean','Most active mean'}, ...
            'Location','best');
        xline(0, '--k');
    else
        text(0.1,0.5,'No valid xcorr curves','FontSize',12);
        axis off;
    end

    sgtitle('Interhemispheric cross-correlation: tail-defined inactive vs active windows');

    if saveOutputs
        saveas(figX, fullfile(outFolder, 'summary_interhem_xcorr_active_vs_inactive.png'));
    end
else
    warning('No valid recordings available for xcorr summary figure.');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% =========================
% LOCAL FUNCTIONS
% =========================

function [dff, fps, pos] = extract_dff_fps_positions(S)
% Robust extraction of:
%   dff : neurons x frames (or frames x neurons; caller can transpose)
%   fps : scalar
%   pos : neurons x >=1, with first column = left-right coordinate

    dff = [];
    fps = [];
    pos = [];

    vars = fieldnames(S);

    % ---------------- dff ----------------
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

    % ---------------- fps ----------------
    preferredFpsNames = {'fps','Fs','fs','frameRate','framerate','samplingRate'};
    for i = 1:numel(preferredFpsNames)
        nm = preferredFpsNames{i};
        if isfield(S, nm) && isnumeric(S.(nm)) && isscalar(S.(nm))
            fps = double(S.(nm));
            break;
        end
    end

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
                            break;
                        end
                    end
                end
            end
            if ~isempty(fps)
                break;
            end
        end
    end

    % ---------------- positions ----------------
    candidatePosNames = { ...
        'positions', 'position', 'xyz', 'coords', 'coordinates', ...
        'allPos', 'centroids', 'cellPos', 'cellPositions', 'roi_pos'};

    for i = 1:numel(candidatePosNames)
        nm = candidatePosNames{i};
        if isfield(S, nm) && isnumeric(S.(nm)) && ismatrix(S.(nm)) && size(S.(nm),2) >= 1
            pos = double(S.(nm));
            break;
        end
    end

    % try struct-style fields
    if isempty(pos)
        for i = 1:numel(vars)
            v = S.(vars{i});
            if isstruct(v)
                fns = fieldnames(v);
                for j = 1:numel(fns)
                    nm = lower(fns{j});
                    if contains(nm, 'pos') || contains(nm, 'coord') || contains(nm, 'xyz')
                        vv = v.(fns{j});
                        if isnumeric(vv) && ismatrix(vv) && size(vv,2) >= 1
                            pos = double(vv);
                            break;
                        end
                    end
                end
            end
            if ~isempty(pos)
                break;
            end
        end
    end

    % try suite2p stat-style centroid extraction if present
    if isempty(pos) && isfield(S, 'stat')
        st = S.stat;
        if isstruct(st) && ~isempty(st)
            try
                n = numel(st);
                x = nan(n,1);
                y = nan(n,1);
                for k = 1:n
                    if isfield(st(k),'med') && numel(st(k).med) >= 2
                        y(k) = st(k).med(1);
                        x(k) = st(k).med(2);
                    elseif isfield(st(k),'xpix') && isfield(st(k),'ypix')
                        x(k) = mean(st(k).xpix);
                        y(k) = mean(st(k).ypix);
                    end
                end
                if any(isfinite(x))
                    pos = [x y];
                end
            catch
            end
        end
    end
end

function [boutOnsetsFrames, tailActivity] = extract_bout_onsets_from_tail(S_tail, nFrames, fps)
% Robustly extract bout onsets from tail_quick.mat

    boutOnsetsFrames = [];
    tailActivity = [];

    vars = fieldnames(S_tail);

    % -------- direct bout onset variables in frames
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
                if ~isempty(x) && max(x) <= nFrames*1.2
                    boutOnsetsFrames = round(x);
                    return;
                end
            end
        end
    end

    % -------- direct bout onset variables in seconds
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

    % -------- binary bout train
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

    % -------- continuous tail activity variables
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

    thr = median(tailActivity, 'omitnan') + 2 * mad(tailActivity, 1);
    active = tailActivity > thr;

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

    [~, lowIdx]  = min(winTable.boutCount);
    [~, highIdx] = max(winTable.boutCount);
end

function M = compute_interhem_xcorr_metrics(leftTrace, rightTrace, fps, maxLagFrames, aucHalfWidthFrames, doZscore)

    leftTrace  = leftTrace(:);
    rightTrace = rightTrace(:);

    if doZscore
        leftTrace  = safe_zscore(leftTrace);
        rightTrace = safe_zscore(rightTrace);
    end

    [xc, lags] = xcorr(leftTrace, rightTrace, maxLagFrames, 'coeff');
    lagsSec = lags / fps;

    % zero-lag
    zeroIdx = find(lags == 0, 1, 'first');
    zeroLag = xc(zeroIdx);

    % peak xcorr
    [peakCorr, peakIdx] = max(xc);
    peakLagSec = lagsSec(peakIdx);

    % area near zero
    keep = abs(lags) <= aucHalfWidthFrames;
    aucNearZero = trapz(lagsSec(keep), xc(keep));

    M = struct();
    M.xc          = xc(:);
    M.lags        = lags(:);
    M.lagsSec     = lagsSec(:);
    M.zeroLag     = zeroLag;
    M.peakCorr    = peakCorr;
    M.peakLagSec  = peakLagSec;
    M.aucNearZero = aucNearZero;
end

function z = safe_zscore(x)
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd <= 0
        sd = 1;
    end
    z = (x - mu) / sd;
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

    hold(ax,'on');
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

function pos = extract_positions_only(S, nNeurons)
    pos = [];
    vars = fieldnames(S);

    candidatePosNames = { ...
        'positions', 'pos', 'xyz', 'coords', 'coordinates', ...
        'allPos', 'centroids', 'cellPos', 'cellPositions', 'roi_pos'};

    for i = 1:numel(candidatePosNames)
        nm = candidatePosNames{i};
        if isfield(S, nm) && isnumeric(S.(nm)) && ismatrix(S.(nm)) && size(S.(nm),2) >= 1
            tmp = double(S.(nm));
            if size(tmp,1) == nNeurons
                pos = tmp;
                return;
            end
        end
    end

    for i = 1:numel(vars)
        v = S.(vars{i});
        if isstruct(v)
            fns = fieldnames(v);
            for j = 1:numel(fns)
                nm = lower(fns{j});
                if contains(nm, 'pos') || contains(nm, 'coord') || contains(nm, 'xyz')
                    vv = v.(fns{j});
                    if isnumeric(vv) && ismatrix(vv) && size(vv,1) == nNeurons && size(vv,2) >= 1
                        pos = double(vv);
                        return;
                    end
                end
            end
        end
    end

    if isfield(S, 'stat')
        st = S.stat;
        if isstruct(st) && numel(st) == nNeurons
            try
                x = nan(nNeurons,1);
                y = nan(nNeurons,1);
                for k = 1:nNeurons
                    if isfield(st(k),'med') && numel(st(k).med) >= 2
                        y(k) = st(k).med(1);
                        x(k) = st(k).med(2);
                    elseif isfield(st(k),'xpix') && isfield(st(k),'ypix')
                        x(k) = mean(st(k).xpix);
                        y(k) = mean(st(k).ypix);
                    end
                end
                if any(isfinite(x))
                    pos = [x y];
                    return;
                end
            catch
            end
        end
    end
end