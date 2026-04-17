%% =========================================================================
% compare_active_vs_inactive_pc_loops_from_tail_v2.m
%
% Tail-defined least-active vs most-active windows:
% detect and compare loops in shared PC space.
%
% Main logic
%   1) Find least-active and most-active windows from tail bouts
%   2) Extract neural activity in those windows
%   3) Build ONE shared PCA space from [lowWindow , highWindow]
%   4) Segment loops separately in low and high windows
%   5) Plot loops in PC space, colored by condition (low vs high)
%   6) Quantify loop size and orientation
%
% Outputs
%   - per-recording loop metrics table
%   - paired signrank statistics
%   - per-recording figures
%   - summary figure
%
% Notes
%   - "Orientation" is quantified in two ways:
%       (a) signed area in the PC1-PC2 plane:
%             > 0 = CCW, < 0 = CW
%       (b) 3D loop plane normal from SVD
%   - "Size" is quantified as:
%       durationSec, meanRadius, maxRadius, pathLength
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

% PCA / loop settings
nPC = 3;

% Radius-shell loop segmentation
smoothRadiusSec = 0.5;   % moving average on radius
thrStdFactor    = 0.5;   % shell threshold = median + thrStdFactor*std
peakStdFactor   = 1.0;   % peak threshold = mean + peakStdFactor*std
minLoopDurSec   = 1.0;   % keep loops at least this long

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_pc_loops_results_v2');

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
R.recordingLabel          = cell(nR,1);
R.folder                  = cell(nR,1);

R.nNeurons_total          = nan(nR,1);
R.nNeurons_used_shared    = nan(nR,1);

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

R.nLoops_low              = nan(nR,1);
R.nLoops_high             = nan(nR,1);

R.low_loopDur_median      = nan(nR,1);
R.high_loopDur_median     = nan(nR,1);

R.low_loopRadius_median   = nan(nR,1);
R.high_loopRadius_median  = nan(nR,1);

R.low_loopPath_median     = nan(nR,1);
R.high_loopPath_median    = nan(nR,1);

R.low_fracCCW             = nan(nR,1);
R.high_fracCCW            = nan(nR,1);

allLoopsLow               = cell(nR,1);
allLoopsHigh              = cell(nR,1);

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
    % 3.4 Select shared usable neurons
    % ---------------------------------------------------------------------
    dff_low  = dff(:, low_startFrame:low_endFrame);
    dff_high = dff(:, high_startFrame:high_endFrame);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;
    keepUse  = keepLow & keepHigh & ~all(isnan(dff),2);

    dff_use = dff(keepUse,:);
    if size(dff_use,1) < 5
        warning('Too few shared usable neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    dff_use = fill_nans_rowwise(dff_use);
    R.nNeurons_used_shared(r) = size(dff_use,1);

    % ---------------------------------------------------------------------
    % 3.5 Extract low/high windows and build shared PCA space
    % ---------------------------------------------------------------------
    Xlow  = dff_use(:, low_startFrame:low_endFrame);
    Xhigh = dff_use(:, high_startFrame:high_endFrame);

    Xcat = [Xlow, Xhigh];
    XcatZ = zscore(Xcat, 0, 2);   % neuronwise across concatenated windows
    XcatZ(isnan(XcatZ)) = 0;

    [~, scoreCat, ~, ~, explained] = pca(XcatZ.', 'NumComponents', nPC);

    nLowT = size(Xlow,2);
    scoreLow  = scoreCat(1:nLowT,  1:3);
    scoreHigh = scoreCat(nLowT+1:end, 1:3);

    fprintf('  PCA explained variance PC1-3: %.1f / %.1f / %.1f %%\n', ...
        explained(1), explained(2), explained(3));

    % ---------------------------------------------------------------------
    % 3.6 Segment loops separately in low and high windows
    % ---------------------------------------------------------------------
    loopsLow = segment_loops_by_condition(scoreLow, fps, minLoopDurSec, ...
        smoothRadiusSec, thrStdFactor, peakStdFactor, 'low');

    loopsHigh = segment_loops_by_condition(scoreHigh, fps, minLoopDurSec, ...
        smoothRadiusSec, thrStdFactor, peakStdFactor, 'high');

    allLoopsLow{r}  = loopsLow;
    allLoopsHigh{r} = loopsHigh;

    % ---------------------------------------------------------------------
    % 3.7 Recording-level summaries
    % ---------------------------------------------------------------------
    [R.nLoops_low(r), ...
     R.low_loopDur_median(r), ...
     R.low_loopRadius_median(r), ...
     R.low_loopPath_median(r), ...
     R.low_fracCCW(r)] = summarize_loop_set(loopsLow);

    [R.nLoops_high(r), ...
     R.high_loopDur_median(r), ...
     R.high_loopRadius_median(r), ...
     R.high_loopPath_median(r), ...
     R.high_fracCCW(r)] = summarize_loop_set(loopsHigh);

    fprintf('  Low loops : n=%d | median dur=%.2f s | median radius=%.3f | median path=%.3f | frac CCW=%.2f\n', ...
        R.nLoops_low(r), R.low_loopDur_median(r), R.low_loopRadius_median(r), ...
        R.low_loopPath_median(r), R.low_fracCCW(r));
    fprintf('  High loops: n=%d | median dur=%.2f s | median radius=%.3f | median path=%.3f | frac CCW=%.2f\n', ...
        R.nLoops_high(r), R.high_loopDur_median(r), R.high_loopRadius_median(r), ...
        R.high_loopPath_median(r), R.high_fracCCW(r));

    % ---------------------------------------------------------------------
    % 3.8 Per-recording figure
    % ---------------------------------------------------------------------
    try
        figL = figure('Name',['PC loops ' R.recordingLabel{r}], ...
            'Position',[100 100 1300 900], 'Color','w');
        tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

        % PC1-PC2
        nexttile; hold on;
        plot(scoreLow(:,1),  scoreLow(:,2),  '-', 'Color',[0.80 0.85 1.00], 'LineWidth',0.5);
        plot(scoreHigh(:,1), scoreHigh(:,2), '-', 'Color',[1.00 0.85 0.85], 'LineWidth',0.5);
        plot_loop_set_2d(loopsLow,  [0 0.35 1]);
        plot_loop_set_2d(loopsHigh, [0.95 0.2 0.2]);
        xlabel('PC1'); ylabel('PC2'); axis equal; box off;
        title('Loops in PC1-PC2');
        legend({'Low traj','High traj','Low loops','High loops'}, 'Location','best');

        % PC1-PC2-PC3
        nexttile; hold on;
        plot3(scoreLow(:,1), scoreLow(:,2), scoreLow(:,3), '-', 'Color',[0.80 0.85 1.00], 'LineWidth',0.4);
        plot3(scoreHigh(:,1),scoreHigh(:,2),scoreHigh(:,3),'-', 'Color',[1.00 0.85 0.85], 'LineWidth',0.4);
        plot_loop_set_3d(loopsLow,  [0 0.35 1]);
        plot_loop_set_3d(loopsHigh, [0.95 0.2 0.2]);
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        grid on; view(40,25); box off;
        title('Loops in PC1-PC2-PC3');

        % Loop durations
        nexttile; hold on;
        if ~isempty(loopsLow)
            histogram([loopsLow.durationSec], 'FaceColor',[0 0.35 1], 'FaceAlpha',0.35);
        end
        if ~isempty(loopsHigh)
            histogram([loopsHigh.durationSec], 'FaceColor',[0.95 0.2 0.2], 'FaceAlpha',0.35);
        end
        xlabel('Loop duration (s)'); ylabel('Count'); box off;
        title('Loop durations');
        legend({'Least active','Most active'}, 'Location','best');

        % Signed area / orientation
        nexttile; hold on;
        yLow = [loopsLow.signedArea12]';
        yHigh = [loopsHigh.signedArea12]';
        if ~isempty(yLow) || ~isempty(yHigh)
            if ~isempty(yLow),  boxchart(ones(numel(yLow),1), yLow); end
            if ~isempty(yHigh), boxchart(ones(numel(yHigh),1)*2, yHigh); end
            xlim([0.5 2.5]);
            set(gca,'XTick',[1 2],'XTickLabel',{'Least active','Most active'});
            ylabel('Signed area in PC1-PC2');
        else
            text(0.5,0.5,'No loops detected','Units','normalized','HorizontalAlignment','center');
            axis off;
        end
        title('Loop orientation metric');

        sgtitle(sprintf('%s | PCA loops | PC var %.1f / %.1f / %.1f %%', ...
            R.recordingLabel{r}, explained(1), explained(2), explained(3)), ...
            'Interpreter','none');

        if saveOutputs
            saveas(figL, fullfile(outFolder, ['PC_loops_' sanitize_filename(R.recordingLabel{r}) '.png']));
        end
        close(figL);
    catch ME
        warning('Could not create loop figure for %s\n%s', recordings(r).folder, ME.message);
    end
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used_shared, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.low_startFrame, R.low_endFrame, R.high_startFrame, R.high_endFrame, ...
    R.low_boutCount, R.high_boutCount, R.low_boutRate_perMin, R.high_boutRate_perMin, ...
    R.nLoops_low, R.nLoops_high, ...
    R.low_loopDur_median, R.high_loopDur_median, ...
    R.low_loopRadius_median, R.high_loopRadius_median, ...
    R.low_loopPath_median, R.high_loopPath_median, ...
    R.low_fracCCW, R.high_fracCCW, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used_shared', ...
    'fps','nFrames','durationMin', ...
    'low_startFrame','low_endFrame','high_startFrame','high_endFrame', ...
    'low_boutCount','high_boutCount','low_boutRate_perMin','high_boutRate_perMin', ...
    'nLoops_low','nLoops_high', ...
    'low_loopDur_median','high_loopDur_median', ...
    'low_loopRadius_median','high_loopRadius_median', ...
    'low_loopPath_median','high_loopPath_median', ...
    'low_fracCCW','high_fracCCW'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_pc_loops_active_vs_inactive.csv'));
    save(fullfile(outFolder, 'summary_pc_loops_active_vs_inactive.mat'), ...
        'R', 'T', 'allLoopsLow', 'allLoopsHigh', '-v7.3');
end

%% =========================
% 5) PAIRED STATISTICS
% =========================

statsSummary = struct();
statsSummary.nLoops  = run_paired_signrank(R.nLoops_low, R.nLoops_high);
statsSummary.dur     = run_paired_signrank(R.low_loopDur_median, R.high_loopDur_median);
statsSummary.radius  = run_paired_signrank(R.low_loopRadius_median, R.high_loopRadius_median);
statsSummary.path    = run_paired_signrank(R.low_loopPath_median, R.high_loopPath_median);
statsSummary.fracCCW = run_paired_signrank(R.low_fracCCW, R.high_fracCCW);

disp('================ PAIRED SIGNRANK TESTS: PCA LOOPS ================')
fprintf('Number of loops:   n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.nLoops.n, statsSummary.nLoops.p, statsSummary.nLoops.deltaMedian);
fprintf('Loop duration:     n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.dur.n, statsSummary.dur.p, statsSummary.dur.deltaMedian);
fprintf('Loop mean radius:  n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.radius.n, statsSummary.radius.p, statsSummary.radius.deltaMedian);
fprintf('Loop path length:  n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.path.n, statsSummary.path.p, statsSummary.path.deltaMedian);
fprintf('Fraction CCW:      n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.fracCCW.n, statsSummary.fracCCW.p, statsSummary.fracCCW.deltaMedian);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_pc_loops.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURE
% =========================

validRec = isfinite(R.low_loopDur_median) & isfinite(R.high_loopDur_median) & ...
           isfinite(R.low_loopRadius_median) & isfinite(R.high_loopRadius_median) & ...
           isfinite(R.low_loopPath_median) & isfinite(R.high_loopPath_median);

if any(validRec)

    idxV = find(validRec(:))';
    xPair = [1 2];

    figS = figure('Name','PC loop summary', 'Position',[100 100 1400 600], 'Color','w');
    tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact');

    % Number of loops
    nexttile; hold on;
    yLow = R.nLoops_low(validRec);
    yHigh = R.nLoops_high(validRec);
    boxchart(ones(numel(yLow),1), yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.nLoops_low(i), R.nLoops_high(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Number of loops');
    title(sprintf('Loop count\nn=%d, p=%.3g %s', ...
        statsSummary.nLoops.n, statsSummary.nLoops.p, p_to_stars(statsSummary.nLoops.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.nLoops.p);

    % Duration
    nexttile; hold on;
    yLow = R.low_loopDur_median(validRec);
    yHigh = R.high_loopDur_median(validRec);
    boxchart(ones(numel(yLow),1), yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_loopDur_median(i), R.high_loopDur_median(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Median loop duration (s)');
    title(sprintf('Loop duration\nn=%d, p=%.3g %s', ...
        statsSummary.dur.n, statsSummary.dur.p, p_to_stars(statsSummary.dur.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.dur.p);

    % Radius
    nexttile; hold on;
    yLow = R.low_loopRadius_median(validRec);
    yHigh = R.high_loopRadius_median(validRec);
    boxchart(ones(numel(yLow),1), yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_loopRadius_median(i), R.high_loopRadius_median(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Median loop mean radius');
    title(sprintf('Loop size\nn=%d, p=%.3g %s', ...
        statsSummary.radius.n, statsSummary.radius.p, p_to_stars(statsSummary.radius.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.radius.p);

    % Path length
    nexttile; hold on;
    yLow = R.low_loopPath_median(validRec);
    yHigh = R.high_loopPath_median(validRec);
    boxchart(ones(numel(yLow),1), yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_loopPath_median(i), R.high_loopPath_median(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Median loop path length');
    title(sprintf('Loop path length\nn=%d, p=%.3g %s', ...
        statsSummary.path.n, statsSummary.path.p, p_to_stars(statsSummary.path.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.path.p);

    sgtitle('Loops in shared PC space: least active vs most active');

    if saveOutputs
        saveas(figS, fullfile(outFolder, 'summary_pc_loops_active_vs_inactive.png'));
    end
else
    warning('No valid recordings available for PC loop summary figure.');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% ================== HELPERS ==================

function loops = segment_loops_by_condition(score, fps, minLoopDurSec, smoothRadiusSec, thrStdFactor, peakStdFactor, condLabel)
% Detect loops from one condition in PCA score space.
%
% score : T x 3
% loops(k) fields:
%   idx, xyz, timeSec, durationSec, meanRadius, maxRadius, pathLength,
%   signedArea12, orientation12, normal3D, condition

    loops = struct('idx',{},'xyz',{},'timeSec',{},'durationSec',{}, ...
                   'meanRadius',{},'maxRadius',{},'pathLength',{}, ...
                   'signedArea12',{},'orientation12',{},'normal3D',{}, ...
                   'condition',{});

    if isempty(score) || size(score,1) < 10 || size(score,2) < 3
        return;
    end

    rt = vecnorm(score(:,1:3), 2, 2);
    w = max(3, round(smoothRadiusSec * fps));
    rt_smooth = movmean(rt, w);

    th = median(rt_smooth, 'omitnan') + thrStdFactor * std(rt_smooth, 'omitnan');
    inShell = rt_smooth < th;

    shellLbl = bwlabel(inShell);
    reEntry  = find(diff(shellLbl) > 0) + 1;

    [~,locs] = findpeaks(rt_smooth, ...
        'MinPeakHeight', mean(rt_smooth,'omitnan') + peakStdFactor * std(rt_smooth,'omitnan'));

    if numel(reEntry) < 2
        return;
    end

    keepSegments = false(numel(reEntry)-1,1);
    for k = 1:numel(reEntry)-1
        idx = reEntry(k):reEntry(k+1);
        keepSegments(k) = any(locs >= idx(1) & locs <= idx(end));
    end

    c = 0;
    for k = 1:numel(reEntry)-1
        idx = reEntry(k):reEntry(k+1);
        durSec = numel(idx) / fps;

        if durSec < minLoopDurSec || ~keepSegments(k)
            continue;
        end

        xyz = score(idx,1:3);
        ctr = mean(xyz, 1, 'omitnan');
        rel = xyz - ctr;

        meanRadius = mean(vecnorm(rel,2,2), 'omitnan');
        maxRadius  = max(vecnorm(rel,2,2), [], 'omitnan');
        pathLength = sum(vecnorm(diff(xyz,1,1),2,2), 'omitnan');

        signedArea12 = polygon_signed_area(rel(:,1), rel(:,2));
        if signedArea12 > 0
            orientation12 = "CCW";
        elseif signedArea12 < 0
            orientation12 = "CW";
        else
            orientation12 = "flat";
        end

        normal3D = estimate_loop_normal(rel);

        c = c + 1;
        loops(c).idx          = idx;
        loops(c).xyz          = xyz;
        loops(c).timeSec      = (idx(:)-1) / fps;
        loops(c).durationSec  = durSec;
        loops(c).meanRadius   = meanRadius;
        loops(c).maxRadius    = maxRadius;
        loops(c).pathLength   = pathLength;
        loops(c).signedArea12 = signedArea12;
        loops(c).orientation12 = orientation12;
        loops(c).normal3D     = normal3D;
        loops(c).condition    = condLabel;
    end
end

function [nLoops, durMed, radiusMed, pathMed, fracCCW] = summarize_loop_set(loops)
    if isempty(loops)
        nLoops = 0;
        durMed = NaN;
        radiusMed = NaN;
        pathMed = NaN;
        fracCCW = NaN;
        return;
    end

    nLoops    = numel(loops);
    durMed    = median([loops.durationSec], 'omitnan');
    radiusMed = median([loops.meanRadius], 'omitnan');
    pathMed   = median([loops.pathLength], 'omitnan');

    ori = string({loops.orientation12});
    fracCCW = mean(ori == "CCW", 'omitnan');
end

function plot_loop_set_2d(loops, col)
    for k = 1:numel(loops)
        xyz = loops(k).xyz;
        plot(xyz(:,1), xyz(:,2), '-', 'Color', col, 'LineWidth', 1.5);
    end
end

function plot_loop_set_3d(loops, col)
    for k = 1:numel(loops)
        xyz = loops(k).xyz;
        plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-', 'Color', col, 'LineWidth', 1.5);
    end
end

function A = polygon_signed_area(x, y)
    x = x(:); y = y(:);
    if numel(x) < 3
        A = 0;
        return;
    end
    x2 = [x; x(1)];
    y2 = [y; y(1)];
    A = 0.5 * sum(x2(1:end-1).*y2(2:end) - x2(2:end).*y2(1:end-1), 'omitnan');
end

function n = estimate_loop_normal(rel)
    if size(rel,1) < 3
        n = [NaN NaN NaN];
        return;
    end
    [~,~,V] = svd(rel, 'econ');
    n = V(:,end)';
    n = n / norm(n);
end

function [dff, fps] = extract_dff_and_fps(S)
    dff = [];
    fps = [];

    vars = fieldnames(S);

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
                            return;
                        end
                    end
                end
            end
        end
    end
end

function [boutOnsetsFrames, tailActivity] = extract_bout_onsets_from_tail(S_tail, nFrames, fps)
    boutOnsetsFrames = [];
    tailActivity = [];

    vars = fieldnames(S_tail);

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

    if numel(tailActivity) ~= nFrames
        tailActivity = try_resample_vector(tailActivity, nFrames);
    end

    tailActivity = fillmissing(tailActivity, 'linear', 'EndValues', 'nearest');
    tailActivity = smoothdata(tailActivity, 'gaussian', 5);

    thr = median(tailActivity, 'omitnan') + 2 * mad(tailActivity, 1);
    active = tailActivity > thr;

    minBoutDurFrames = max(2, round(0.10 * fps));
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

function s = sanitize_filename(s)
    s = char(s);
    s = regexprep(s, '[<>:"/\\|?*\s]+', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_+|_+$', '');
end