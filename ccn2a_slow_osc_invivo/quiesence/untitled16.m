%% =========================================================================
% compare_active_vs_inactive_events_from_tail_lalamethod.m
%
% Tail-defined least-active vs most-active windows:
% compare neural event statistics using the LL significant-transient method
% (dynamic threshold, Gaussian baseline-noise model, 95% confidence)
%
% Requires:
%   - tail_quick.mat
%   - dffs_repact_respcells.mat
%   - Statistics and Image Processing toolboxes
%
% Output:
%   - per-recording summary table
%   - paired signrank statistics
%   - one summary figure
%
% Notes:
%   - This implements the LL event detector logic discussed earlier.
%   - Event detection is run on the full recording, after resampling traces
%     to a common fps (default 4 Hz), then window metrics are extracted from
%     the tail-defined inactive vs active windows.
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

% LL detector parameters
LL = struct();
LL.targetFps = 4;                % resample traces before event detection
LL.tauDecay = 2;                 % seconds
LL.confCutOff = 95;              % confidence cutoff
LL.cri_dur_sec = 0.5;            % merge/delete threshold in seconds
LL.tracesArePercent = true;      % IMPORTANT: set false if dff is already fractional, not percent
LL.plotFlag = 0;
makePerFishRasterQC = true;

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_tail_events_results_LL');

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
        recordings(end+1).folder  = thisFolder; %#ok<SAGROW>
        recordings(end).tailPath  = fullfile(tailFiles(i).folder, tailFiles(i).name);
        recordings(end).dffPath   = dffMap(thisFolder);
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
R.nNeurons_used_low       = nan(nR,1);
R.nNeurons_used_high      = nan(nR,1);
R.nNeurons_used_shared    = nan(nR,1);

R.fps_orig                = nan(nR,1);
R.fps_event               = nan(nR,1);
R.nFrames_orig            = nan(nR,1);
R.nFrames_event           = nan(nR,1);
R.durationMin             = nan(nR,1);

R.low_startFrame_orig     = nan(nR,1);
R.low_endFrame_orig       = nan(nR,1);
R.high_startFrame_orig    = nan(nR,1);
R.high_endFrame_orig      = nan(nR,1);

R.low_startFrame_event    = nan(nR,1);
R.low_endFrame_event      = nan(nR,1);
R.high_startFrame_event   = nan(nR,1);
R.high_endFrame_event     = nan(nR,1);

R.low_boutCount           = nan(nR,1);
R.high_boutCount          = nan(nR,1);
R.low_boutRate_perMin     = nan(nR,1);
R.high_boutRate_perMin    = nan(nR,1);

% Event metrics
R.low_eventRate_mean      = nan(nR,1);
R.high_eventRate_mean     = nan(nR,1);

R.low_eventRate_median    = nan(nR,1);
R.high_eventRate_median   = nan(nR,1);

R.low_totalEventCount     = nan(nR,1);
R.high_totalEventCount    = nan(nR,1);

R.low_fracActiveNeurons   = nan(nR,1);
R.high_fracActiveNeurons  = nan(nR,1);

R.low_eventAmp_mean       = nan(nR,1);
R.high_eventAmp_mean      = nan(nR,1);

R.low_eventDur_mean       = nan(nR,1);
R.high_eventDur_mean      = nan(nR,1);

R.low_eventAUC_mean       = nan(nR,1);
R.high_eventAUC_mean      = nan(nR,1);

allLowRates               = cell(nR,1);
allHighRates              = cell(nR,1);
allLowEventAmp            = cell(nR,1);
allHighEventAmp           = cell(nR,1);
allLowEventDur            = cell(nR,1);
allHighEventDur           = cell(nR,1);
allLowEventAUC            = cell(nR,1);
allHighEventAUC           = cell(nR,1);

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
    R.fps_orig(r)       = fps;
    R.nFrames_orig(r)   = nFrames;
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
    % 3.3 Rank windows by tail activity (original fps / original frames)
    % ---------------------------------------------------------------------
    winFrames_orig  = max(1, round(winMin  * 60 * fps));
    stepFrames_orig = max(1, round(stepMin * 60 * fps));

    if nFrames < winFrames_orig
        warning('Recording shorter than one analysis window. Skipping.');
        continue;
    end

    [winTable, lowIdx, highIdx] = rank_windows_by_bout_activity( ...
        boutOnsetsFrames, nFrames, winFrames_orig, stepFrames_orig, fps);

    if isempty(winTable)
        warning('Could not build sliding windows for %s. Skipping.', recordings(r).folder);
        continue;
    end

    low_startFrame_orig  = winTable.startFrame(lowIdx);
    low_endFrame_orig    = winTable.endFrame(lowIdx);
    high_startFrame_orig = winTable.startFrame(highIdx);
    high_endFrame_orig   = winTable.endFrame(highIdx);

    R.low_startFrame_orig(r)  = low_startFrame_orig;
    R.low_endFrame_orig(r)    = low_endFrame_orig;
    R.high_startFrame_orig(r) = high_startFrame_orig;
    R.high_endFrame_orig(r)   = high_endFrame_orig;

    R.low_boutCount(r)        = winTable.boutCount(lowIdx);
    R.high_boutCount(r)       = winTable.boutCount(highIdx);
    R.low_boutRate_perMin(r)  = winTable.boutRatePerMin(lowIdx);
    R.high_boutRate_perMin(r) = winTable.boutRatePerMin(highIdx);

    fprintf('  Least active window: frames %d-%d | bouts = %d | %.3f bouts/min\n', ...
        low_startFrame_orig, low_endFrame_orig, R.low_boutCount(r), R.low_boutRate_perMin(r));
    fprintf('  Most active window : frames %d-%d | bouts = %d | %.3f bouts/min\n', ...
        high_startFrame_orig, high_endFrame_orig, R.high_boutCount(r), R.high_boutRate_perMin(r));

    % ---------------------------------------------------------------------
    % 3.4 Select usable neurons
    % ---------------------------------------------------------------------
    dff_low  = dff(:, low_startFrame_orig:low_endFrame_orig);
    dff_high = dff(:, high_startFrame_orig:high_endFrame_orig);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;

    keepUse = keepLow & keepHigh & ~all(isnan(dff),2);

    dff_use = dff(keepUse,:);
    R.nNeurons_used_low(r)    = sum(keepLow);
    R.nNeurons_used_high(r)   = sum(keepHigh);
    R.nNeurons_used_shared(r) = size(dff_use,1);

    if size(dff_use,1) < 1
        warning('No shared usable neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    dff_use = fill_nans_rowwise(dff_use);

    % ---------------------------------------------------------------------
    % 3.5 LL event detection on full recording
    % ---------------------------------------------------------------------
    fprintf('  Running LL event detection...\n');

    try
        [eve_fra_cel, dff_dow_fra_cel, parLL] = ext_eve_LL(LL, fps, dff_use);  % frames x cells, frames x cells
    catch ME
        warning('LL event detection failed in %s\n%s', recordings(r).folder, ME.message);
        continue;
    end

    if isempty(eve_fra_cel) || isempty(dff_dow_fra_cel)
        warning('Empty LL outputs in %s. Skipping.', recordings(r).folder);
        continue;
    end

    if size(eve_fra_cel,1) ~= size(dff_dow_fra_cel,1) || size(eve_fra_cel,2) ~= size(dff_dow_fra_cel,2)
        warning('LL raster / resampled DFF size mismatch in %s. Skipping.', recordings(r).folder);
        continue;
    end

    fps_event = parLL.fps;
    R.fps_event(r)    = fps_event;
    R.nFrames_event(r)= size(eve_fra_cel,1);

    fprintf('  LL raster: %d frames x %d neurons | event fps = %.4f\n', ...
        size(eve_fra_cel,1), size(eve_fra_cel,2), fps_event);

    % ---------------------------------------------------------------------
    % 3.6 Convert low/high windows to event-frame indices
    % ---------------------------------------------------------------------
    low_startFrame_event  = max(1, round((low_startFrame_orig  - 1) * fps_event / fps) + 1);
    low_endFrame_event    = min(size(eve_fra_cel,1), round((low_endFrame_orig  - 1) * fps_event / fps) + 1);
    high_startFrame_event = max(1, round((high_startFrame_orig - 1) * fps_event / fps) + 1);
    high_endFrame_event   = min(size(eve_fra_cel,1), round((high_endFrame_orig - 1) * fps_event / fps) + 1);

    if low_endFrame_event <= low_startFrame_event || high_endFrame_event <= high_startFrame_event
        warning('Event-frame mapped window invalid in %s. Skipping.', recordings(r).folder);
        continue;
    end

    R.low_startFrame_event(r)  = low_startFrame_event;
    R.low_endFrame_event(r)    = low_endFrame_event;
    R.high_startFrame_event(r) = high_startFrame_event;
    R.high_endFrame_event(r)   = high_endFrame_event;

    low_idx_event  = low_startFrame_event:low_endFrame_event;
    high_idx_event = high_startFrame_event:high_endFrame_event;

    % ---------------------------------------------------------------------
    % 3.7 Quantify events in low/high windows
    % ---------------------------------------------------------------------
    [lowMetrics, lowAmp, lowDur, lowAUC]   = summarize_LL_events_in_window(eve_fra_cel, dff_dow_fra_cel, low_idx_event, fps_event);
    [highMetrics, highAmp, highDur, highAUC] = summarize_LL_events_in_window(eve_fra_cel, dff_dow_fra_cel, high_idx_event, fps_event);

    R.low_eventRate_mean(r)     = lowMetrics.rate_mean;
    R.high_eventRate_mean(r)    = highMetrics.rate_mean;

    R.low_eventRate_median(r)   = lowMetrics.rate_median;
    R.high_eventRate_median(r)  = highMetrics.rate_median;

    R.low_totalEventCount(r)    = lowMetrics.total_count;
    R.high_totalEventCount(r)   = highMetrics.total_count;

    R.low_fracActiveNeurons(r)  = lowMetrics.frac_active_neurons;
    R.high_fracActiveNeurons(r) = highMetrics.frac_active_neurons;

    R.low_eventAmp_mean(r)      = lowMetrics.amp_mean;
    R.high_eventAmp_mean(r)     = highMetrics.amp_mean;

    R.low_eventDur_mean(r)      = lowMetrics.dur_mean;
    R.high_eventDur_mean(r)     = highMetrics.dur_mean;

    R.low_eventAUC_mean(r)      = lowMetrics.auc_mean;
    R.high_eventAUC_mean(r)     = highMetrics.auc_mean;

    allLowRates{r}      = lowMetrics.rates_per_neuron;
    allHighRates{r}     = highMetrics.rates_per_neuron;
    allLowEventAmp{r}   = lowAmp;
    allHighEventAmp{r}  = highAmp;
    allLowEventDur{r}   = lowDur;
    allHighEventDur{r}  = highDur;
    allLowEventAUC{r}   = lowAUC;
    allHighEventAUC{r}  = highAUC;

    fprintf('  Low  events: total=%d, mean rate=%.4g, frac active=%.3f, mean dur=%.4g s, mean AUC=%.4g\n', ...
        lowMetrics.total_count, lowMetrics.rate_mean, lowMetrics.frac_active_neurons, ...
        lowMetrics.dur_mean, lowMetrics.auc_mean);
    fprintf('  High events: total=%d, mean rate=%.4g, frac active=%.3f, mean dur=%.4g s, mean AUC=%.4g\n', ...
        highMetrics.total_count, highMetrics.rate_mean, highMetrics.frac_active_neurons, ...
        highMetrics.dur_mean, highMetrics.auc_mean);

    % ---------------------------------------------------------------------
    % 3.8 QC raster figure for this fish: low vs high window
    % ---------------------------------------------------------------------
    if makePerFishRasterQC
        try
            figQC = make_LL_window_raster_QC( ...
                eve_fra_cel, dff_dow_fra_cel, ...
                low_idx_event, high_idx_event, ...
                fps_event, R.recordingLabel{r});

            if saveOutputs
                qcName = sprintf('QC_raster_%s.png', sanitize_filename(R.recordingLabel{r}));
                saveas(figQC, fullfile(outFolder, qcName));
            end
            close(figQC);
        catch ME
            warning('Could not create QC raster for %s\n%s', recordings(r).folder, ME.message);
        end
    end
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used_low, R.nNeurons_used_high, R.nNeurons_used_shared, ...
    R.fps_orig, R.fps_event, R.nFrames_orig, R.nFrames_event, R.durationMin, ...
    R.low_startFrame_orig, R.low_endFrame_orig, R.high_startFrame_orig, R.high_endFrame_orig, ...
    R.low_startFrame_event, R.low_endFrame_event, R.high_startFrame_event, R.high_endFrame_event, ...
    R.low_boutCount, R.high_boutCount, R.low_boutRate_perMin, R.high_boutRate_perMin, ...
    R.low_eventRate_mean, R.high_eventRate_mean, ...
    R.low_eventRate_median, R.high_eventRate_median, ...
    R.low_totalEventCount, R.high_totalEventCount, ...
    R.low_fracActiveNeurons, R.high_fracActiveNeurons, ...
    R.low_eventAmp_mean, R.high_eventAmp_mean, ...
    R.low_eventDur_mean, R.high_eventDur_mean, ...
    R.low_eventAUC_mean, R.high_eventAUC_mean, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used_low','nNeurons_used_high','nNeurons_used_shared', ...
    'fps_orig','fps_event','nFrames_orig','nFrames_event','durationMin', ...
    'low_startFrame_orig','low_endFrame_orig','high_startFrame_orig','high_endFrame_orig', ...
    'low_startFrame_event','low_endFrame_event','high_startFrame_event','high_endFrame_event', ...
    'low_boutCount','high_boutCount','low_boutRate_perMin','high_boutRate_perMin', ...
    'low_eventRate_mean','high_eventRate_mean', ...
    'low_eventRate_median','high_eventRate_median', ...
    'low_totalEventCount','high_totalEventCount', ...
    'low_fracActiveNeurons','high_fracActiveNeurons', ...
    'low_eventAmp_mean','high_eventAmp_mean', ...
    'low_eventDur_mean','high_eventDur_mean', ...
    'low_eventAUC_mean','high_eventAUC_mean'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_events_active_vs_inactive_LL.csv'));
    save(fullfile(outFolder, 'summary_events_active_vs_inactive_LL.mat'), ...
        'R', 'T', ...
        'allLowRates', 'allHighRates', ...
        'allLowEventAmp', 'allHighEventAmp', ...
        'allLowEventDur', 'allHighEventDur', ...
        'allLowEventAUC', 'allHighEventAUC', ...
        '-v7.3');
end

%% =========================
% 5) PAIRED STATISTICS
% =========================

statsSummary = struct();
statsSummary.rateMean   = run_paired_signrank(R.low_eventRate_mean,    R.high_eventRate_mean);
statsSummary.rateMedian = run_paired_signrank(R.low_eventRate_median,  R.high_eventRate_median);
statsSummary.totalCount = run_paired_signrank(R.low_totalEventCount,   R.high_totalEventCount);
statsSummary.fracActive = run_paired_signrank(R.low_fracActiveNeurons, R.high_fracActiveNeurons);
statsSummary.ampMean    = run_paired_signrank(R.low_eventAmp_mean,     R.high_eventAmp_mean);
statsSummary.durMean    = run_paired_signrank(R.low_eventDur_mean,     R.high_eventDur_mean);
statsSummary.aucMean    = run_paired_signrank(R.low_eventAUC_mean,     R.high_eventAUC_mean);

disp('================ PAIRED SIGNRANK TESTS ================')
fprintf('Mean event rate:      n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.rateMean.n, statsSummary.rateMean.p, statsSummary.rateMean.deltaMedian);
fprintf('Median event rate:    n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.rateMedian.n, statsSummary.rateMedian.p, statsSummary.rateMedian.deltaMedian);
fprintf('Total event count:    n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.totalCount.n, statsSummary.totalCount.p, statsSummary.totalCount.deltaMedian);
fprintf('Frac active neurons:  n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.fracActive.n, statsSummary.fracActive.p, statsSummary.fracActive.deltaMedian);
fprintf('Mean event amplitude: n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.ampMean.n, statsSummary.ampMean.p, statsSummary.ampMean.deltaMedian);
fprintf('Mean event duration:  n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.durMean.n, statsSummary.durMean.p, statsSummary.durMean.deltaMedian);
fprintf('Mean event AUC:       n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.aucMean.n, statsSummary.aucMean.p, statsSummary.aucMean.deltaMedian);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_events_LL.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURE
% =========================

validRec = isfinite(R.low_eventRate_mean) & isfinite(R.high_eventRate_mean) & ...
           isfinite(R.low_totalEventCount) & isfinite(R.high_totalEventCount) & ...
           isfinite(R.low_fracActiveNeurons) & isfinite(R.high_fracActiveNeurons) & ...
           isfinite(R.low_eventAmp_mean) & isfinite(R.high_eventAmp_mean) & ...
           isfinite(R.low_eventDur_mean) & isfinite(R.high_eventDur_mean) & ...
           isfinite(R.low_eventAUC_mean) & isfinite(R.high_eventAUC_mean);

if any(validRec)

    xPair = [1 2];
    idxV  = find(validRec(:))';

    figE = figure('Name','Active vs inactive LL-event comparison', ...
        'Position',[100 100 1900 700]);
    tiledlayout(1,6, 'Padding','compact', 'TileSpacing','compact');

    % Mean rate
    nexttile; hold on;
    yLow = R.low_eventRate_mean(validRec);
    yHigh = R.high_eventRate_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_eventRate_mean(i), R.high_eventRate_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean event rate (events/s)');
    title(sprintf('Mean rate\nn=%d, p=%.3g %s', ...
        statsSummary.rateMean.n, statsSummary.rateMean.p, p_to_stars(statsSummary.rateMean.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.rateMean.p);

    % Total count
    nexttile; hold on;
    yLow = R.low_totalEventCount(validRec);
    yHigh = R.high_totalEventCount(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_totalEventCount(i), R.high_totalEventCount(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Total event count');
    title(sprintf('Total count\nn=%d, p=%.3g %s', ...
        statsSummary.totalCount.n, statsSummary.totalCount.p, p_to_stars(statsSummary.totalCount.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.totalCount.p);

    % Fraction active neurons
    nexttile; hold on;
    yLow = R.low_fracActiveNeurons(validRec);
    yHigh = R.high_fracActiveNeurons(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_fracActiveNeurons(i), R.high_fracActiveNeurons(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Fraction active neurons');
    title(sprintf('Frac active\nn=%d, p=%.3g %s', ...
        statsSummary.fracActive.n, statsSummary.fracActive.p, p_to_stars(statsSummary.fracActive.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.fracActive.p);

    % Mean event amplitude
    nexttile; hold on;
    yLow = R.low_eventAmp_mean(validRec);
    yHigh = R.high_eventAmp_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_eventAmp_mean(i), R.high_eventAmp_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean event amplitude (dF/F)');
    title(sprintf('Amplitude\nn=%d, p=%.3g %s', ...
        statsSummary.ampMean.n, statsSummary.ampMean.p, p_to_stars(statsSummary.ampMean.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.ampMean.p);

    % Mean event duration
    nexttile; hold on;
    yLow = R.low_eventDur_mean(validRec);
    yHigh = R.high_eventDur_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_eventDur_mean(i), R.high_eventDur_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean event duration (s)');
    title(sprintf('Duration\nn=%d, p=%.3g %s', ...
        statsSummary.durMean.n, statsSummary.durMean.p, p_to_stars(statsSummary.durMean.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.durMean.p);

    % Mean event AUC
    nexttile; hold on;
    yLow = R.low_eventAUC_mean(validRec);
    yHigh = R.high_eventAUC_mean(validRec);
    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);
    for i = idxV
        plot(xPair, [R.low_eventAUC_mean(i), R.high_eventAUC_mean(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end
    xlim([0.5 2.5]); set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Mean event AUC');
    title(sprintf('AUC\nn=%d, p=%.3g %s', ...
        statsSummary.aucMean.n, statsSummary.aucMean.p, p_to_stars(statsSummary.aucMean.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.aucMean.p);

    sgtitle('LL-event comparison: tail-defined inactive vs active windows');

    if saveOutputs
        saveas(figE, fullfile(outFolder, 'summary_events_active_vs_inactive_LL.png'));
    end
else
    warning('No valid recordings available for LL event summary figure.');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% ================== HELPERS ==================

function [eve_fra_cel, dff_dow_fra_cel, par] = ext_eve_LL(LL, fra_rat_cal, dff_cel_fra)
% dff_cel_fra: cells x frames
% outputs:
%   eve_fra_cel     : frames x cells
%   dff_dow_fra_cel : frames x cells

    par = [];
    par.fps = LL.targetFps;
    par.tauDecay = LL.tauDecay;
    par.BaselineNoiseMethod = 'Gaussian model';
    par.methodSignificatTransients = 'Dynamic threshold';
    par.confCutOff = LL.confCutOff;
    par.tracesArePercent = LL.tracesArePercent;
    par.plotFlag = LL.plotFlag;

    par.rasterQuickFixMergingFrames  = round(par.fps * LL.cri_dur_sec);
    par.rasterQuickFixDeletingFrames = round(par.fps * LL.cri_dur_sec);

    if isempty(dff_cel_fra)
        eve_fra_cel = [];
        dff_dow_fra_cel = [];
        return;
    end

    keep = ~all(isnan(dff_cel_fra), 2);
    dff_cel_fra = dff_cel_fra(keep, :);

    for i = 1:size(dff_cel_fra,1)
        xi = dff_cel_fra(i,:);
        if any(isnan(xi))
            xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
            if any(isnan(xi))
                xi = fillmissing(xi, 'constant', 0);
            end
        end
        dff_cel_fra(i,:) = xi;
    end

    freqRatio = par.fps / fra_rat_cal;
    [p, q] = rat(freqRatio);

    dff_dow_cel_fra = resample(dff_cel_fra', p, q)';   % cells x frames
    eve_fra_cel = identifySignifEvents_LL_fixed(dff_dow_cel_fra, par); % frames x cells
    eve_fra_cel = rasterQuickFix(eve_fra_cel, ...
        par.rasterQuickFixMergingFrames, ...
        par.rasterQuickFixDeletingFrames);

    dff_dow_fra_cel = dff_dow_cel_fra'; % frames x cells
end

function raster = identifySignifEvents_LL_fixed(traces, params)
% traces: cells x frames

    traces = double(traces);

    if isfield(params, 'tracesArePercent') && params.tracesArePercent
        deltaFoF = (traces ./ 100)';
    else
        deltaFoF = traces';
    end

    [deltaFoF, ~, sigma] = EstimateBaselineNoiseLL_fixed(deltaFoF);

    movements = 0;
    plotFlag = 0;
    if isfield(params, 'plotFlag')
        plotFlag = params.plotFlag;
    end

    [densityData, densityNoise, xev, yev] = NoiseModelLL(deltaFoF, sigma, movements, plotFlag);
    mapOfOdds = SignificantOdds(deltaFoF, sigma, movements, densityData, densityNoise, xev, params, plotFlag);
    [raster, ~] = Rasterize(deltaFoF, sigma, movements, mapOfOdds, xev, yev, params);
end

function [deltaFoF, mu, sigma] = EstimateBaselineNoiseLL_fixed(deltaFoF)
% deltaFoF: frames x cells

    numCells = size(deltaFoF,2);
    numFrames = size(deltaFoF,1);

    mu = nan(1, numCells);
    sigma = nan(1, numCells);

    for numcell = 1:numCells
        dataCell = deltaFoF(:,numcell);
        dataCell = dataCell(isfinite(dataCell));

        if numel(dataCell) < 10
            mu(numcell) = 0;
            sigma(numcell) = 1;
            continue;
        end

        [smoothDist,x] = ksdensity(dataCell);
        [~,indPeak] = max(smoothDist);

        xFit = x(1:indPeak);
        dataToFit = smoothDist(1:indPeak) / numFrames;

        try
            [sigmaFit, muFit, ~] = mygaussfit(xFit', dataToFit);
        catch
            sigmaFit = NaN;
            muFit = NaN;
        end

        if ~isreal(sigmaFit) || ~isfinite(sigmaFit) || sigmaFit <= 0
            dev = nanstd(dataCell);
            outliers = abs(dataCell) > 2 * dev;
            deltaF2 = dataCell;
            deltaF2(outliers) = NaN;
            sigmaFit = nanstd(deltaF2);
            muFit = nanmean(deltaF2);
        end

        if ~isfinite(sigmaFit) || sigmaFit <= 0
            sigmaFit = nanstd(dataCell);
        end
        if ~isfinite(sigmaFit) || sigmaFit <= 0
            sigmaFit = 1;
        end
        if ~isfinite(muFit)
            muFit = nanmean(dataCell);
        end
        if ~isfinite(muFit)
            muFit = 0;
        end

        sigma(numcell) = sigmaFit;
        mu(numcell) = muFit;
    end

    deltaFoF = bsxfun(@minus, deltaFoF, mu);
end

function [densityData, densityNoise, xev, yev] = NoiseModelLL(deltaFoF, sigma, movements, plotFlag)

    lambda=8;        

    deltaFoF(logical(movements),:)=NaN;
    transfDataMatrix=bsxfun(@rdivide,deltaFoF,sigma);
    points(:,1)=reshape(transfDataMatrix(1:end-1,:),[],1);
    points(:,2)=reshape(transfDataMatrix(2:end,:),[],1);
    points(isnan(points(:,1)),:)=[];
    points(isnan(points(:,2)),:)=[];

    pointsNeg=points(points(:,1)<0 & points(:,2)<0,1:2);
    Sigma=cov(pointsNeg(:,1), pointsNeg(:,2));
    Mu = [0 0];
    Sigma=[1 Sigma(1,2); Sigma(2,1) 1];
    dataGaussianCov = mvnrnd(Mu,Sigma,4*length(pointsNeg));
    clear pointsNeg

    mn=floor(min(min(min(transfDataMatrix)),min(reshape(dataGaussianCov,[],1))));
    mx=ceil(max(max(max(transfDataMatrix)),max(reshape(dataGaussianCov,[],1))));
    nevs=1000;
    [xev yev]=meshgrid(linspace(mn,mx,nevs),linspace(mn,mx,nevs));

    binColIdx=reshape(xev,[],1);
    binRowIdx=repmat(linspace(mn,mx,nevs),1,1000.)';

    [~,dist]=knnsearch([points(:,1) points(:,2)],[binColIdx binRowIdx],'K',100);
    binSpread=reshape(mean(dist,2),nevs,nevs);

    clear dist

    [~,distNoise]=knnsearch([dataGaussianCov(:,1) dataGaussianCov(:,2)],[binColIdx binRowIdx],'K',100);
    binSpreadNoise=reshape(mean(distNoise,2),nevs,nevs);

    clear binColIdx binRowIdx distNoise

    sizeWinFilt=8;
    smoothParam= 1/(sizeWinFilt*min(min(binSpread)));
    smoothParamNoise= 1/(sizeWinFilt*min(min(binSpreadNoise)));

    pointCol=interp1(xev(1,:),1:nevs,points(:,1),'nearest');
    pointRow=interp1(xev(1,:),1:nevs,points(:,2),'nearest');

    hist2dim=accumarray([pointRow pointCol],1,[nevs nevs]);

    pointColNoise=interp1(xev(1,:),1:nevs,dataGaussianCov(:,1) ,'nearest');
    pointRowNoise=interp1(xev(1,:),1:nevs,dataGaussianCov(:,2) ,'nearest');

    hist2dimNoise=accumarray([pointRowNoise pointColNoise],1,[nevs nevs]);

    densityData=zeros(size(xev));
    [binRows,binCols]=find(hist2dim);
    for i=1:length(binRows)
        sigmaFilt=smoothParam*binSpread(binRows(i),binCols(i));
        if mod(ceil(sizeWinFilt*sigmaFilt),2) == 0
            binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt)+1 ceil(sizeWinFilt*sigmaFilt)+1]),double(sigmaFilt));
        else
            binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt) ceil(sizeWinFilt*sigmaFilt)]),double(sigmaFilt));
        end

        widthRect=(size(binFilter,1)-1)/2;
        centerRect=widthRect+1;

        rectRows=binRows(i)-min(binRows(i)-1,widthRect):binRows(i)+min(size(densityData,1)-binRows(i),widthRect);
        rectCols=binCols(i)-min(binCols(i)-1,widthRect):binCols(i)+min(size(densityData,2)-binCols(i),widthRect);
        rectFiltRows=centerRect-min(binRows(i)-1,widthRect):centerRect+min(size(densityData,1)-binRows(i),widthRect);
        rectFiltCols=centerRect-min(binCols(i)-1,widthRect):centerRect+min(size(densityData,2)-binCols(i),widthRect);

        densityData(rectRows,rectCols)=densityData(rectRows,rectCols) + hist2dim(binRows(i),binCols(i))*binFilter(rectFiltRows,rectFiltCols);
    end

    densityNoise=zeros(size(xev));
    [binRowsNoise,binColsNoise]=find(hist2dimNoise);
    for i=1:length(binRowsNoise)
        sigmaFilt=smoothParamNoise*binSpreadNoise(binRowsNoise(i),binColsNoise(i));
        if mod(ceil(sizeWinFilt*sigmaFilt),2) == 0
            binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt)+1 ceil(sizeWinFilt*sigmaFilt)+1]),double(sigmaFilt));
        else
            binFilter=fspecial('gaussian',double([ceil(sizeWinFilt*sigmaFilt) ceil(sizeWinFilt*sigmaFilt)]),double(sigmaFilt));
        end

        widthRect=(size(binFilter,1)-1)/2;
        centerRect=widthRect+1;

        rectRows=binRowsNoise(i)-min(binRowsNoise(i)-1,widthRect):binRowsNoise(i)+min(size(densityNoise,1)-binRowsNoise(i),widthRect);
        rectCols=binColsNoise(i)-min(binColsNoise(i)-1,widthRect):binColsNoise(i)+min(size(densityNoise,2)-binColsNoise(i),widthRect);
        rectFiltRows=centerRect-min(binRowsNoise(i)-1,widthRect):centerRect+min(size(densityNoise,1)-binRowsNoise(i),widthRect);
        rectFiltCols=centerRect-min(binColsNoise(i)-1,widthRect):centerRect+min(size(densityNoise,2)-binColsNoise(i),widthRect);

        densityNoise(rectRows,rectCols)=densityNoise(rectRows,rectCols) + hist2dimNoise(binRowsNoise(i),binColsNoise(i))*binFilter(rectFiltRows,rectFiltCols);
    end

    densityData = Smooth1D(densityData,lambda);
    densityData = Smooth1D(densityData',lambda)';
    densityData = densityData./(sum(sum(densityData)));
    densityNoise = Smooth1D(densityNoise,lambda);
    densityNoise = Smooth1D(densityNoise',lambda)';
    densityNoise = densityNoise./(sum(sum(densityNoise)));

    if plotFlag
        figure
        z = log10(densityData);
        contour(xev,yev,z,linspace(max(max(z))-4,max(max(z)),20));
        hold on
        z = log10(densityNoise);
        contour(xev,yev,z,linspace(max(max(z))-4,max(max(z)),20));
        axis tight
        axis square
        ext=get(gca,'XLim');
        xlim([ext(1)-1 ext(2)+1]); ylim([ext(1)-1 ext(2)+1])
        hC=colorbar;
        tcks=get(hC,'Ytick');
        set(hC,'Ytick',unique(round(tcks)),'YTicklabel',10.^unique(round(tcks)));
        xlabel('Baseline noise normalized dFoF @ sample i');
        ylabel('Baseline noise normalized dFoF @ sample i+1');
        set(gcf,'color','w');
        set(gca,'FontSize',14);
    end
end

function Z = Smooth1D(Y,lambda)
    [m,~] = size(Y);
    E = eye(m);
    D1 = diff(E,1);
    D2 = diff(D1,1);
    P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
    Z = (E + P) \ Y;
end

function [mapOfOdds] = SignificantOdds(deltaFoF, sigma, movements, densityData, densityNoise, xev, params, plotFlag)

    deltaFoF(logical(movements),:)=NaN;
    transfDataMatrix=bsxfun(@rdivide,deltaFoF,sigma);
    points(:,1)=reshape(transfDataMatrix(1:end-1,:),[],1);
    points(:,2)=reshape(transfDataMatrix(2:end,:),[],1);
    points(isnan(points(:,1)),:)=[];
    points(isnan(points(:,2)),:)=[];

    pCutOff=(100-params.confCutOff)/100;
    mapOfOdds=densityNoise<=pCutOff*densityData;

    if plotFlag
        figure; plot(points(:,1),points(:,2),'k.'); axis equal; hold on
        h=imagesc(xev(1,:),xev(1,:), mapOfOdds); axis xy; axis tight
        alpha(h,0.6);
        xlabel('z-transformed value @ sample i');
        ylabel('z-transformed value @ sample i+1');
        set(gcf,'color','w');
        set(gca,'FontSize',14);
    end
end

function [raster, mapOfOddsJoint]=Rasterize(deltaFoF, sigma, movements, mapOfOdds, xev, yev, params)

    temp=deltaFoF;
    temp(logical(movements),:)=0;
    transfDataMatrix=bsxfun(@rdivide,temp,sigma);

    cc=bwconncomp(~mapOfOdds);
    stats=regionprops(cc,'PixelList');
    indZero=find(xev(1,:)>0,1,'first');

    for i=1:cc.NumObjects
        if ismember(indZero,stats(i).PixelList(:,1))
            break
        end
    end
    mapOfOddsCorrected=ones(size(mapOfOdds));
    mapOfOddsCorrected(cc.PixelIdxList{i})=0;

    noiseBias=1.5;
    factorDecay=exp(-1/(params.fps*params.tauDecay));
    decayMap=ones(size(mapOfOdds));
    rowsDecay=1:size(decayMap,1);
    for i=1:size(decayMap,2)
        decayMap(rowsDecay(yev(:,1)<factorDecay*(xev(1,i)-noiseBias)-noiseBias),i)=0;
    end

    riseMap=ones(size(mapOfOdds));
    riseMap(end,:)=0;

    mapOfOddsJoint = mapOfOddsCorrected & riseMap & decayMap;

    raster=zeros(size(transfDataMatrix,1),size(transfDataMatrix,2));
    for numNeuron=1:size(transfDataMatrix,2)
        [~,bins] = histc(transfDataMatrix(1:end,numNeuron),xev(1,:));
        for numFrame=3:size(transfDataMatrix,1)-2
            optA= (mapOfOddsJoint(bins(numFrame+1),bins(numFrame)) & mapOfOddsJoint(bins(numFrame),bins(numFrame-1)));
            optB= (mapOfOddsJoint(bins(numFrame),bins(numFrame-1)) & mapOfOddsJoint(bins(numFrame-1),bins(numFrame-2)));
            optC= (mapOfOddsJoint(bins(numFrame+2),bins(numFrame+1)) & mapOfOddsJoint(bins(numFrame+1),bins(numFrame)));
            if optA || optB || optC
                raster(numFrame,numNeuron)=1;
            end
        end
    end
end

function rasterFixed = rasterQuickFix(raster, N, M)
% raster: frames x cells

    rasterFixed=zeros(size(raster,1), size(raster, 2));

    for i=1:size(raster, 2)
        eventsOneCell=raster(:, i);
        eventsFixed=eventsOneCell;

        indexesEvents=find(eventsOneCell);
        gapBetwEvents=diff(indexesEvents);
        for k=1:N
            toFix=find(gapBetwEvents==1+k);
            for p=1:length(toFix)
                eventsFixed(indexesEvents(toFix(p))+1:indexesEvents(toFix(p))+k)=1;
            end
        end

        indexesNoEvent=find(not(eventsFixed));
        gapBetwNoEvents=diff(indexesNoEvent);
        for k=1:M
            toFix=find(gapBetwNoEvents==1+k);
            for p=1:length(toFix)
                eventsFixed(indexesNoEvent(toFix(p))+1:indexesNoEvent(toFix(p))+k)=0;
            end
        end

        rasterFixed(:, i)=eventsFixed;
    end
end

function [sigma,mu,A]=mygaussfit(x,y,h,optionStab)

    if nargin==2, h=0.2; optionStab=1; end

    ymax=max(y);
    xnew=[];
    ynew=[];
    for n=1:length(x)
        if y(n)>ymax*h
            xnew=[xnew,x(n)]; %#ok<AGROW>
            ynew=[ynew,y(n)]; %#ok<AGROW>
        end
    end

    ylog=log(ynew);
    xlog=xnew;

    if optionStab==1
        p=polyfit(xlog,ylog,2);
        A2=p(1);
        A1=p(2);
        A0=p(3);
        sigma=sqrt(-1/(2*A2));
        mu=A1*sigma^2;
        A=exp(A0+mu^2/(2*sigma^2));
    end

    if optionStab==2
        [p,~,s2]=polyfit(xlog,ylog,2);
        A2=p(1);
        A1=p(2);
        A0=p(3);
        A1=A1/s2(2)-2*A2*s2(1)/s2(2);
        A2=A2/s2(2)^2;
        sigma=sqrt(-1/(2*A2));
        mu=A1*sigma^2;
        A=exp(A0+mu^2/(2*sigma^2));
    end
end

function [M, ampVals, durVals, aucVals] = summarize_LL_events_in_window(eve_fra_cel, dff_fra_cel, idx, fps)
% eve_fra_cel : frames x cells logical
% dff_fra_cel : frames x cells
% idx         : frame indices for one window

    idx = idx(:);
    idx = idx(idx >= 1 & idx <= size(eve_fra_cel,1));

    evWin  = logical(eve_fra_cel(idx,:));   % Tw x N
    dffWin = dff_fra_cel(idx,:);            % Tw x N
    durSecWin = numel(idx) / fps;

    N = size(evWin,2);
    counts = zeros(N,1);
    rates  = zeros(N,1);
    isActive = false(N,1);
    meanAmpPerNeuron = nan(N,1);
    meanDurPerNeuron = nan(N,1);
    meanAUCPerNeuron = nan(N,1);

    ampVals = [];
    durVals = [];
    aucVals = [];

    totalCount = 0;

    for n = 1:N
        r = evWin(:,n);
        x = dffWin(:,n);

        d = diff([false; r; false]);
        starts = find(d == 1);
        ends   = find(d == -1) - 1;

        nEvents = numel(starts);
        counts(n) = nEvents;
        rates(n)  = nEvents / durSecWin;
        isActive(n) = nEvents > 0;
        totalCount = totalCount + nEvents;

        if nEvents == 0
            continue;
        end

        ampN = nan(nEvents,1);
        durN = nan(nEvents,1);
        aucN = nan(nEvents,1);

        for k = 1:nEvents
            segIdx = starts(k):ends(k);
            seg = x(segIdx);

            durN(k) = numel(segIdx) / fps;
            ampN(k) = max(seg, [], 'omitnan');
            aucN(k) = trapz(seg) / fps;
        end

        meanAmpPerNeuron(n) = mean(ampN, 'omitnan');
        meanDurPerNeuron(n) = mean(durN, 'omitnan');
        meanAUCPerNeuron(n) = mean(aucN, 'omitnan');

        ampVals = [ampVals; ampN]; %#ok<AGROW>
        durVals = [durVals; durN]; %#ok<AGROW>
        aucVals = [aucVals; aucN]; %#ok<AGROW>
    end

    M = struct();
    M.counts_per_neuron   = counts;
    M.rates_per_neuron    = rates;
    M.rate_mean           = mean(rates, 'omitnan');
    M.rate_median         = median(rates, 'omitnan');
    M.total_count         = totalCount;
    M.frac_active_neurons = mean(isActive, 'omitnan');
    M.amp_mean            = mean(ampVals, 'omitnan');
    M.dur_mean            = mean(durVals, 'omitnan');
    M.auc_mean            = mean(aucVals, 'omitnan');
    M.amp_mean_per_neuron = meanAmpPerNeuron;
    M.dur_mean_per_neuron = meanDurPerNeuron;
    M.auc_mean_per_neuron = meanAUCPerNeuron;
end

function [dff, fps] = extract_dff_and_fps(S)
% Try to robustly extract dff matrix and fps from loaded dff mat file

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

function fig = make_LL_window_raster_QC(eve_fra_cel, dff_fra_cel, low_idx, high_idx, fps, recLabel)
% eve_fra_cel : frames x cells logical
% dff_fra_cel : frames x cells
% low_idx     : frame indices for least-active window
% high_idx    : frame indices for most-active window

    low_idx  = low_idx(:);
    high_idx = high_idx(:);

    low_idx  = low_idx(low_idx >= 1 & low_idx <= size(eve_fra_cel,1));
    high_idx = high_idx(high_idx >= 1 & high_idx <= size(eve_fra_cel,1));

    eveLow = logical(eve_fra_cel(low_idx, :));     % Tw x N
    eveHigh = logical(eve_fra_cel(high_idx, :));   % Tw x N

    dffLow = dff_fra_cel(low_idx, :);              % Tw x N
    dffHigh = dff_fra_cel(high_idx, :);            % Tw x N

    % Sort neurons by event occupancy within each window
    lowSortMetric  = sum(eveLow, 1, 'omitnan');
    highSortMetric = sum(eveHigh, 1, 'omitnan');

    [~, ordLow]  = sort(lowSortMetric,  'ascend');
    [~, ordHigh] = sort(highSortMetric, 'ascend');

    eveLowS  = eveLow(:, ordLow)';
    eveHighS = eveHigh(:, ordHigh)';

    dffLowS  = dffLow(:, ordLow)';
    dffHighS = dffHigh(:, ordHigh)';

    % Time vectors
    tLow  = (0:numel(low_idx)-1) / fps;
    tHigh = (0:numel(high_idx)-1) / fps;

    % Bottom summary traces:
    % proportion of neurons active at each frame
    popLow  = mean(eveLow,  2, 'omitnan');
    popHigh = mean(eveHigh, 2, 'omitnan');

    % Optional: mean dF/F across neurons
    meanDffLow  = mean(dffLow,  2, 'omitnan');
    meanDffHigh = mean(dffHigh, 2, 'omitnan');

    fig = figure('Name', ['QC raster ' recLabel], ...
        'Position', [100 80 1200 850], ...
        'Color', 'w');

    tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

    % ---------------- LEFT TOP: low raster
    nexttile; hold on;
    imagesc(tLow, 1:size(eveLowS,1), eveLowS);
    colormap(gca, flipud(gray));
    caxis([0 1]);
    set(gca, 'YDir','normal');
    xlabel('Time (s)');
    ylabel('Neurons');
    title(sprintf('Least active window\n%d neurons, %d events', ...
        size(eveLowS,1), sum(eveLowS(:))));
    box off;

    % scale bar-ish y label hint
    text(0.01, 0.98, sprintf('Window = %.1f min', numel(low_idx)/fps/60), ...
        'Units','normalized', 'HorizontalAlignment','left', ...
        'VerticalAlignment','top', 'FontSize',10, 'BackgroundColor','w');

    % ---------------- RIGHT TOP: high raster
    nexttile; hold on;
    imagesc(tHigh, 1:size(eveHighS,1), eveHighS);
    colormap(gca, flipud(gray));
    caxis([0 1]);
    set(gca, 'YDir','normal');
    xlabel('Time (s)');
    ylabel('Neurons');
    title(sprintf('Most active window\n%d neurons, %d events', ...
        size(eveHighS,1), sum(eveHighS(:))));
    box off;

    text(0.01, 0.98, sprintf('Window = %.1f min', numel(high_idx)/fps/60), ...
        'Units','normalized', 'HorizontalAlignment','left', ...
        'VerticalAlignment','top', 'FontSize',10, 'BackgroundColor','w');

    % ---------------- LEFT BOTTOM: low summary
    nexttile; hold on;
    yyaxis left
    plot(tLow, popLow, 'k-', 'LineWidth', 1.2);
    ylabel('Active neuron fraction');

    yyaxis right
    plot(tLow, meanDffLow, '-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5]);
    ylabel('Mean dF/F');

    xlabel('Time (s)');
    title('Least active window summary');
    box off;

    % ---------------- RIGHT BOTTOM: high summary
    nexttile; hold on;
    yyaxis left
    plot(tHigh, popHigh, 'k-', 'LineWidth', 1.2);
    ylabel('Active neuron fraction');

    yyaxis right
    plot(tHigh, meanDffHigh, '-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5]);
    ylabel('Mean dF/F');

    xlabel('Time (s)');
    title('Most active window summary');
    box off;

    sgtitle(sprintf('LL event raster QC: %s', recLabel), 'FontWeight', 'bold');
end