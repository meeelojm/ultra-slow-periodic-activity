%% =========================================================================
% bout_onset_phase_reset_ultraslow_from_tail_PER_FISH.m
%
% SELF-CONTAINED SCRIPT
%
% What it does:
%   1) Recursively finds recordings that contain:
%        - tail_quick.mat
%        - dffs_repact_respcells.mat
%   2) Loads neural dff and fps
%   3) Extracts tail bout onsets
%   4) Bandpass-filters each neuron in an ultraslow band
%   5) Uses Hilbert transform to get phase and instantaneous frequency
%   6) Aligns neurons to bout onsets
%   7) Quantifies:
%        - ITPC (phase concentration across bouts)
%        - pre vs post ITPC
%        - pre vs post instantaneous frequency
%        - reset strength = max post-bout ITPC minus pre-bout baseline
%   8) Summarizes across neurons within each recording
%   9) Builds summary table, paired signrank stats, and summary figures
%  10) Builds pooled PSTH-like figures for:
%        - ITPC
%        - bandpassed neural activity
%  11) NEW: pools all valid neurons PER FISH across recordings and creates:
%        - per-fish ITPC PSTH-like figures
%        - per-fish activity PSTH-like figures
%        - per-fish neuron-level box/scatter statistical figures
%
% IMPORTANT:
%   This script tests bout-triggered phase/frequency reset of ultraslow
%   neural activity. It is not using least-vs-most-active windows.
%
% Toolbox requirements:
%   - Signal Processing Toolbox (butter, filtfilt, hilbert)
% =========================================================================

clear; clc; close all;

%% =========================
% 0) USER SETTINGS
% =========================

% Root folder
rootFolder = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low';

% File names
tailFileName = 'tail_quick.mat';
dffFileName  = 'dffs_repact_respcells.mat';

% Ultraslow band of interest
frequencyRange = [0.01 0.05];   % Hz

% Bout-triggered analysis windows
periSec     = 120;              % total peri-bout half window, seconds
preStatSec  = [-40 -5];         % pre-bout statistics window
postStatSec = [5 40];           % post-bout statistics window
postPeakSec = [0 40];           % window in which to look for post-bout ITPC peak

% Bout inclusion
minBoutCount    = 8;            % minimum usable bouts per recording
minInterBoutSec = 20;           % reject bouts too close together

% Neural inclusion
maxNanFracPerNeuron = 0.20;     % max NaN fraction in whole recording

% Filtering
butterOrder = 3;

% Optional per-recording figures
makePerRecordingITPCFigure     = true;
makePerRecordingActivityFigure = true;
makePerRecordingStatsFigure    = true;
makeRecordingSummaryFigure = false;   % old experiment-level summary across recordings

% NEW: per-fish figures
makePerFishITPCFigure     = true;
makePerFishActivityFigure = true;
makePerFishStatsFigure    = true;

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'bout_reset_ultraslow_results_per_fish');

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

dffMap = containers.Map;
for i = 1:numel(dffFiles)
    dffMap(dffFiles(i).folder) = fullfile(dffFiles(i).folder, dffFiles(i).name);
end

recordings = struct('folder', {}, 'tailPath', {}, 'dffPath', {}, 'fishID', {});
for i = 1:numel(tailFiles)
    thisFolder = tailFiles(i).folder;
    if isKey(dffMap, thisFolder)
        recordings(end+1).folder = thisFolder; %#ok<SAGROW>
        recordings(end).tailPath = fullfile(tailFiles(i).folder, tailFiles(i).name);
        recordings(end).dffPath  = dffMap(thisFolder);
        recordings(end).fishID   = extract_fish_id(thisFolder);
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
R.recordingLabel           = cell(nR,1);
R.folder                   = cell(nR,1);
R.fishID                   = cell(nR,1);

R.nNeurons_total           = nan(nR,1);
R.nNeurons_used            = nan(nR,1);

R.fps                      = nan(nR,1);
R.nFrames                  = nan(nR,1);
R.durationMin              = nan(nR,1);

R.nBouts_total             = nan(nR,1);
R.nBouts_used              = nan(nR,1);

R.pre_ITPC_mean            = nan(nR,1);
R.post_ITPC_mean           = nan(nR,1);
R.delta_ITPC_mean          = nan(nR,1);
R.resetStrength_mean       = nan(nR,1);

R.pre_instFreq_mean        = nan(nR,1);
R.post_instFreq_mean       = nan(nR,1);
R.delta_instFreq_mean      = nan(nR,1);

R.pre_ITPC_median          = nan(nR,1);
R.post_ITPC_median         = nan(nR,1);
R.delta_ITPC_median        = nan(nR,1);
R.resetStrength_median     = nan(nR,1);

R.pre_instFreq_median      = nan(nR,1);
R.post_instFreq_median     = nan(nR,1);
R.delta_instFreq_median    = nan(nR,1);

R.sessionLabel             = cell(nR,1);

all_ITPC_curves              = cell(nR,1);   % mean ITPC curve per recording
all_freq_curves              = cell(nR,1);   % mean inst freq curve per recording
all_meanFilt_curves          = cell(nR,1);   % mean bandpassed activity curve per recording
all_tSec                     = cell(nR,1);   % time axis per recording

all_ITPC_curves_perNeuron    = cell(nR,1);   % full ITPC matrix: neurons x time
all_freq_curves_perNeuron    = cell(nR,1);   % full freq matrix: neurons x time
all_filtCurves_perNeuron     = cell(nR,1);   % full filtered activity matrix: neurons x time

all_neuron_preITPC           = cell(nR,1);
all_neuron_postITPC          = cell(nR,1);
all_neuron_deltaITPC         = cell(nR,1);
all_neuron_resetStrength     = cell(nR,1);

all_neuron_preFreq           = cell(nR,1);
all_neuron_postFreq          = cell(nR,1);
all_neuron_deltaFreq         = cell(nR,1);

%% =========================
% 3) PROCESS EACH RECORDING
% =========================

for r = 1:nR

    fprintf('\n============================================================\n');
    fprintf('Recording %d / %d\n', r, nR);
    fprintf('%s\n', recordings(r).folder);

    R.folder{r}          = recordings(r).folder;
    R.recordingLabel{r}  = make_recording_label(recordings(r).folder);
    R.fishID{r}          = recordings(r).fishID;
    R.sessionLabel{r} = extract_session_label(recordings(r).folder);

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

    % ---------------------------------------------------------------------
    % 3.2 Load tail and extract bout onsets
    % ---------------------------------------------------------------------
    S_tail = load(recordings(r).tailPath);
    [boutOnsetsFrames, ~] = extract_bout_onsets_from_tail(S_tail, nFrames, fps);

    if isempty(boutOnsetsFrames)
        warning('No bout onsets could be extracted for %s. Skipping.', recordings(r).tailPath);
        continue;
    end

    boutOnsetsFrames = unique(round(boutOnsetsFrames(:)));
    boutOnsetsFrames = boutOnsetsFrames(isfinite(boutOnsetsFrames));
    boutOnsetsFrames = boutOnsetsFrames(boutOnsetsFrames >= 1 & boutOnsetsFrames <= nFrames);

    R.nBouts_total(r) = numel(boutOnsetsFrames);

    fprintf('  Extracted %d raw bout onsets.\n', numel(boutOnsetsFrames));

    % ---------------------------------------------------------------------
    % 3.3 Keep only bouts far from edges and from each other
    % ---------------------------------------------------------------------
    periFrames = round(periSec * fps);
    minInterBoutFrames = round(minInterBoutSec * fps);

    keepBout = boutOnsetsFrames > periFrames & boutOnsetsFrames < (nFrames - periFrames);
    boutOnsetsFrames = boutOnsetsFrames(keepBout);

    if isempty(boutOnsetsFrames)
        warning('No bouts remain after edge exclusion in %s. Skipping.', recordings(r).folder);
        continue;
    end

    boutKeep = true(size(boutOnsetsFrames));
    for b = 2:numel(boutOnsetsFrames)
        if boutOnsetsFrames(b) - boutOnsetsFrames(b-1) < minInterBoutFrames
            boutKeep(b) = false;
        end
    end
    boutOnsetsFrames = boutOnsetsFrames(boutKeep);

    R.nBouts_used(r) = numel(boutOnsetsFrames);

    fprintf('  Using %d bouts after edge/spacing exclusion.\n', numel(boutOnsetsFrames));

    if numel(boutOnsetsFrames) < minBoutCount
        warning('Too few usable bouts in %s. Skipping.', recordings(r).folder);
        continue;
    end

    % ---------------------------------------------------------------------
    % 3.4 Keep usable neurons and fill NaNs
    % ---------------------------------------------------------------------
    keepNeuron = mean(isnan(dff), 2) <= maxNanFracPerNeuron;
    dff_use = dff(keepNeuron, :);
    dff_use = dff_use(~all(isnan(dff_use),2), :);

    R.nNeurons_used(r) = size(dff_use,1);

    if size(dff_use,1) < 3
        warning('Too few usable neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    dff_use = fill_nans_rowwise(dff_use);

    % ---------------------------------------------------------------------
    % 3.5 Bout-triggered reset metrics
    % ---------------------------------------------------------------------
    try
        M = compute_bout_reset_metrics( ...
            dff_use, boutOnsetsFrames, fps, ...
            frequencyRange, butterOrder, periSec, ...
            preStatSec, postStatSec, postPeakSec);
    catch ME
        warning('Bout reset analysis failed in %s\n%s', recordings(r).folder, ME.message);
        continue;
    end

    requiredFields = { ...
        'tSec', 'itpcCurves', 'freqCurves', 'meanFiltCurves', ...
        'meanITPC', 'meanInstFreq', 'meanPopulationFilt', ...
        'neuron_preITPC', 'neuron_postITPC', 'neuron_deltaITPC', 'neuron_resetStrength', ...
        'neuron_preFreq', 'neuron_postFreq', 'neuron_deltaFreq'};

    missingFields = requiredFields(~isfield(M, requiredFields));
    if ~isempty(missingFields)
        warning('compute_bout_reset_metrics is missing fields in %s: %s', ...
            recordings(r).folder, strjoin(missingFields, ', '));
        continue;
    end

    % Store recording-level summaries
    R.pre_ITPC_mean(r)         = mean(M.neuron_preITPC, 'omitnan');
    R.post_ITPC_mean(r)        = mean(M.neuron_postITPC, 'omitnan');
    R.delta_ITPC_mean(r)       = mean(M.neuron_deltaITPC, 'omitnan');
    R.resetStrength_mean(r)    = mean(M.neuron_resetStrength, 'omitnan');

    R.pre_instFreq_mean(r)     = mean(M.neuron_preFreq, 'omitnan');
    R.post_instFreq_mean(r)    = mean(M.neuron_postFreq, 'omitnan');
    R.delta_instFreq_mean(r)   = mean(M.neuron_deltaFreq, 'omitnan');

    R.pre_ITPC_median(r)       = median(M.neuron_preITPC, 'omitnan');
    R.post_ITPC_median(r)      = median(M.neuron_postITPC, 'omitnan');
    R.delta_ITPC_median(r)     = median(M.neuron_deltaITPC, 'omitnan');
    R.resetStrength_median(r)  = median(M.neuron_resetStrength, 'omitnan');

    R.pre_instFreq_median(r)   = median(M.neuron_preFreq, 'omitnan');
    R.post_instFreq_median(r)  = median(M.neuron_postFreq, 'omitnan');
    R.delta_instFreq_median(r) = median(M.neuron_deltaFreq, 'omitnan');

    % Store curves
    all_ITPC_curves{r}           = M.meanITPC(:);
    all_freq_curves{r}           = M.meanInstFreq(:);
    all_meanFilt_curves{r}       = M.meanPopulationFilt(:);
    all_tSec{r}                  = M.tSec(:);

    all_ITPC_curves_perNeuron{r} = M.itpcCurves;
    all_freq_curves_perNeuron{r} = M.freqCurves;
    all_filtCurves_perNeuron{r}  = M.meanFiltCurves;

    % Store per-neuron summaries
    all_neuron_preITPC{r}        = M.neuron_preITPC(:);
    all_neuron_postITPC{r}       = M.neuron_postITPC(:);
    all_neuron_deltaITPC{r}      = M.neuron_deltaITPC(:);
    all_neuron_resetStrength{r}  = M.neuron_resetStrength(:);

    all_neuron_preFreq{r}        = M.neuron_preFreq(:);
    all_neuron_postFreq{r}       = M.neuron_postFreq(:);
    all_neuron_deltaFreq{r}      = M.neuron_deltaFreq(:);

    fprintf('  Mean ITPC pre/post: %.4f -> %.4f | delta = %.4f\n', ...
        R.pre_ITPC_mean(r), R.post_ITPC_mean(r), R.delta_ITPC_mean(r));
    fprintf('  Median ITPC pre/post: %.4f -> %.4f | delta = %.4f\n', ...
        R.pre_ITPC_median(r), R.post_ITPC_median(r), R.delta_ITPC_median(r));
    fprintf('  Mean inst. freq pre/post: %.4f -> %.4f Hz | delta = %.6f Hz\n', ...
        R.pre_instFreq_mean(r), R.post_instFreq_mean(r), R.delta_instFreq_mean(r));
    fprintf('  Median inst. freq pre/post: %.4f -> %.4f Hz | delta = %.6f Hz\n', ...
        R.pre_instFreq_median(r), R.post_instFreq_median(r), R.delta_instFreq_median(r));
    fprintf('  Mean reset strength: %.4f | Median reset strength: %.4f\n', ...
        R.resetStrength_mean(r), R.resetStrength_median(r));

    % Optional per-recording ITPC figure
    if makePerRecordingITPCFigure
        try
            figTmp = plot_itpc_psth_like( ...
                M.itpcCurves, M.tSec, M.neuron_resetStrength, R.recordingLabel{r});
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('itpc_psth_like_%02d_%s.png', r, sanitize_filename(R.recordingLabel{r}))));
            end
        catch ME
            warning('Could not make ITPC PSTH-like figure for %s\n%s', ...
                recordings(r).folder, ME.message);
        end
    end

    % Optional per-recording activity figure
    if makePerRecordingActivityFigure
        try
            figTmp = plot_activity_psth_like( ...
                M.meanFiltCurves, M.tSec, M.neuron_resetStrength, ...
                R.recordingLabel{r}, 'Bandpassed \DeltaF/F');
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('activity_psth_like_%02d_%s.png', r, sanitize_filename(R.recordingLabel{r}))));
            end
        catch ME
            warning('Could not make activity PSTH-like figure for %s\n%s', ...
                recordings(r).folder, ME.message);
        end
    end

        % Optional per-recording ITPC figure
    if makePerRecordingITPCFigure
        try
            figTmp = plot_itpc_psth_like( ...
                M.itpcCurves, M.tSec, M.neuron_resetStrength, R.recordingLabel{r});
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('itpc_psth_like_%02d_%s.png', r, sanitize_filename(R.recordingLabel{r}))));
            end
        catch ME
            warning('Could not make ITPC PSTH-like figure for %s\n%s', ...
                recordings(r).folder, ME.message);
        end
    end

    % Optional per-recording activity figure
    if makePerRecordingActivityFigure
        try
            figTmp = plot_activity_psth_like( ...
                M.meanFiltCurves, M.tSec, M.neuron_resetStrength, ...
                R.recordingLabel{r}, 'Bandpassed \DeltaF/F');
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('activity_psth_like_%02d_%s.png', r, sanitize_filename(R.recordingLabel{r}))));
            end
        catch ME
            warning('Could not make activity PSTH-like figure for %s\n%s', ...
                recordings(r).folder, ME.message);
        end
    end

    % Optional per-recording neuron-stat figure
    if makePerRecordingStatsFigure
        try
            figTmp = plot_per_entity_stats( ...
                M.neuron_preITPC, M.neuron_postITPC, M.neuron_deltaITPC, M.neuron_resetStrength, ...
                M.neuron_preFreq, M.neuron_postFreq, M.neuron_deltaFreq, ...
                R.recordingLabel{r});
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('stats_neurons_%02d_%s.png', r, sanitize_filename(R.recordingLabel{r}))));
            end
        catch ME
            warning('Could not make stats figure for %s\n%s', ...
                recordings(r).folder, ME.message);
        end
    end
end

%% =========================
% 4) BUILD RECORDING SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, R.fishID, ...
    R.nNeurons_total, R.nNeurons_used, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.nBouts_total, R.nBouts_used, ...
    R.pre_ITPC_mean, R.post_ITPC_mean, R.delta_ITPC_mean, R.resetStrength_mean, ...
    R.pre_instFreq_mean, R.post_instFreq_mean, R.delta_instFreq_mean, ...
    R.pre_ITPC_median, R.post_ITPC_median, R.delta_ITPC_median, R.resetStrength_median, ...
    R.pre_instFreq_median, R.post_instFreq_median, R.delta_instFreq_median, ...
    'VariableNames', { ...
    'recordingLabel','folder','fishID', ...
    'nNeurons_total','nNeurons_used', ...
    'fps','nFrames','durationMin', ...
    'nBouts_total','nBouts_used', ...
    'pre_ITPC_mean','post_ITPC_mean','delta_ITPC_mean','resetStrength_mean', ...
    'pre_instFreq_mean','post_instFreq_mean','delta_instFreq_mean', ...
    'pre_ITPC_median','post_ITPC_median','delta_ITPC_median','resetStrength_median', ...
    'pre_instFreq_median','post_instFreq_median','delta_instFreq_median'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_bout_reset_ultraslow_recordings.csv'));
end

%% =========================
% 5) RECORDING-LEVEL PAIRED STATISTICS
% =========================

statsSummary = struct();

statsSummary.ITPC_pre_post      = run_paired_signrank(R.pre_ITPC_median, R.post_ITPC_median);
statsSummary.instFreq_pre_post  = run_paired_signrank(R.pre_instFreq_median, R.post_instFreq_median);

statsSummary.deltaITPC_vs_zero      = run_signrank_vs_zero(R.delta_ITPC_median);
statsSummary.resetStrength_vs_zero  = run_signrank_vs_zero(R.resetStrength_median);
statsSummary.deltaFreq_vs_zero      = run_signrank_vs_zero(R.delta_instFreq_median);

disp('================ RECORDING-LEVEL BOUT RESET STATISTICS ================')
fprintf('ITPC pre vs post:            n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.ITPC_pre_post.n, statsSummary.ITPC_pre_post.p, statsSummary.ITPC_pre_post.deltaMedian);
fprintf('InstFreq pre vs post:        n = %d, p = %.4g, median Δ = %.4g Hz\n', ...
    statsSummary.instFreq_pre_post.n, statsSummary.instFreq_pre_post.p, statsSummary.instFreq_pre_post.deltaMedian);
fprintf('Delta ITPC vs 0:             n = %d, p = %.4g, median = %.4g\n', ...
    statsSummary.deltaITPC_vs_zero.n, statsSummary.deltaITPC_vs_zero.p, statsSummary.deltaITPC_vs_zero.medianValue);
fprintf('Reset strength vs 0:         n = %d, p = %.4g, median = %.4g\n', ...
    statsSummary.resetStrength_vs_zero.n, statsSummary.resetStrength_vs_zero.p, statsSummary.resetStrength_vs_zero.medianValue);
fprintf('Delta instFreq vs 0:         n = %d, p = %.4g, median = %.4g Hz\n', ...
    statsSummary.deltaFreq_vs_zero.n, statsSummary.deltaFreq_vs_zero.p, statsSummary.deltaFreq_vs_zero.medianValue);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_bout_reset_ultraslow_recordings.mat'), 'statsSummary');
end

%% =========================
% 6) RECORDING-LEVEL SUMMARY FIGURE
% =========================

validRec = isfinite(R.pre_ITPC_median) & isfinite(R.post_ITPC_median) & ...
           isfinite(R.pre_instFreq_median) & isfinite(R.post_instFreq_median);

if any(validRec)

    figR = figure('Name','Bout-triggered ultraslow reset | recording summary', ...
        'Position',[100 100 1400 700]);
    tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact');

    xPair = [1 2];
    idxV  = find(validRec(:))';

    % ITPC pre vs post
    nexttile; hold on;
    yLow = R.pre_ITPC_median(validRec);
    yHigh = R.post_ITPC_median(validRec);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxV
        plot(xPair, [R.pre_ITPC_median(i), R.post_ITPC_median(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Pre-bout','Post-bout'});
    ylabel('Median ITPC');
    title(sprintf('Phase concentration\nn=%d, p=%.3g %s', ...
        statsSummary.ITPC_pre_post.n, statsSummary.ITPC_pre_post.p, p_to_stars(statsSummary.ITPC_pre_post.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.ITPC_pre_post.p);

    % Inst frequency pre vs post
    nexttile; hold on;
    yLow = R.pre_instFreq_median(validRec);
    yHigh = R.post_instFreq_median(validRec);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxV
        plot(xPair, [R.pre_instFreq_median(i), R.post_instFreq_median(i)], '-o', ...
            'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Pre-bout','Post-bout'});
    ylabel('Median instantaneous frequency (Hz)');
    title(sprintf('Instantaneous frequency\nn=%d, p=%.3g %s', ...
        statsSummary.instFreq_pre_post.n, statsSummary.instFreq_pre_post.p, p_to_stars(statsSummary.instFreq_pre_post.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.instFreq_pre_post.p);

    % Mean ITPC timecourse
    nexttile; hold on;
    validCurves = ~cellfun(@isempty, all_ITPC_curves) & ~cellfun(@isempty, all_tSec);
    if any(validCurves)
        refIdx = find(validCurves, 1, 'first');
        tSec = all_tSec{refIdx}(:);

        keepIdx = [];
        for i = find(validCurves(:))'
            if numel(all_ITPC_curves{i}) == numel(tSec)
                keepIdx(end+1) = i; %#ok<SAGROW>
            end
        end

        itpcMat = nan(numel(tSec), numel(keepIdx));
        for k = 1:numel(keepIdx)
            itpcMat(:,k) = all_ITPC_curves{keepIdx(k)}(:);
        end

        m = mean(itpcMat, 2, 'omitnan');
        s = std(itpcMat, 0, 2, 'omitnan') ./ sqrt(size(itpcMat,2));

        fill([tSec; flipud(tSec)], [m-s; flipud(m+s)], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(tSec, m, 'b', 'LineWidth', 1.8);

        xline(0, '--k');
        xline(preStatSec(1), ':k');
        xline(preStatSec(2), ':k');
        xline(postStatSec(1), ':r');
        xline(postStatSec(2), ':r');

        xlabel('Time from bout onset (s)');
        ylabel('ITPC');
        title('Mean ITPC timecourse');
        xlim([-periSec periSec]);
    else
        text(0.1,0.5,'No valid ITPC curves','FontSize',12);
        axis off;
    end

    % Mean instantaneous frequency timecourse
    nexttile; hold on;
    validCurves = ~cellfun(@isempty, all_freq_curves) & ~cellfun(@isempty, all_tSec);
    if any(validCurves)
        refIdx = find(validCurves, 1, 'first');
        tSec = all_tSec{refIdx}(:);

        keepIdx = [];
        for i = find(validCurves(:))'
            if numel(all_freq_curves{i}) == numel(tSec)
                keepIdx(end+1) = i; %#ok<SAGROW>
            end
        end

        fMat = nan(numel(tSec), numel(keepIdx));
        for k = 1:numel(keepIdx)
            fMat(:,k) = all_freq_curves{keepIdx(k)}(:);
        end

        m = mean(fMat, 2, 'omitnan');
        s = std(fMat, 0, 2, 'omitnan') ./ sqrt(size(fMat,2));

        fill([tSec; flipud(tSec)], [m-s; flipud(m+s)], ...
            'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(tSec, m, 'r', 'LineWidth', 1.8);

        xline(0, '--k');
        xline(preStatSec(1), ':k');
        xline(preStatSec(2), ':k');
        xline(postStatSec(1), ':r');
        xline(postStatSec(2), ':r');

        xlabel('Time from bout onset (s)');
        ylabel('Instantaneous frequency (Hz)');
        title('Mean instantaneous frequency');
        xlim([-periSec periSec]);
    else
        text(0.1,0.5,'No valid frequency curves','FontSize',12);
        axis off;
    end

    sgtitle('Bout-triggered reset of ultraslow neural activity | recordings');

    if saveOutputs
        saveas(figR, fullfile(outFolder, 'summary_bout_reset_ultraslow_recordings.png'));
    end
else
    warning('No valid recordings available for recording summary figure.');
end

%% =========================
% 7) POOLED PSTH-LIKE ITPC FIGURE (ALL NEURONS)
% =========================

if any(~cellfun(@isempty, all_tSec))
    refIdx  = find(~cellfun(@isempty, all_tSec), 1, 'first');
    tCommon = all_tSec{refIdx}(:)';

    allITPC  = [];
    allReset = [];

    for r = 1:numel(all_ITPC_curves_perNeuron)
        if isempty(all_ITPC_curves_perNeuron{r}) || isempty(all_tSec{r}) || isempty(all_neuron_resetStrength{r})
            continue;
        end

        tThis = all_tSec{r}(:)';
        X     = all_ITPC_curves_perNeuron{r};
        rs    = all_neuron_resetStrength{r}(:);

        if size(X,2) ~= numel(tCommon) || numel(tThis) ~= numel(tCommon) || any(abs(tThis - tCommon) > 1e-9)
            Xrs = nan(size(X,1), numel(tCommon));
            for n = 1:size(X,1)
                Xrs(n,:) = interp1(tThis, X(n,:), tCommon, 'linear', NaN);
            end
            X = Xrs;
        end

        allITPC  = [allITPC; X]; %#ok<AGROW>
        allReset = [allReset; rs]; %#ok<AGROW>
    end

    if ~isempty(allITPC)
        figPoolITPC = plot_itpc_psth_like(allITPC, tCommon, allReset, ...
            sprintf('All recordings | bout-triggered phase concentration %.02f - %.02f Hz', ...
            frequencyRange(1), frequencyRange(2)));

        if saveOutputs
            saveas(figPoolITPC, fullfile(outFolder, 'pooled_itpc_psth_like_all_neurons.png'));
        end
    end
end

%% =========================
% 8) POOLED PSTH-LIKE ACTIVITY FIGURE (ALL NEURONS)
% =========================

if any(~cellfun(@isempty, all_tSec))
    refIdx  = find(~cellfun(@isempty, all_tSec), 1, 'first');
    tCommon = all_tSec{refIdx}(:)';

    allAct   = [];
    allReset = [];

    for r = 1:numel(all_filtCurves_perNeuron)
        if isempty(all_filtCurves_perNeuron{r}) || isempty(all_tSec{r}) || isempty(all_neuron_resetStrength{r})
            continue;
        end

        tThis = all_tSec{r}(:)';
        X     = all_filtCurves_perNeuron{r};
        rs    = all_neuron_resetStrength{r}(:);

        if size(X,2) ~= numel(tCommon) || numel(tThis) ~= numel(tCommon) || any(abs(tThis - tCommon) > 1e-9)
            Xrs = nan(size(X,1), numel(tCommon));
            for n = 1:size(X,1)
                Xrs(n,:) = interp1(tThis, X(n,:), tCommon, 'linear', NaN);
            end
            X = Xrs;
        end

        allAct   = [allAct; X]; %#ok<AGROW>
        allReset = [allReset; rs]; %#ok<AGROW>
    end

    if ~isempty(allAct)
        figPoolAct = plot_activity_psth_like(allAct, tCommon, allReset, ...
            sprintf('All recordings | bandpassed activity %.02f - %.02f Hz', ...
            frequencyRange(1), frequencyRange(2)), ...
            'Bandpassed \DeltaF/F');
        try
            clim = prctile(allAct(:), [2 98]);
            if all(isfinite(clim)) && clim(1) < clim(2)
                caxis(clim);
            end
        catch
        end

        if saveOutputs
            saveas(figPoolAct, fullfile(outFolder, 'pooled_activity_psth_like_all_neurons.png'));
        end
    end
end

%% =========================
% 9) NEW: BUILD PER-FISH POOLING
% =========================

fishIDsAll = R.fishID;
validFishID = ~cellfun(@isempty, fishIDsAll);
fishList = unique(fishIDsAll(validFishID));

F = struct();
F.fishID = fishList(:);

nFish = numel(fishList);

F.nRecordings            = nan(nFish,1);
F.nNeurons_total         = nan(nFish,1);
F.nNeurons_used          = nan(nFish,1);
F.nBouts_used_total      = nan(nFish,1);

F.pre_ITPC_medianNeuron  = nan(nFish,1);
F.post_ITPC_medianNeuron = nan(nFish,1);
F.delta_ITPC_medianNeuron = nan(nFish,1);
F.resetStrength_medianNeuron = nan(nFish,1);

F.pre_Freq_medianNeuron  = nan(nFish,1);
F.post_Freq_medianNeuron = nan(nFish,1);
F.delta_Freq_medianNeuron = nan(nFish,1);

% store pooled neuron-level data per fish
fish_ITPC_curves         = cell(nFish,1);
fish_freq_curves         = cell(nFish,1);
fish_act_curves          = cell(nFish,1);
fish_tSec                = cell(nFish,1);

fish_preITPC             = cell(nFish,1);
fish_postITPC            = cell(nFish,1);
fish_deltaITPC           = cell(nFish,1);
fish_resetStrength       = cell(nFish,1);

fish_preFreq             = cell(nFish,1);
fish_postFreq            = cell(nFish,1);
fish_deltaFreq           = cell(nFish,1);

fish_stats = struct([]);

for f = 1:nFish
    thisFish = fishList{f};
    idxRecFish = find(strcmp(R.fishID, thisFish));

    F.nRecordings(f) = numel(idxRecFish);
    F.nNeurons_total(f) = nansum(R.nNeurons_total(idxRecFish));
    F.nNeurons_used(f)  = nansum(R.nNeurons_used(idxRecFish));
    F.nBouts_used_total(f) = nansum(R.nBouts_used(idxRecFish));

    refIdx = [];
    for rr = idxRecFish(:)'
        if ~isempty(all_tSec{rr})
            refIdx = rr;
            break;
        end
    end

    if isempty(refIdx)
        warning('Fish %s has no valid time axis. Skipping pooled fish analysis.', thisFish);
        continue;
    end

    tCommon = all_tSec{refIdx}(:)';

    Xitpc = [];
    Xfreq = [];
    Xact  = [];

    v_preITPC   = [];
    v_postITPC  = [];
    v_deltaITPC = [];
    v_reset     = [];

    v_preFreq   = [];
    v_postFreq  = [];
    v_deltaFreq = [];

    for rr = idxRecFish(:)'
        if isempty(all_tSec{rr})
            continue;
        end

        tThis = all_tSec{rr}(:)';

        % --- ITPC curves per neuron
        if ~isempty(all_ITPC_curves_perNeuron{rr})
            X = all_ITPC_curves_perNeuron{rr};
            X = resample_rows_to_common_time(X, tThis, tCommon);
            Xitpc = [Xitpc; X]; %#ok<AGROW>
        end

        % --- freq curves per neuron
        if ~isempty(all_freq_curves_perNeuron{rr})
            X = all_freq_curves_perNeuron{rr};
            X = resample_rows_to_common_time(X, tThis, tCommon);
            Xfreq = [Xfreq; X]; %#ok<AGROW>
        end

        % --- activity curves per neuron
        if ~isempty(all_filtCurves_perNeuron{rr})
            X = all_filtCurves_perNeuron{rr};
            X = resample_rows_to_common_time(X, tThis, tCommon);
            Xact = [Xact; X]; %#ok<AGROW>
        end

        % --- neuron-level summaries
        if ~isempty(all_neuron_preITPC{rr}),   v_preITPC   = [v_preITPC;   all_neuron_preITPC{rr}(:)];   end %#ok<AGROW>
        if ~isempty(all_neuron_postITPC{rr}),  v_postITPC  = [v_postITPC;  all_neuron_postITPC{rr}(:)];  end %#ok<AGROW>
        if ~isempty(all_neuron_deltaITPC{rr}), v_deltaITPC = [v_deltaITPC; all_neuron_deltaITPC{rr}(:)]; end %#ok<AGROW>
        if ~isempty(all_neuron_resetStrength{rr}), v_reset = [v_reset; all_neuron_resetStrength{rr}(:)]; end %#ok<AGROW>

        if ~isempty(all_neuron_preFreq{rr}),   v_preFreq   = [v_preFreq;   all_neuron_preFreq{rr}(:)];   end %#ok<AGROW>
        if ~isempty(all_neuron_postFreq{rr}),  v_postFreq  = [v_postFreq;  all_neuron_postFreq{rr}(:)];  end %#ok<AGROW>
        if ~isempty(all_neuron_deltaFreq{rr}), v_deltaFreq = [v_deltaFreq; all_neuron_deltaFreq{rr}(:)]; end %#ok<AGROW>
    end

    fish_tSec{f}          = tCommon(:);
    fish_ITPC_curves{f}   = Xitpc;
    fish_freq_curves{f}   = Xfreq;
    fish_act_curves{f}    = Xact;

    fish_preITPC{f}       = v_preITPC;
    fish_postITPC{f}      = v_postITPC;
    fish_deltaITPC{f}     = v_deltaITPC;
    fish_resetStrength{f} = v_reset;

    fish_preFreq{f}       = v_preFreq;
    fish_postFreq{f}      = v_postFreq;
    fish_deltaFreq{f}     = v_deltaFreq;

    F.pre_ITPC_medianNeuron(f)      = median(v_preITPC, 'omitnan');
    F.post_ITPC_medianNeuron(f)     = median(v_postITPC, 'omitnan');
    F.delta_ITPC_medianNeuron(f)    = median(v_deltaITPC, 'omitnan');
    F.resetStrength_medianNeuron(f) = median(v_reset, 'omitnan');

    F.pre_Freq_medianNeuron(f)      = median(v_preFreq, 'omitnan');
    F.post_Freq_medianNeuron(f)     = median(v_postFreq, 'omitnan');
    F.delta_Freq_medianNeuron(f)    = median(v_deltaFreq, 'omitnan');

    fish_stats(f).fishID = thisFish;
    fish_stats(f).ITPC_pre_post     = run_paired_signrank(v_preITPC, v_postITPC);
    fish_stats(f).Freq_pre_post     = run_paired_signrank(v_preFreq, v_postFreq);
    fish_stats(f).deltaITPC_vs_zero = run_signrank_vs_zero(v_deltaITPC);
    fish_stats(f).reset_vs_zero     = run_signrank_vs_zero(v_reset);
    fish_stats(f).deltaFreq_vs_zero = run_signrank_vs_zero(v_deltaFreq);

    fprintf('\n---------------- FISH %s ----------------\n', thisFish);
    fprintf('Recordings pooled: %d\n', F.nRecordings(f));
    fprintf('Neurons pooled:    %d\n', numel(v_preITPC));
    fprintf('Bouts pooled sum:  %d\n', F.nBouts_used_total(f));
    fprintf('Neuron-level ITPC pre vs post:   n=%d, p=%.4g, medianΔ=%.4g\n', ...
        fish_stats(f).ITPC_pre_post.n, fish_stats(f).ITPC_pre_post.p, fish_stats(f).ITPC_pre_post.deltaMedian);
    fprintf('Neuron-level Freq pre vs post:   n=%d, p=%.4g, medianΔ=%.4g Hz\n', ...
        fish_stats(f).Freq_pre_post.n, fish_stats(f).Freq_pre_post.p, fish_stats(f).Freq_pre_post.deltaMedian);
    fprintf('Neuron-level reset vs 0:         n=%d, p=%.4g, median=%.4g\n', ...
        fish_stats(f).reset_vs_zero.n, fish_stats(f).reset_vs_zero.p, fish_stats(f).reset_vs_zero.medianValue);

    % -------- PER-FISH ITPC PSTH FIGURE
    if makePerFishITPCFigure && ~isempty(Xitpc)
        try
            figTmp = plot_itpc_psth_like( ...
                Xitpc, tCommon, v_reset, sprintf('%s | pooled neurons', thisFish));
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('perFish_ITPC_PSTH_%s.png', sanitize_filename(thisFish))));
            end
        catch ME
            warning('Could not make per-fish ITPC figure for %s\n%s', thisFish, ME.message);
        end
    end

    % -------- PER-FISH ACTIVITY PSTH FIGURE
    if makePerFishActivityFigure && ~isempty(Xact)
        try
            figTmp = plot_activity_psth_like( ...
                Xact, tCommon, v_reset, sprintf('%s | pooled neurons', thisFish), 'Bandpassed \DeltaF/F');
            try
                clim = prctile(Xact(:), [2 98]);
                if all(isfinite(clim)) && clim(1) < clim(2)
                    caxis(clim);
                end
            catch
            end
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('perFish_activity_PSTH_%s.png', sanitize_filename(thisFish))));
            end
        catch ME
            warning('Could not make per-fish activity figure for %s\n%s', thisFish, ME.message);
        end
    end

    % -------- PER-FISH NEURON-LEVEL STATS FIGURE
    if makePerFishStatsFigure && ~isempty(v_preITPC)
        try
            figTmp = plot_per_entity_stats( ...
                v_preITPC, v_postITPC, v_deltaITPC, v_reset, ...
                v_preFreq, v_postFreq, v_deltaFreq, ...
                sprintf('%s | neuron-level pooled statistics', thisFish));
            if saveOutputs
                saveas(figTmp, fullfile(outFolder, ...
                    sprintf('perFish_stats_neurons_%s.png', sanitize_filename(thisFish))));
            end
        catch ME
            warning('Could not make per-fish stats figure for %s\n%s', thisFish, ME.message);
        end
    end
end

%% =========================
% 10) PER-FISH SUMMARY TABLE
% =========================

T_fish = table( ...
    F.fishID, F.nRecordings, F.nNeurons_total, F.nNeurons_used, F.nBouts_used_total, ...
    F.pre_ITPC_medianNeuron, F.post_ITPC_medianNeuron, F.delta_ITPC_medianNeuron, F.resetStrength_medianNeuron, ...
    F.pre_Freq_medianNeuron, F.post_Freq_medianNeuron, F.delta_Freq_medianNeuron, ...
    'VariableNames', { ...
    'fishID','nRecordings','nNeurons_total','nNeurons_used','nBouts_used_total', ...
    'pre_ITPC_medianNeuron','post_ITPC_medianNeuron','delta_ITPC_medianNeuron','resetStrength_medianNeuron', ...
    'pre_Freq_medianNeuron','post_Freq_medianNeuron','delta_Freq_medianNeuron'});

disp('================ PER-FISH SUMMARY TABLE ================');
disp(T_fish);

if saveOutputs
    writetable(T_fish, fullfile(outFolder, 'summary_bout_reset_ultraslow_perFish.csv'));
end

%% =========================
% 11) SAVE ALL OUTPUT STRUCTS
% =========================

if saveOutputs
    save(fullfile(outFolder, 'summary_bout_reset_ultraslow_ALL.mat'), ...
        'R', 'T', 'F', 'T_fish', 'fish_stats', ...
        'all_ITPC_curves', 'all_freq_curves', 'all_meanFilt_curves', 'all_tSec', ...
        'all_ITPC_curves_perNeuron', 'all_freq_curves_perNeuron', 'all_filtCurves_perNeuron', ...
        'all_neuron_preITPC', 'all_neuron_postITPC', 'all_neuron_deltaITPC', ...
        'all_neuron_resetStrength', 'all_neuron_preFreq', 'all_neuron_postFreq', ...
        'all_neuron_deltaFreq', ...
        'fish_ITPC_curves', 'fish_freq_curves', 'fish_act_curves', 'fish_tSec', ...
        'fish_preITPC', 'fish_postITPC', 'fish_deltaITPC', 'fish_resetStrength', ...
        'fish_preFreq', 'fish_postFreq', 'fish_deltaFreq');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% =========================
% 12) NEW: WITHIN-FISH PAIRED SESSION COMPARISONS
% =========================
% Assumption:
%   Neuron row order is preserved across sessions within a fish.
%   Therefore neuron #i in highact is the same neuron as neuron #i in lowact.
%
% If this is NOT true for your pipeline, these paired comparisons are invalid.

fishIDsAll = R.fishID;
validFishID = ~cellfun(@isempty, fishIDsAll);
fishList = unique(fishIDsAll(validFishID));

pairedSessionSummary = struct([]);

for f = 1:numel(fishList)

    thisFish = fishList{f};
    idxFish = find(strcmp(R.fishID, thisFish));

    if isempty(idxFish)
        continue;
    end

    fishSessionNames = lower(R.sessionLabel(idxFish));

    idxHigh = idxFish(contains(fishSessionNames, 'highact'));
    idxMid  = idxFish(contains(fishSessionNames, 'midact'));
    idxLow  = idxFish(contains(fishSessionNames, 'lowact'));

    fprintf('\n============================================================\n');
    fprintf('WITHIN-FISH PAIRED SESSION COMPARISONS: %s\n', thisFish);

    pairedSessionSummary(f).fishID = thisFish;
    pairedSessionSummary(f).idxHigh = idxHigh;
    pairedSessionSummary(f).idxMid  = idxMid;
    pairedSessionSummary(f).idxLow  = idxLow;

    % If multiple entries exist for a given state, take the first valid one.
    idxHigh = pick_first_valid_recording(idxHigh, all_neuron_preITPC);
    idxMid  = pick_first_valid_recording(idxMid,  all_neuron_preITPC);
    idxLow  = pick_first_valid_recording(idxLow,  all_neuron_preITPC);

    pairedSessionSummary(f).chosenHigh = idxHigh;
    pairedSessionSummary(f).chosenMid  = idxMid;
    pairedSessionSummary(f).chosenLow  = idxLow;

    % ------------------------------------------------------------
    % highact vs lowact
    % ------------------------------------------------------------
    if ~isempty(idxHigh) && ~isempty(idxLow)
        P = build_paired_session_comparison(idxHigh, idxLow, R, ...
            all_neuron_preITPC, all_neuron_postITPC, all_neuron_deltaITPC, all_neuron_resetStrength, ...
            all_neuron_preFreq, all_neuron_postFreq, all_neuron_deltaFreq);

        pairedSessionSummary(f).high_vs_low = P;

        fprintf('  highact vs lowact: matched neurons = %d\n', P.nMatched);

        if P.nMatched >= 2
            figH = plot_paired_session_comparison(P, ...
                sprintf('%s | highact vs lowact | matched neurons', thisFish));
            if saveOutputs
                saveas(figH, fullfile(outFolder, ...
                    sprintf('pairedSessions_%s_high_vs_low.png', sanitize_filename(thisFish))));
            end
        end
    else
        pairedSessionSummary(f).high_vs_low = [];
        fprintf('  highact vs lowact: missing one or both sessions.\n');
    end

    % ------------------------------------------------------------
    % highact vs midact
    % ------------------------------------------------------------
    if ~isempty(idxHigh) && ~isempty(idxMid)
        P = build_paired_session_comparison(idxHigh, idxMid, R, ...
            all_neuron_preITPC, all_neuron_postITPC, all_neuron_deltaITPC, all_neuron_resetStrength, ...
            all_neuron_preFreq, all_neuron_postFreq, all_neuron_deltaFreq);

        pairedSessionSummary(f).high_vs_mid = P;

        fprintf('  highact vs midact: matched neurons = %d\n', P.nMatched);

        if P.nMatched >= 2
            figH = plot_paired_session_comparison(P, ...
                sprintf('%s | highact vs midact | matched neurons', thisFish));
            if saveOutputs
                saveas(figH, fullfile(outFolder, ...
                    sprintf('pairedSessions_%s_high_vs_mid.png', sanitize_filename(thisFish))));
            end
        end
    else
        pairedSessionSummary(f).high_vs_mid = [];
        fprintf('  highact vs midact: missing one or both sessions.\n');
    end

    % ------------------------------------------------------------
    % midact vs lowact
    % ------------------------------------------------------------
    if ~isempty(idxMid) && ~isempty(idxLow)
        P = build_paired_session_comparison(idxMid, idxLow, R, ...
            all_neuron_preITPC, all_neuron_postITPC, all_neuron_deltaITPC, all_neuron_resetStrength, ...
            all_neuron_preFreq, all_neuron_postFreq, all_neuron_deltaFreq);

        pairedSessionSummary(f).mid_vs_low = P;

        fprintf('  midact vs lowact: matched neurons = %d\n', P.nMatched);

        if P.nMatched >= 2
            figH = plot_paired_session_comparison(P, ...
                sprintf('%s | midact vs lowact | matched neurons', thisFish));
            if saveOutputs
                saveas(figH, fullfile(outFolder, ...
                    sprintf('pairedSessions_%s_mid_vs_low.png', sanitize_filename(thisFish))));
            end
        end
    else
        pairedSessionSummary(f).mid_vs_low = [];
        fprintf('  midact vs lowact: missing one or both sessions.\n');
    end
end

if saveOutputs
    save(fullfile(outFolder, 'paired_session_summary_perFish.mat'), 'pairedSessionSummary');
end

%% =========================
% LOCAL FUNCTIONS
% =========================

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
        'boutOnsetFrames', 'boutStartFrames', 'bou_ons_fra', ...
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

function M = compute_bout_reset_metrics(dff, boutOnsetsFrames, fps, frequencyRange, butterOrder, periSec, preStatSec, postStatSec, postPeakSec)
% dff: neurons x time

    [nNeurons, ~] = size(dff);
    periFrames = round(periSec * fps);
    tFrames = (-periFrames:periFrames);
    tSec = tFrames / fps;

    preMask      = tSec >= preStatSec(1)  & tSec <= preStatSec(2);
    postMask     = tSec >= postStatSec(1) & tSec <= postStatSec(2);
    postPeakMask = tSec >= postPeakSec(1) & tSec <= postPeakSec(2);

    if ~any(preMask) || ~any(postMask) || ~any(postPeakMask)
        error('One of the requested time windows is empty. Check periSec and stats windows.');
    end

    Wn = frequencyRange / (fps/2);
    if any(Wn <= 0) || any(Wn >= 1) || Wn(1) >= Wn(2)
        error('Invalid frequencyRange for this fps.');
    end
    [b,a] = butter(butterOrder, Wn, 'bandpass');

    neuron_preITPC       = nan(nNeurons,1);
    neuron_postITPC      = nan(nNeurons,1);
    neuron_deltaITPC     = nan(nNeurons,1);
    neuron_resetStrength = nan(nNeurons,1);

    neuron_preFreq       = nan(nNeurons,1);
    neuron_postFreq      = nan(nNeurons,1);
    neuron_deltaFreq     = nan(nNeurons,1);

    itpcCurves     = nan(nNeurons, numel(tSec));
    freqCurves     = nan(nNeurons, numel(tSec));
    meanFiltCurves = nan(nNeurons, numel(tSec));

    nB = numel(boutOnsetsFrames);

    for n = 1:nNeurons
        x = double(dff(n,:));
        x = filtfilt(b, a, x);

        z = hilbert(x);
        ph = angle(z);

        phUn = unwrap(ph);
        instFreq = [nan diff(phUn)] * fps / (2*pi);  % Hz
        instFreq = smoothdata(instFreq, 'gaussian', max(3, round(fps)));

        phaseMat = nan(nB, numel(tSec));
        freqMat  = nan(nB, numel(tSec));
        filtMat  = nan(nB, numel(tSec));

        for bIdx = 1:nB
            c = boutOnsetsFrames(bIdx);
            idx = c + tFrames;

            phaseMat(bIdx,:) = ph(idx);
            freqMat(bIdx,:)  = instFreq(idx);
            filtMat(bIdx,:)  = x(idx);
        end

        itpc = abs(mean(exp(1i*phaseMat), 1, 'omitnan'));

        preITPC   = mean(itpc(preMask), 'omitnan');
        postITPC  = mean(itpc(postMask), 'omitnan');
        resetStr  = max(itpc(postPeakMask), [], 'omitnan') - preITPC;

        preFreq   = mean(freqMat(:,preMask),  'all', 'omitnan');
        postFreq  = mean(freqMat(:,postMask), 'all', 'omitnan');

        neuron_preITPC(n)       = preITPC;
        neuron_postITPC(n)      = postITPC;
        neuron_deltaITPC(n)     = postITPC - preITPC;
        neuron_resetStrength(n) = resetStr;

        neuron_preFreq(n)       = preFreq;
        neuron_postFreq(n)      = postFreq;
        neuron_deltaFreq(n)     = postFreq - preFreq;

        itpcCurves(n,:)         = itpc;
        freqCurves(n,:)         = mean(freqMat, 1, 'omitnan');
        meanFiltCurves(n,:)     = mean(filtMat, 1, 'omitnan');
    end

    M = struct();
    M.tSec                 = tSec(:);

    M.itpcCurves           = itpcCurves;
    M.freqCurves           = freqCurves;
    M.meanFiltCurves       = meanFiltCurves;

    M.meanITPC             = mean(itpcCurves, 1, 'omitnan').';
    M.meanInstFreq         = mean(freqCurves, 1, 'omitnan').';
    M.meanPopulationFilt   = mean(meanFiltCurves, 1, 'omitnan').';

    M.neuron_preITPC       = neuron_preITPC;
    M.neuron_postITPC      = neuron_postITPC;
    M.neuron_deltaITPC     = neuron_deltaITPC;
    M.neuron_resetStrength = neuron_resetStrength;

    M.neuron_preFreq       = neuron_preFreq;
    M.neuron_postFreq      = neuron_postFreq;
    M.neuron_deltaFreq     = neuron_deltaFreq;
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

function out = run_signrank_vs_zero(x)
    x = x(isfinite(x));

    out = struct('p',NaN,'h',NaN,'signedrank',NaN,'n',numel(x), ...
                 'medianValue',NaN);

    if numel(x) >= 2
        [p,h,stats] = signrank(x);
        out.p = p;
        out.h = h;
        out.signedrank = stats.signedrank;
        out.medianValue = median(x, 'omitnan');
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

function fishID = extract_fish_id(folderPath)
    tok = regexp(lower(folderPath), '(f\d+)', 'tokens', 'once');
    if ~isempty(tok)
        fishID = tok{1};
    else
        parts = regexp(folderPath, '[\\/]', 'split');
        fishID = parts{max(1, numel(parts)-1)};
    end
end

function out = sanitize_filename(in)
    out = regexprep(in, '[^\w\-]', '_');
end

function Xout = resample_rows_to_common_time(X, tThis, tCommon)
    if isempty(X)
        Xout = X;
        return;
    end

    if size(X,2) == numel(tCommon) && numel(tThis) == numel(tCommon) && all(abs(tThis - tCommon) < 1e-9)
        Xout = X;
        return;
    end

    Xout = nan(size(X,1), numel(tCommon));
    for i = 1:size(X,1)
        Xout(i,:) = interp1(tThis, X(i,:), tCommon, 'linear', NaN);
    end
end

function figH = plot_itpc_psth_like(itpcCurves, tSec, sortMetric, recordingLabel)

    if nargin < 4
        recordingLabel = '';
    end

    tSec = tSec(:)';
    nNeurons = size(itpcCurves,1);

    if nargin < 3 || isempty(sortMetric) || numel(sortMetric) ~= nNeurons
        preMask  = tSec >= -40 & tSec <= -5;
        postMask = tSec >= 5   & tSec <= 40;
        sortMetric = mean(itpcCurves(:,postMask),2,'omitnan') - ...
                     mean(itpcCurves(:,preMask),2,'omitnan');
    end

    [~, ord] = sort(sortMetric(:), 'descend');
    X = itpcCurves(ord,:);

    mu  = mean(itpcCurves, 1, 'omitnan');
    sem = std(itpcCurves, 0, 1, 'omitnan') ./ sqrt(size(itpcCurves,1));

    figH = figure('Color','w', 'Position', [100 100 560 760]);
    tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

    nexttile; hold on;
    fill([tSec fliplr(tSec)], ...
         [mu-sem fliplr(mu+sem)], ...
         [0.8 0.8 0.8], ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.7);
    plot(tSec, mu, 'k', 'LineWidth', 2);
    xline(0, '--k', 'LineWidth', 1);

    ylabel('ITPC');
    title('Bout-triggered phase concentration');
    xlim([min(tSec) max(tSec)]);
    box off;

    nexttile;
    imagesc(tSec, 1:size(X,1), X);
    set(gca, 'YDir', 'normal');
    xline(0, '--k', 'LineWidth', 1);

    xlabel('Time from bout onset (s)');
    ylabel('Neuron #');
    title('Neurons sorted by reset ITPC');

    colormap(gca, hot);
    cb = colorbar;
    ylabel(cb, 'ITPC');

    xlim([min(tSec) max(tSec)]);

    if ~isempty(recordingLabel)
        sgtitle(recordingLabel, 'Interpreter', 'none');
    end
end

function figH = plot_activity_psth_like(activityCurves, tSec, sortMetric, recordingLabel, ylabTxt)

    if nargin < 5 || isempty(ylabTxt)
        ylabTxt = 'Filtered \DeltaF/F';
    end
    if nargin < 4
        recordingLabel = '';
    end

    tSec = tSec(:)';
    nNeurons = size(activityCurves,1);

    if nargin < 3 || isempty(sortMetric) || numel(sortMetric) ~= nNeurons
        sortMetric = nanmean(activityCurves, 2);
    end

    [~, ord] = sort(sortMetric(:), 'descend');
    X = activityCurves(ord,:);

    mu  = mean(activityCurves, 1, 'omitnan');
    sem = std(activityCurves, 0, 1, 'omitnan') ./ sqrt(size(activityCurves,1));

    figH = figure('Color','w', 'Position', [100 100 560 760]);
    tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

    nexttile; hold on;
    fill([tSec fliplr(tSec)], ...
         [mu-sem fliplr(mu+sem)], ...
         [0.8 0.8 0.8], ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.7);
    plot(tSec, mu, 'k', 'LineWidth', 2);
    xline(0, '--k', 'LineWidth', 1);

    ylabel(ylabTxt);
    title('Bout-triggered neural activity');
    xlim([min(tSec) max(tSec)]);
    box off;

    nexttile;
    imagesc(tSec, 1:size(X,1), X);
    set(gca, 'YDir', 'normal');
    xline(0, '--k', 'LineWidth', 1);

    xlabel('Time from bout onset (s)');
    ylabel('Neuron #');
    title('Neurons sorted by reset ITPC');

    colormap(gca, hot);
    cb = colorbar;
    ylabel(cb, ylabTxt);

    xlim([min(tSec) max(tSec)]);

    if ~isempty(recordingLabel)
        sgtitle(recordingLabel, 'Interpreter', 'none');
    end
end

function figH = plot_per_entity_stats(preITPC, postITPC, deltaITPC, resetStrength, preFreq, postFreq, deltaFreq, ttl)

    figH = figure('Color','w', 'Position', [100 100 1350 780]);
    tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

    statsITPC  = run_paired_signrank(preITPC, postITPC);
    statsFreq  = run_paired_signrank(preFreq, postFreq);
    statsDItpc = run_signrank_vs_zero(deltaITPC);
    statsReset = run_signrank_vs_zero(resetStrength);
    statsDFreq = run_signrank_vs_zero(deltaFreq);

    % 1) ITPC pre/post
    nexttile; hold on;
    plot_box_scatter_pair(preITPC, postITPC, {'Pre','Post'});
    ylabel('ITPC');
    title(sprintf('ITPC pre vs post\nn=%d, p=%.3g %s', ...
        statsITPC.n, statsITPC.p, p_to_stars(statsITPC.p)));
    add_sig_bar(gca, [1 2], [preITPC(:); postITPC(:)], statsITPC.p);

    % 2) Freq pre/post
    nexttile; hold on;
    plot_box_scatter_pair(preFreq, postFreq, {'Pre','Post'});
    ylabel('Instantaneous frequency (Hz)');
    title(sprintf('Frequency pre vs post\nn=%d, p=%.3g %s', ...
        statsFreq.n, statsFreq.p, p_to_stars(statsFreq.p)));
    add_sig_bar(gca, [1 2], [preFreq(:); postFreq(:)], statsFreq.p);

    % 3) Delta ITPC
    nexttile; hold on;
    plot_box_scatter_single(deltaITPC, 1);
    yline(0,'--k');
    xlim([0.5 1.5]);
    set(gca,'XTick',1,'XTickLabel',{'Post-Pre'});
    ylabel('\Delta ITPC');
    title(sprintf('\x0394ITPC vs 0\nn=%d, p=%.3g %s', ...
        statsDItpc.n, statsDItpc.p, p_to_stars(statsDItpc.p)));

    % 4) Reset strength
    nexttile; hold on;
    plot_box_scatter_single(resetStrength, 1);
    yline(0,'--k');
    xlim([0.5 1.5]);
    set(gca,'XTick',1,'XTickLabel',{'Reset'});
    ylabel('Reset strength');
    title(sprintf('Reset strength vs 0\nn=%d, p=%.3g %s', ...
        statsReset.n, statsReset.p, p_to_stars(statsReset.p)));

    % 5) Delta frequency
    nexttile; hold on;
    plot_box_scatter_single(deltaFreq, 1);
    yline(0,'--k');
    xlim([0.5 1.5]);
    set(gca,'XTick',1,'XTickLabel',{'Post-Pre'});
    ylabel('\Delta instantaneous frequency (Hz)');
    title(sprintf('\x0394Freq vs 0\nn=%d, p=%.3g %s', ...
        statsDFreq.n, statsDFreq.p, p_to_stars(statsDFreq.p)));

    % 6) text panel
    nexttile; axis off;
    text(0, 0.95, ttl, 'FontWeight','bold', 'FontSize', 13, 'Interpreter','none');
    text(0, 0.78, sprintf('ITPC pre-post p = %.4g', statsITPC.p), 'FontSize', 11);
    text(0, 0.66, sprintf('Freq pre-post p = %.4g', statsFreq.p), 'FontSize', 11);
    text(0, 0.54, sprintf('Delta ITPC vs 0 p = %.4g', statsDItpc.p), 'FontSize', 11);
    text(0, 0.42, sprintf('Reset vs 0 p = %.4g', statsReset.p), 'FontSize', 11);
    text(0, 0.30, sprintf('Delta Freq vs 0 p = %.4g', statsDFreq.p), 'FontSize', 11);

    sgtitle(ttl, 'Interpreter','none');
end

function plot_box_scatter_pair(x1, x2, labels)
    x1 = x1(isfinite(x1));
    x2 = x2(isfinite(x2));

    boxchart(ones(numel(x1),1)*1, x1, 'BoxFaceAlpha', 0.3);
    boxchart(ones(numel(x2),1)*2, x2, 'BoxFaceAlpha', 0.3);

    jitter1 = 1 + 0.08*(rand(numel(x1),1)-0.5);
    jitter2 = 2 + 0.08*(rand(numel(x2),1)-0.5);

    scatter(jitter1, x1, 18, 'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);
    scatter(jitter2, x2, 18, 'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);

    set(gca,'XTick',[1 2],'XTickLabel',labels);
    xlim([0.5 2.5]);
    box off;
end

function plot_box_scatter_single(x, xpos)
    x = x(isfinite(x));

    boxchart(ones(numel(x),1)*xpos, x, 'BoxFaceAlpha', 0.3);
    jitter = xpos + 0.08*(rand(numel(x),1)-0.5);
    scatter(jitter, x, 18, 'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);

    box off;
end
function sessionLabel = extract_session_label(folderPath)
    parts = regexp(folderPath, '[\\/]', 'split');
    parts = parts(~cellfun(@isempty, parts));

    sessionLabel = '';

    for i = numel(parts):-1:1
        p = lower(parts{i});
        if contains(p, 'highact') || contains(p, 'midact') || contains(p, 'lowact')
            sessionLabel = p;
            return;
        end
    end

    if numel(parts) >= 1
        sessionLabel = lower(parts{end});
    end
end

function idx = pick_first_valid_recording(idxList, metricCell)
    idx = [];
    for k = 1:numel(idxList)
        ii = idxList(k);
        if ii >= 1 && ii <= numel(metricCell) && ~isempty(metricCell{ii})
            idx = ii;
            return;
        end
    end
end

function P = build_paired_session_comparison(idxA, idxB, R, ...
    all_neuron_preITPC, all_neuron_postITPC, all_neuron_deltaITPC, all_neuron_resetStrength, ...
    all_neuron_preFreq, all_neuron_postFreq, all_neuron_deltaFreq)

    P = struct();

    P.idxA = idxA;
    P.idxB = idxB;
    P.labelA = R.sessionLabel{idxA};
    P.labelB = R.sessionLabel{idxB};
    P.fishID = R.fishID{idxA};

    % pull vectors
    preITPC_A   = all_neuron_preITPC{idxA}(:);
    postITPC_A  = all_neuron_postITPC{idxA}(:);
    deltaITPC_A = all_neuron_deltaITPC{idxA}(:);
    reset_A     = all_neuron_resetStrength{idxA}(:);

    preFreq_A   = all_neuron_preFreq{idxA}(:);
    postFreq_A  = all_neuron_postFreq{idxA}(:);
    deltaFreq_A = all_neuron_deltaFreq{idxA}(:);

    preITPC_B   = all_neuron_preITPC{idxB}(:);
    postITPC_B  = all_neuron_postITPC{idxB}(:);
    deltaITPC_B = all_neuron_deltaITPC{idxB}(:);
    reset_B     = all_neuron_resetStrength{idxB}(:);

    preFreq_B   = all_neuron_preFreq{idxB}(:);
    postFreq_B  = all_neuron_postFreq{idxB}(:);
    deltaFreq_B = all_neuron_deltaFreq{idxB}(:);

    % match by common neuron index
    nCommon = min([ ...
        numel(preITPC_A), numel(preITPC_B), ...
        numel(postITPC_A), numel(postITPC_B), ...
        numel(deltaITPC_A), numel(deltaITPC_B), ...
        numel(reset_A), numel(reset_B), ...
        numel(preFreq_A), numel(preFreq_B), ...
        numel(postFreq_A), numel(postFreq_B), ...
        numel(deltaFreq_A), numel(deltaFreq_B)]);

    if isempty(nCommon) || nCommon < 1
        P.nMatched = 0;
        return;
    end

    preITPC_A   = preITPC_A(1:nCommon);
    postITPC_A  = postITPC_A(1:nCommon);
    deltaITPC_A = deltaITPC_A(1:nCommon);
    reset_A     = reset_A(1:nCommon);
    preFreq_A   = preFreq_A(1:nCommon);
    postFreq_A  = postFreq_A(1:nCommon);
    deltaFreq_A = deltaFreq_A(1:nCommon);

    preITPC_B   = preITPC_B(1:nCommon);
    postITPC_B  = postITPC_B(1:nCommon);
    deltaITPC_B = deltaITPC_B(1:nCommon);
    reset_B     = reset_B(1:nCommon);
    preFreq_B   = preFreq_B(1:nCommon);
    postFreq_B  = postFreq_B(1:nCommon);
    deltaFreq_B = deltaFreq_B(1:nCommon);

    % valid paired masks per metric
    valid_preITPC   = isfinite(preITPC_A)   & isfinite(preITPC_B);
    valid_postITPC  = isfinite(postITPC_A)  & isfinite(postITPC_B);
    valid_deltaITPC = isfinite(deltaITPC_A) & isfinite(deltaITPC_B);
    valid_reset     = isfinite(reset_A)     & isfinite(reset_B);
    valid_preFreq   = isfinite(preFreq_A)   & isfinite(preFreq_B);
    valid_postFreq  = isfinite(postFreq_A)  & isfinite(postFreq_B);
    valid_deltaFreq = isfinite(deltaFreq_A) & isfinite(deltaFreq_B);

    P.nMatched = nCommon;

    P.preITPC_A   = preITPC_A(valid_preITPC);
    P.preITPC_B   = preITPC_B(valid_preITPC);

    P.postITPC_A  = postITPC_A(valid_postITPC);
    P.postITPC_B  = postITPC_B(valid_postITPC);

    P.deltaITPC_A = deltaITPC_A(valid_deltaITPC);
    P.deltaITPC_B = deltaITPC_B(valid_deltaITPC);

    P.reset_A     = reset_A(valid_reset);
    P.reset_B     = reset_B(valid_reset);

    P.preFreq_A   = preFreq_A(valid_preFreq);
    P.preFreq_B   = preFreq_B(valid_preFreq);

    P.postFreq_A  = postFreq_A(valid_postFreq);
    P.postFreq_B  = postFreq_B(valid_postFreq);

    P.deltaFreq_A = deltaFreq_A(valid_deltaFreq);
    P.deltaFreq_B = deltaFreq_B(valid_deltaFreq);

    P.stats_preITPC   = run_paired_signrank(P.preITPC_A,   P.preITPC_B);
    P.stats_postITPC  = run_paired_signrank(P.postITPC_A,  P.postITPC_B);
    P.stats_deltaITPC = run_paired_signrank(P.deltaITPC_A, P.deltaITPC_B);
    P.stats_reset     = run_paired_signrank(P.reset_A,     P.reset_B);

    P.stats_preFreq   = run_paired_signrank(P.preFreq_A,   P.preFreq_B);
    P.stats_postFreq  = run_paired_signrank(P.postFreq_A,  P.postFreq_B);
    P.stats_deltaFreq = run_paired_signrank(P.deltaFreq_A, P.deltaFreq_B);
end

function figH = plot_paired_session_comparison(P, figTitle)

    figH = figure('Color','w', 'Position', [100 100 1500 820]);
    tiledlayout(2,4, 'TileSpacing','compact', 'Padding','compact');

    labs = {P.labelA, P.labelB};

    % 1 pre ITPC
    nexttile; hold on;
    plot_paired_points(P.preITPC_A, P.preITPC_B, labs);
    ylabel('Pre ITPC');
    title(sprintf('Pre ITPC\nn=%d, p=%.3g %s', ...
        P.stats_preITPC.n, P.stats_preITPC.p, p_to_stars(P.stats_preITPC.p)));
    add_sig_bar(gca, [1 2], [P.preITPC_A(:); P.preITPC_B(:)], P.stats_preITPC.p);

    % 2 post ITPC
    nexttile; hold on;
    plot_paired_points(P.postITPC_A, P.postITPC_B, labs);
    ylabel('Post ITPC');
    title(sprintf('Post ITPC\nn=%d, p=%.3g %s', ...
        P.stats_postITPC.n, P.stats_postITPC.p, p_to_stars(P.stats_postITPC.p)));
    add_sig_bar(gca, [1 2], [P.postITPC_A(:); P.postITPC_B(:)], P.stats_postITPC.p);

    % 3 delta ITPC
    nexttile; hold on;
    plot_paired_points(P.deltaITPC_A, P.deltaITPC_B, labs);
    ylabel('\Delta ITPC');
    title(sprintf('\x0394 ITPC\nn=%d, p=%.3g %s', ...
        P.stats_deltaITPC.n, P.stats_deltaITPC.p, p_to_stars(P.stats_deltaITPC.p)));
    add_sig_bar(gca, [1 2], [P.deltaITPC_A(:); P.deltaITPC_B(:)], P.stats_deltaITPC.p);

    % 4 reset strength
    nexttile; hold on;
    plot_paired_points(P.reset_A, P.reset_B, labs);
    ylabel('Reset strength');
    title(sprintf('Reset strength\nn=%d, p=%.3g %s', ...
        P.stats_reset.n, P.stats_reset.p, p_to_stars(P.stats_reset.p)));
    add_sig_bar(gca, [1 2], [P.reset_A(:); P.reset_B(:)], P.stats_reset.p);

    % 5 pre freq
    nexttile; hold on;
    plot_paired_points(P.preFreq_A, P.preFreq_B, labs);
    ylabel('Pre inst. freq (Hz)');
    title(sprintf('Pre frequency\nn=%d, p=%.3g %s', ...
        P.stats_preFreq.n, P.stats_preFreq.p, p_to_stars(P.stats_preFreq.p)));
    add_sig_bar(gca, [1 2], [P.preFreq_A(:); P.preFreq_B(:)], P.stats_preFreq.p);

    % 6 post freq
    nexttile; hold on;
    plot_paired_points(P.postFreq_A, P.postFreq_B, labs);
    ylabel('Post inst. freq (Hz)');
    title(sprintf('Post frequency\nn=%d, p=%.3g %s', ...
        P.stats_postFreq.n, P.stats_postFreq.p, p_to_stars(P.stats_postFreq.p)));
    add_sig_bar(gca, [1 2], [P.postFreq_A(:); P.postFreq_B(:)], P.stats_postFreq.p);

    % 7 delta freq
    nexttile; hold on;
    plot_paired_points(P.deltaFreq_A, P.deltaFreq_B, labs);
    ylabel('\Delta inst. freq (Hz)');
    title(sprintf('\x0394 frequency\nn=%d, p=%.3g %s', ...
        P.stats_deltaFreq.n, P.stats_deltaFreq.p, p_to_stars(P.stats_deltaFreq.p)));
    add_sig_bar(gca, [1 2], [P.deltaFreq_A(:); P.deltaFreq_B(:)], P.stats_deltaFreq.p);

    % 8 text panel
    nexttile; axis off;
    text(0, 0.95, figTitle, 'FontWeight','bold', 'FontSize', 13, 'Interpreter','none');
    text(0, 0.82, sprintf('Fish: %s', P.fishID), 'FontSize', 11, 'Interpreter','none');
    text(0, 0.72, sprintf('A = %s', P.labelA), 'FontSize', 11, 'Interpreter','none');
    text(0, 0.64, sprintf('B = %s', P.labelB), 'FontSize', 11, 'Interpreter','none');
    text(0, 0.54, sprintf('Matched neurons by row index: %d', P.nMatched), 'FontSize', 11);
    text(0, 0.42, 'Pairing assumption: same neuron order across sessions', 'FontSize', 10);

    sgtitle(figTitle, 'Interpreter','none');
end

function plot_paired_points(xA, xB, labels)

    idx = isfinite(xA) & isfinite(xB);
    xA = xA(idx);
    xB = xB(idx);

    if isempty(xA)
        return;
    end

    boxchart(ones(numel(xA),1)*1, xA, 'BoxFaceAlpha', 0.25);
    boxchart(ones(numel(xB),1)*2, xB, 'BoxFaceAlpha', 0.25);

    for i = 1:numel(xA)
        plot([1 2], [xA(i) xB(i)], '-', 'Color', [0.7 0.7 0.7]);
    end

    scatter(ones(numel(xA),1)*1, xA, 18, 'filled', 'MarkerFaceAlpha',0.7);
    scatter(ones(numel(xB),1)*2, xB, 18, 'filled', 'MarkerFaceAlpha',0.7);

    set(gca, 'XTick', [1 2], 'XTickLabel', labels);
    xlim([0.5 2.5]);
    box off;
end

function figH = plot_per_recording_neuron_stats(preITPC, postITPC, preFreq, postFreq, deltaITPC, resetStrength, deltaFreq, recordingLabel)

    statsITPC  = run_paired_signrank(preITPC, postITPC);
    statsFreq  = run_paired_signrank(preFreq, postFreq);
    statsDItpc = run_signrank_vs_zero(deltaITPC);
    statsReset = run_signrank_vs_zero(resetStrength);
    statsDFreq = run_signrank_vs_zero(deltaFreq);

    figH = figure('Color','w', 'Position', [100 100 1400 760]);
    tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

    % -------------------------------------------------
    % 1) Phase concentration: pre vs post
    % -------------------------------------------------
    nexttile; hold on;

    x1 = preITPC(:);
    x2 = postITPC(:);
    valid = isfinite(x1) & isfinite(x2);
    x1 = x1(valid);
    x2 = x2(valid);

    boxchart(ones(numel(x1),1)*1, x1, 'BoxFaceAlpha', 0.25);
    boxchart(ones(numel(x2),1)*2, x2, 'BoxFaceAlpha', 0.25);

    jitter1 = 1 + 0.08*(rand(numel(x1),1)-0.5);
    jitter2 = 2 + 0.08*(rand(numel(x2),1)-0.5);

    for i = 1:numel(x1)
        plot([jitter1(i) jitter2(i)], [x1(i) x2(i)], '-', 'Color', [0.8 0.8 0.8]);
    end

    scatter(jitter1, x1, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    scatter(jitter2, x2, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);

    set(gca,'XTick',[1 2],'XTickLabel',{'Pre-bout','Post-bout'});
    xlim([0.5 2.5]);
    ylabel('ITPC');
    title(sprintf('Phase concentration\nn=%d neurons, p=%.3g %s', ...
        statsITPC.n, statsITPC.p, p_to_stars(statsITPC.p)));
    add_sig_bar(gca, [1 2], [x1; x2], statsITPC.p);
    box off;

    % -------------------------------------------------
    % 2) Inst frequency: pre vs post
    % -------------------------------------------------
    nexttile; hold on;

    x1 = preFreq(:);
    x2 = postFreq(:);
    valid = isfinite(x1) & isfinite(x2);
    x1 = x1(valid);
    x2 = x2(valid);

    boxchart(ones(numel(x1),1)*1, x1, 'BoxFaceAlpha', 0.25);
    boxchart(ones(numel(x2),1)*2, x2, 'BoxFaceAlpha', 0.25);

    jitter1 = 1 + 0.08*(rand(numel(x1),1)-0.5);
    jitter2 = 2 + 0.08*(rand(numel(x2),1)-0.5);

    for i = 1:numel(x1)
        plot([jitter1(i) jitter2(i)], [x1(i) x2(i)], '-', 'Color', [0.8 0.8 0.8]);
    end

    scatter(jitter1, x1, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    scatter(jitter2, x2, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);

    set(gca,'XTick',[1 2],'XTickLabel',{'Pre-bout','Post-bout'});
    xlim([0.5 2.5]);
    ylabel('Instantaneous frequency (Hz)');
    title(sprintf('Instantaneous frequency\nn=%d neurons, p=%.3g %s', ...
        statsFreq.n, statsFreq.p, p_to_stars(statsFreq.p)));
    add_sig_bar(gca, [1 2], [x1; x2], statsFreq.p);
    box off;

    % -------------------------------------------------
    % 3) Delta ITPC vs 0
    % -------------------------------------------------
    nexttile; hold on;

    x = deltaITPC(:);
    x = x(isfinite(x));

    boxchart(ones(numel(x),1), x, 'BoxFaceAlpha', 0.25);
    jitter = 1 + 0.08*(rand(numel(x),1)-0.5);
    scatter(jitter, x, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    yline(0, '--k');

    set(gca,'XTick',1,'XTickLabel',{'Post-Pre'});
    xlim([0.5 1.5]);
    ylabel('\Delta ITPC');
    title(sprintf('\x0394 ITPC vs 0\nn=%d neurons, p=%.3g %s', ...
        statsDItpc.n, statsDItpc.p, p_to_stars(statsDItpc.p)));
    box off;

    % -------------------------------------------------
    % 4) Reset strength vs 0
    % -------------------------------------------------
    nexttile; hold on;

    x = resetStrength(:);
    x = x(isfinite(x));

    boxchart(ones(numel(x),1), x, 'BoxFaceAlpha', 0.25);
    jitter = 1 + 0.08*(rand(numel(x),1)-0.5);
    scatter(jitter, x, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    yline(0, '--k');

    set(gca,'XTick',1,'XTickLabel',{'Reset'});
    xlim([0.5 1.5]);
    ylabel('Reset strength');
    title(sprintf('Reset strength vs 0\nn=%d neurons, p=%.3g %s', ...
        statsReset.n, statsReset.p, p_to_stars(statsReset.p)));
    box off;

    % -------------------------------------------------
    % 5) Delta frequency vs 0
    % -------------------------------------------------
    nexttile; hold on;

    x = deltaFreq(:);
    x = x(isfinite(x));

    boxchart(ones(numel(x),1), x, 'BoxFaceAlpha', 0.25);
    jitter = 1 + 0.08*(rand(numel(x),1)-0.5);
    scatter(jitter, x, 18, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
    yline(0, '--k');

    set(gca,'XTick',1,'XTickLabel',{'Post-Pre'});
    xlim([0.5 1.5]);
    ylabel('\Delta frequency (Hz)');
    title(sprintf('\x0394 frequency vs 0\nn=%d neurons, p=%.3g %s', ...
        statsDFreq.n, statsDFreq.p, p_to_stars(statsDFreq.p)));
    box off;

    % -------------------------------------------------
    % 6) Text panel
    % -------------------------------------------------
    nexttile; axis off;
    text(0, 0.95, recordingLabel, 'FontWeight','bold', 'FontSize', 13, 'Interpreter','none');
    text(0, 0.80, sprintf('Each scatter point = one neuron'), 'FontSize', 11);
    text(0, 0.68, sprintf('Phase concentration pre-post p = %.4g', statsITPC.p), 'FontSize', 11);
    text(0, 0.56, sprintf('Frequency pre-post p = %.4g', statsFreq.p), 'FontSize', 11);
    text(0, 0.44, sprintf('Delta ITPC vs 0 p = %.4g', statsDItpc.p), 'FontSize', 11);
    text(0, 0.32, sprintf('Reset strength vs 0 p = %.4g', statsReset.p), 'FontSize', 11);
    text(0, 0.20, sprintf('Delta frequency vs 0 p = %.4g', statsDFreq.p), 'FontSize', 11);

    sgtitle(sprintf('%s | neuron-level bout-reset statistics', recordingLabel), 'Interpreter','none');
end