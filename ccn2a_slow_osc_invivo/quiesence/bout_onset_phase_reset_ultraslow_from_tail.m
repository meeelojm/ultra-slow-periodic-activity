%% =========================================================================
% bout_onset_phase_reset_ultraslow_from_tail.m
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
makePerRecordingITPCFigure     = false;
makePerRecordingActivityFigure = false;

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'bout_reset_ultraslow_results');

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
R.recordingLabel           = cell(nR,1);
R.folder                   = cell(nR,1);

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
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.nBouts_total, R.nBouts_used, ...
    R.pre_ITPC_mean, R.post_ITPC_mean, R.delta_ITPC_mean, R.resetStrength_mean, ...
    R.pre_instFreq_mean, R.post_instFreq_mean, R.delta_instFreq_mean, ...
    R.pre_ITPC_median, R.post_ITPC_median, R.delta_ITPC_median, R.resetStrength_median, ...
    R.pre_instFreq_median, R.post_instFreq_median, R.delta_instFreq_median, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used', ...
    'fps','nFrames','durationMin', ...
    'nBouts_total','nBouts_used', ...
    'pre_ITPC_mean','post_ITPC_mean','delta_ITPC_mean','resetStrength_mean', ...
    'pre_instFreq_mean','post_instFreq_mean','delta_instFreq_mean', ...
    'pre_ITPC_median','post_ITPC_median','delta_ITPC_median','resetStrength_median', ...
    'pre_instFreq_median','post_instFreq_median','delta_instFreq_median'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_bout_reset_ultraslow.csv'));
    save(fullfile(outFolder, 'summary_bout_reset_ultraslow.mat'), ...
        'R', 'T', ...
        'all_ITPC_curves', 'all_freq_curves', 'all_meanFilt_curves', 'all_tSec', ...
        'all_ITPC_curves_perNeuron', 'all_freq_curves_perNeuron', 'all_filtCurves_perNeuron', ...
        'all_neuron_preITPC', 'all_neuron_postITPC', 'all_neuron_deltaITPC', ...
        'all_neuron_resetStrength', 'all_neuron_preFreq', 'all_neuron_postFreq', ...
        'all_neuron_deltaFreq');
end

%% =========================
% 5) PAIRED STATISTICS
% =========================

statsSummary = struct();

statsSummary.ITPC_pre_post      = run_paired_signrank(R.pre_ITPC_median, R.post_ITPC_median);
statsSummary.instFreq_pre_post  = run_paired_signrank(R.pre_instFreq_median, R.post_instFreq_median);

statsSummary.deltaITPC_vs_zero      = run_signrank_vs_zero(R.delta_ITPC_median);
statsSummary.resetStrength_vs_zero  = run_signrank_vs_zero(R.resetStrength_median);
statsSummary.deltaFreq_vs_zero      = run_signrank_vs_zero(R.delta_instFreq_median);

disp('================ BOUT RESET STATISTICS ================')
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
    save(fullfile(outFolder, 'stats_summary_bout_reset_ultraslow.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURE
% =========================

validRec = isfinite(R.pre_ITPC_median) & isfinite(R.post_ITPC_median) & ...
           isfinite(R.pre_instFreq_median) & isfinite(R.post_instFreq_median);

if any(validRec)

    figR = figure('Name','Bout-triggered ultraslow reset', ...
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

    sgtitle('Bout-triggered reset of ultraslow neural activity');

    if saveOutputs
        saveas(figR, fullfile(outFolder, 'summary_bout_reset_ultraslow.png'));
    end
else
    warning('No valid recordings available for summary figure.');
end

%% =========================
% 7) POOLED PSTH-LIKE ITPC FIGURE
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
            sprintf('All recordings | bout-triggered phase concentration %.02f - %.02f', ...
            frequencyRange(1), frequencyRange(2)));

        if saveOutputs
            saveas(figPoolITPC, fullfile(outFolder, 'pooled_itpc_psth_like.png'));
        end
    end
end

%% =========================
% 8) POOLED PSTH-LIKE ACTIVITY FIGURE
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
            'Bandpassed \DeltaF/F'); caxis([-10 50])

        if saveOutputs
            saveas(figPoolAct, fullfile(outFolder, 'pooled_activity_psth_like.png'));
        end
    end
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
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

function out = sanitize_filename(in)
    out = regexprep(in, '[^\w\-]', '_');
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

    figH = figure('Color','w', 'Position', [100 100 520 720]);
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

    figH = figure('Color','w', 'Position', [100 100 520 720]);
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