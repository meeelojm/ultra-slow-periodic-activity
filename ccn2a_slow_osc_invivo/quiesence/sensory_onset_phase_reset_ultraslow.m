%% =========================================================================
% sensory_onset_phase_reset_ultraslow.m
%
% SELF-CONTAINED SCRIPT
%
% What it does:
%   1) Recursively finds recordings that contain:
%        - dffs_repact_respcells.mat
%   2) Loads neural dff and fps
%   3) Extracts sensory onset frames from variables related to:
%        - eve_two_pho
%        - optional stimulus type labels
%   4) Splits events into 3 categories:
%        - light
%        - tap
%        - lighttap
%   5) Bandpass-filters each neuron in an ultraslow band
%   6) Uses Hilbert transform to get phase and instantaneous frequency
%   7) Aligns neurons to sensory onsets
%   8) Quantifies:
%        - ITPC
%        - pre vs post ITPC
%        - pre vs post instantaneous frequency
%        - reset strength
%   9) Builds summary tables, paired stats, and pooled PSTH-like figures
%
% IMPORTANT:
%   This script assumes sensory onset information is stored in the same
%   .mat file or in a parseable format accessible from that file.
%
% =========================================================================

clear; clc; close all;

%% =========================
% 0) USER SETTINGS
% =========================

rootFolder = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low';
dffFileName = 'dffs_repact_respcells.mat';

frequencyRange = [0.01 0.05];   % Hz

periSec     = 120;
preStatSec  = [-40 -5];
postStatSec = [5 40];
postPeakSec = [0 40];

minEventCount     = 24;
minInterEventSec  = 20;
maxNanFracPerNeuron = 0.20;
butterOrder = 3;

saveOutputs = true;
outFolder   = fullfile(rootFolder, 'sensory_reset_ultraslow_results');

eventTypes = {'light','tap','lighttap'};

set(0, 'DefaultFigureColor', [1 1 1]);

if saveOutputs && ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%% =========================
% 1) FIND RECORDINGS
% =========================

fprintf('\nSearching for recordings under:\n%s\n\n', rootFolder);

dffFiles = dir(fullfile(rootFolder, '**', dffFileName));
if isempty(dffFiles)
    error('No %s files found under rootFolder.', dffFileName);
end

recordings = struct('folder', {}, 'dffPath', {});
for i = 1:numel(dffFiles)
    recordings(end+1).folder = dffFiles(i).folder; %#ok<SAGROW>
    recordings(end).dffPath  = fullfile(dffFiles(i).folder, dffFiles(i).name);
end

nR = numel(recordings);
fprintf('Found %d candidate recordings.\n', nR);

%% =========================
% 2) PREALLOCATE RESULTS
% =========================

Res = struct();

for e = 1:numel(eventTypes)
    ev = eventTypes{e};

    Res.(ev).R = struct();
    Res.(ev).R.recordingLabel           = cell(nR,1);
    Res.(ev).R.folder                   = cell(nR,1);

    Res.(ev).R.nNeurons_total           = nan(nR,1);
    Res.(ev).R.nNeurons_used            = nan(nR,1);

    Res.(ev).R.fps                      = nan(nR,1);
    Res.(ev).R.nFrames                  = nan(nR,1);
    Res.(ev).R.durationMin              = nan(nR,1);

    Res.(ev).R.nEvents_total            = nan(nR,1);
    Res.(ev).R.nEvents_used             = nan(nR,1);

    Res.(ev).R.pre_ITPC_mean            = nan(nR,1);
    Res.(ev).R.post_ITPC_mean           = nan(nR,1);
    Res.(ev).R.delta_ITPC_mean          = nan(nR,1);
    Res.(ev).R.resetStrength_mean       = nan(nR,1);

    Res.(ev).R.pre_instFreq_mean        = nan(nR,1);
    Res.(ev).R.post_instFreq_mean       = nan(nR,1);
    Res.(ev).R.delta_instFreq_mean      = nan(nR,1);

    Res.(ev).R.pre_ITPC_median          = nan(nR,1);
    Res.(ev).R.post_ITPC_median         = nan(nR,1);
    Res.(ev).R.delta_ITPC_median        = nan(nR,1);
    Res.(ev).R.resetStrength_median     = nan(nR,1);

    Res.(ev).R.pre_instFreq_median      = nan(nR,1);
    Res.(ev).R.post_instFreq_median     = nan(nR,1);
    Res.(ev).R.delta_instFreq_median    = nan(nR,1);

    Res.(ev).all_ITPC_curves            = cell(nR,1);
    Res.(ev).all_freq_curves            = cell(nR,1);
    Res.(ev).all_meanFilt_curves        = cell(nR,1);
    Res.(ev).all_tSec                   = cell(nR,1);

    Res.(ev).all_ITPC_curves_perNeuron  = cell(nR,1);
    Res.(ev).all_freq_curves_perNeuron  = cell(nR,1);
    Res.(ev).all_filtCurves_perNeuron   = cell(nR,1);

    Res.(ev).all_neuron_resetStrength   = cell(nR,1);
end

%% =========================
% 3) MAIN LOOP
% =========================

for r = 1:nR

    fprintf('\n============================================================\n');
    fprintf('Recording %d / %d\n', r, nR);
    fprintf('%s\n', recordings(r).folder);
    % Sensory events only exist in highact folders
    if ~contains(lower(recordings(r).folder), 'highact')
        fprintf('  Skipping: sensory events only analyzed in highact folders.\n');
        continue;
    end

    S = load(recordings(r).dffPath);
    [dff, fps] = extract_dff_and_fps(S);

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
    recDurMin = nFrames / fps / 60;

    recLabel = make_recording_label(recordings(r).folder);

    fprintf('  DFF: %d neurons x %d frames | fps = %.4f | duration = %.2f min\n', ...
        nNeurons, nFrames, fps, recDurMin);

    % -------- extract sensory events --------
    stim = extract_stim_onsets_by_type(S, nFrames);

    for e = 1:numel(eventTypes)
        ev = eventTypes{e};

        Res.(ev).R.folder{r} = recordings(r).folder;
        Res.(ev).R.recordingLabel{r} = recLabel;

        Res.(ev).R.nNeurons_total(r) = nNeurons;
        Res.(ev).R.fps(r)            = fps;
        Res.(ev).R.nFrames(r)        = nFrames;
        Res.(ev).R.durationMin(r)    = recDurMin;

        if ~isfield(stim, ev) || isempty(stim.(ev))
            fprintf('  [%s] No events found.\n', ev);
            continue;
        end

        eventFrames = unique(round(stim.(ev)(:)));
        eventFrames = eventFrames(isfinite(eventFrames));
        eventFrames = eventFrames(eventFrames >= 1 & eventFrames <= nFrames);

        Res.(ev).R.nEvents_total(r) = numel(eventFrames);

        fprintf('  [%s] Raw events: %d\n', ev, numel(eventFrames));

        % -------- edge and spacing exclusion --------
        periFrames = round(periSec * fps);
        minInterEventFrames = round(minInterEventSec * fps);

        keepEvent = eventFrames > periFrames & eventFrames < (nFrames - periFrames);
        eventFrames = eventFrames(keepEvent);

        if isempty(eventFrames)
            fprintf('  [%s] No events after edge exclusion.\n', ev);
            continue;
        end

        evKeep = true(size(eventFrames));
        for k = 2:numel(eventFrames)
            if eventFrames(k) - eventFrames(k-1) < minInterEventFrames
                evKeep(k) = false;
            end
        end
        eventFrames = eventFrames(evKeep);

        Res.(ev).R.nEvents_used(r) = numel(eventFrames);

        fprintf('  [%s] Used events: %d\n', ev, numel(eventFrames));

        % if numel(eventFrames) < minEventCount
        %     fprintf('  [%s] Too few usable events. Skipping.\n', ev);
        %     continue;
        % end

        % -------- neuron filtering --------
        keepNeuron = mean(isnan(dff), 2) <= maxNanFracPerNeuron;
        dff_use = dff(keepNeuron, :);
        dff_use = dff_use(~all(isnan(dff_use),2), :);

        Res.(ev).R.nNeurons_used(r) = size(dff_use,1);

        if size(dff_use,1) < 3
            fprintf('  [%s] Too few usable neurons. Skipping.\n', ev);
            continue;
        end

        dff_use = fill_nans_rowwise(dff_use);

        % -------- reset metrics --------
        try
            M = compute_event_reset_metrics( ...
                dff_use, eventFrames, fps, ...
                frequencyRange, butterOrder, periSec, ...
                preStatSec, postStatSec, postPeakSec);
        catch ME
            warning('[%s] reset analysis failed in %s\n%s', ev, recordings(r).folder, ME.message);
            continue;
        end

        Res.(ev).R.pre_ITPC_mean(r)         = mean(M.neuron_preITPC, 'omitnan');
        Res.(ev).R.post_ITPC_mean(r)        = mean(M.neuron_postITPC, 'omitnan');
        Res.(ev).R.delta_ITPC_mean(r)       = mean(M.neuron_deltaITPC, 'omitnan');
        Res.(ev).R.resetStrength_mean(r)    = mean(M.neuron_resetStrength, 'omitnan');

        Res.(ev).R.pre_instFreq_mean(r)     = mean(M.neuron_preFreq, 'omitnan');
        Res.(ev).R.post_instFreq_mean(r)    = mean(M.neuron_postFreq, 'omitnan');
        Res.(ev).R.delta_instFreq_mean(r)   = mean(M.neuron_deltaFreq, 'omitnan');

        Res.(ev).R.pre_ITPC_median(r)       = median(M.neuron_preITPC, 'omitnan');
        Res.(ev).R.post_ITPC_median(r)      = median(M.neuron_postITPC, 'omitnan');
        Res.(ev).R.delta_ITPC_median(r)     = median(M.neuron_deltaITPC, 'omitnan');
        Res.(ev).R.resetStrength_median(r)  = median(M.neuron_resetStrength, 'omitnan');

        Res.(ev).R.pre_instFreq_median(r)   = median(M.neuron_preFreq, 'omitnan');
        Res.(ev).R.post_instFreq_median(r)  = median(M.neuron_postFreq, 'omitnan');
        Res.(ev).R.delta_instFreq_median(r) = median(M.neuron_deltaFreq, 'omitnan');

        Res.(ev).all_ITPC_curves{r}           = M.meanITPC(:);
        Res.(ev).all_freq_curves{r}           = M.meanInstFreq(:);
        Res.(ev).all_meanFilt_curves{r}       = M.meanPopulationFilt(:);
        Res.(ev).all_tSec{r}                  = M.tSec(:);

        Res.(ev).all_ITPC_curves_perNeuron{r} = M.itpcCurves;
        Res.(ev).all_freq_curves_perNeuron{r} = M.freqCurves;
        Res.(ev).all_filtCurves_perNeuron{r}  = M.meanFiltCurves;
        Res.(ev).all_neuron_resetStrength{r}  = M.neuron_resetStrength(:);
    end
end

%% =========================
% 4) TABLES + STATS + FIGURES
% =========================

for e = 1:numel(eventTypes)
    ev = eventTypes{e};
    R = Res.(ev).R;

    fprintf('\n================ %s ================\n', upper(ev));

    T = table( ...
        R.recordingLabel, R.folder, ...
        R.nNeurons_total, R.nNeurons_used, ...
        R.fps, R.nFrames, R.durationMin, ...
        R.nEvents_total, R.nEvents_used, ...
        R.pre_ITPC_mean, R.post_ITPC_mean, R.delta_ITPC_mean, R.resetStrength_mean, ...
        R.pre_instFreq_mean, R.post_instFreq_mean, R.delta_instFreq_mean, ...
        R.pre_ITPC_median, R.post_ITPC_median, R.delta_ITPC_median, R.resetStrength_median, ...
        R.pre_instFreq_median, R.post_instFreq_median, R.delta_instFreq_median, ...
        'VariableNames', { ...
        'recordingLabel','folder', ...
        'nNeurons_total','nNeurons_used', ...
        'fps','nFrames','durationMin', ...
        'nEvents_total','nEvents_used', ...
        'pre_ITPC_mean','post_ITPC_mean','delta_ITPC_mean','resetStrength_mean', ...
        'pre_instFreq_mean','post_instFreq_mean','delta_instFreq_mean', ...
        'pre_ITPC_median','post_ITPC_median','delta_ITPC_median','resetStrength_median', ...
        'pre_instFreq_median','post_instFreq_median','delta_instFreq_median'});

    statsSummary = struct();
    statsSummary.ITPC_pre_post      = run_paired_signrank(R.pre_ITPC_median, R.post_ITPC_median);
    statsSummary.instFreq_pre_post  = run_paired_signrank(R.pre_instFreq_median, R.post_instFreq_median);
    statsSummary.deltaITPC_vs_zero      = run_signrank_vs_zero(R.delta_ITPC_median);
    statsSummary.resetStrength_vs_zero  = run_signrank_vs_zero(R.resetStrength_median);
    statsSummary.deltaFreq_vs_zero      = run_signrank_vs_zero(R.delta_instFreq_median);

    fprintf('ITPC pre vs post:       n = %d, p = %.4g\n', ...
        statsSummary.ITPC_pre_post.n, statsSummary.ITPC_pre_post.p);
    fprintf('InstFreq pre vs post:   n = %d, p = %.4g\n', ...
        statsSummary.instFreq_pre_post.n, statsSummary.instFreq_pre_post.p);

    if saveOutputs
        writetable(T, fullfile(outFolder, sprintf('summary_%s_reset_ultraslow.csv', ev)));
        save(fullfile(outFolder, sprintf('summary_%s_reset_ultraslow.mat', ev)), ...
            'T', 'R', 'statsSummary');
    end

    % -------- summary figure --------
    validRec = isfinite(R.pre_ITPC_median) & isfinite(R.post_ITPC_median) & ...
               isfinite(R.pre_instFreq_median) & isfinite(R.post_instFreq_median);

    if any(validRec)
        figR = figure('Name', sprintf('%s-triggered reset', ev), ...
            'Position',[100 100 1400 700]);
        tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact');

        xPair = [1 2];
        idxV  = find(validRec(:))';

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
        set(gca,'XTick',xPair,'XTickLabel',{'Pre','Post'});
        ylabel('Median ITPC');
        title(sprintf('%s ITPC\nn=%d, p=%.3g %s', ev, ...
            statsSummary.ITPC_pre_post.n, statsSummary.ITPC_pre_post.p, p_to_stars(statsSummary.ITPC_pre_post.p)));
        add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.ITPC_pre_post.p);

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
        set(gca,'XTick',xPair,'XTickLabel',{'Pre','Post'});
        ylabel('Median instantaneous frequency (Hz)');
        title(sprintf('%s inst. freq\nn=%d, p=%.3g %s', ev, ...
            statsSummary.instFreq_pre_post.n, statsSummary.instFreq_pre_post.p, p_to_stars(statsSummary.instFreq_pre_post.p)));
        add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.instFreq_pre_post.p);

        nexttile; hold on;
        plot_mean_curve_with_sem(Res.(ev).all_ITPC_curves, Res.(ev).all_tSec, 'b');
        xline(0, '--k');
        xline(preStatSec(1), ':k'); xline(preStatSec(2), ':k');
        xline(postStatSec(1), ':r'); xline(postStatSec(2), ':r');
        xlabel('Time from event onset (s)');
        ylabel('ITPC');
        title(sprintf('%s mean ITPC', ev));
        xlim([-periSec periSec]);

        nexttile; hold on;
        plot_mean_curve_with_sem(Res.(ev).all_freq_curves, Res.(ev).all_tSec, 'r');
        xline(0, '--k');
        xline(preStatSec(1), ':k'); xline(preStatSec(2), ':k');
        xline(postStatSec(1), ':r'); xline(postStatSec(2), ':r');
        xlabel('Time from event onset (s)');
        ylabel('Instantaneous frequency (Hz)');
        title(sprintf('%s mean inst. frequency', ev));
        xlim([-periSec periSec]);

        sgtitle(sprintf('%s-triggered reset of ultraslow neural activity', ev));

        if saveOutputs
            saveas(figR, fullfile(outFolder, sprintf('summary_%s_reset_ultraslow.png', ev)));
        end
    end

    % -------- pooled PSTH-like ITPC --------
    make_pooled_psth_figures( ...
        Res.(ev).all_ITPC_curves_perNeuron, ...
        Res.(ev).all_filtCurves_perNeuron, ...
        Res.(ev).all_tSec, ...
        Res.(ev).all_neuron_resetStrength, ...
        ev, frequencyRange, outFolder, saveOutputs);
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

function stim = extract_stim_onsets_by_type(S, nFrames)
% Extract stimulus onsets from eve_two_pho
%
% Assumptions:
%   - S.eve_two_pho is a numeric vector of event onset frames
%   - first entry should be discarded
%   - last entry should be discarded
%   - remaining events are ordered as:
%         first third   = light
%         second third  = tap
%         last third    = lighttap

    stim = struct('light', [], 'tap', [], 'lighttap', []);

    if ~isfield(S, 'stims') || ~isnumeric(S.stims)
        return;
    end

    ev = S.stims(:);
    ev = ev(isfinite(ev));

    % keep only plausible frame indices
    ev = ev(ev >= 1 & ev <= nFrames);

    % need at least 3 events after trimming
    if numel(ev) < 5
        warning('eve_two_pho has too few entries after filtering.');
        return;
    end

    
    nEv = numel(ev);
    nPerCond = floor(nEv / 3);

    if nPerCond < 1
        warning('Not enough events in eve_two_pho to split into 3 conditions.');
        return;
    end

    % if not perfectly divisible by 3, discard extras at the end
    nUse = 3 * nPerCond;
    if nUse < nEv
        warning('eve_two_pho length after trimming is not divisible by 3; discarding last %d extra event(s).', ...
            nEv - nUse);
        ev = ev(1:nUse);
    end

    stim.light    = ev(1:nPerCond);
    stim.tap      = ev(nPerCond+1 : 2*nPerCond);
    stim.lighttap = ev(2*nPerCond+1 : 3*nPerCond);
end

function x = get_possible_field(S, names)
    x = [];
    for i = 1:numel(names)
        if isfield(S, names{i}) && isnumeric(S.(names{i}))
            x = S.(names{i})(:);
            return;
        end
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

function M = compute_event_reset_metrics(dff, eventFrames, fps, frequencyRange, butterOrder, periSec, preStatSec, postStatSec, postPeakSec)

    [nNeurons, ~] = size(dff);
    periFrames = round(periSec * fps);
    tFrames = (-periFrames:periFrames);
    tSec = tFrames / fps;

    preMask      = tSec >= preStatSec(1)  & tSec <= preStatSec(2);
    postMask     = tSec >= postStatSec(1) & tSec <= postStatSec(2);
    postPeakMask = tSec >= postPeakSec(1) & tSec <= postPeakSec(2);

    Wn = frequencyRange / (fps/2);
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

    nE = numel(eventFrames);

    for n = 1:nNeurons
        x = double(dff(n,:));
        x = filtfilt(b, a, x);

        z = hilbert(x);
        ph = angle(z);

        phUn = unwrap(ph);
        instFreq = [nan diff(phUn)] * fps / (2*pi);
        instFreq = smoothdata(instFreq, 'gaussian', max(3, round(fps)));

        phaseMat = nan(nE, numel(tSec));
        freqMat  = nan(nE, numel(tSec));
        filtMat  = nan(nE, numel(tSec));

        for k = 1:nE
            idx = eventFrames(k) + tFrames;
            phaseMat(k,:) = ph(idx);
            freqMat(k,:)  = instFreq(idx);
            filtMat(k,:)  = x(idx);
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
    M.tSec               = tSec(:);
    M.itpcCurves         = itpcCurves;
    M.freqCurves         = freqCurves;
    M.meanFiltCurves     = meanFiltCurves;
    M.meanITPC           = mean(itpcCurves, 1, 'omitnan').';
    M.meanInstFreq       = mean(freqCurves, 1, 'omitnan').';
    M.meanPopulationFilt = mean(meanFiltCurves, 1, 'omitnan').';

    M.neuron_preITPC       = neuron_preITPC;
    M.neuron_postITPC      = neuron_postITPC;
    M.neuron_deltaITPC     = neuron_deltaITPC;
    M.neuron_resetStrength = neuron_resetStrength;

    M.neuron_preFreq       = neuron_preFreq;
    M.neuron_postFreq      = neuron_postFreq;
    M.neuron_deltaFreq     = neuron_deltaFreq;
end

function plot_mean_curve_with_sem(curves, tSecCells, lineColor)

    validCurves = ~cellfun(@isempty, curves) & ~cellfun(@isempty, tSecCells);
    if ~any(validCurves)
        text(0.1,0.5,'No valid curves','FontSize',12);
        axis off;
        return;
    end

    refIdx = find(validCurves, 1, 'first');
    tSec = tSecCells{refIdx}(:);

    keepIdx = [];
    for i = find(validCurves(:))'
        if numel(curves{i}) == numel(tSec)
            keepIdx(end+1) = i; %#ok<SAGROW>
        end
    end

    X = nan(numel(tSec), numel(keepIdx));
    for k = 1:numel(keepIdx)
        X(:,k) = curves{keepIdx(k)}(:);
    end

    m = mean(X, 2, 'omitnan');
    s = std(X, 0, 2, 'omitnan') ./ sqrt(size(X,2));

    fill([tSec; flipud(tSec)], [m-s; flipud(m+s)], ...
        lineColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(tSec, m, lineColor, 'LineWidth', 1.8);
end

function make_pooled_psth_figures(itpcCells, actCells, tCells, resetCells, ev, frequencyRange, outFolder, saveOutputs)

    if ~any(~cellfun(@isempty, tCells))
        return;
    end

    refIdx  = find(~cellfun(@isempty, tCells), 1, 'first');
    tCommon = tCells{refIdx}(:)';

    allITPC  = [];
    allAct   = [];
    allReset = [];

    for r = 1:numel(tCells)
        if isempty(tCells{r}) || isempty(resetCells{r})
            continue;
        end

        tThis = tCells{r}(:)';
        rs    = resetCells{r}(:);

        if ~isempty(itpcCells{r})
            X1 = itpcCells{r};
            if size(X1,2) ~= numel(tCommon) || numel(tThis) ~= numel(tCommon) || any(abs(tThis - tCommon) > 1e-9)
                Xrs = nan(size(X1,1), numel(tCommon));
                for n = 1:size(X1,1)
                    Xrs(n,:) = interp1(tThis, X1(n,:), tCommon, 'linear', NaN);
                end
                X1 = Xrs;
            end
            allITPC = [allITPC; X1]; %#ok<AGROW>
        end

        if ~isempty(actCells{r})
            X2 = actCells{r};
            if size(X2,2) ~= numel(tCommon) || numel(tThis) ~= numel(tCommon) || any(abs(tThis - tCommon) > 1e-9)
                Xrs = nan(size(X2,1), numel(tCommon));
                for n = 1:size(X2,1)
                    Xrs(n,:) = interp1(tThis, X2(n,:), tCommon, 'linear', NaN);
                end
                X2 = Xrs;
            end
            allAct = [allAct; X2]; %#ok<AGROW>
        end

        allReset = [allReset; rs]; %#ok<AGROW>
    end

    if ~isempty(allITPC)
        fig1 = plot_itpc_psth_like(allITPC, tCommon, allReset, ...
            sprintf('%s | pooled ITPC %.02f-%.02f Hz', ev, frequencyRange(1), frequencyRange(2)));
        if saveOutputs
            saveas(fig1, fullfile(outFolder, sprintf('pooled_%s_itpc_psth_like.png', ev)));
        end
    end

    if ~isempty(allAct)
        fig2 = plot_activity_psth_like(allAct, tCommon, allReset, ...
            sprintf('%s | pooled bandpassed activity %.02f-%.02f Hz', ev, frequencyRange(1), frequencyRange(2)), ...
            'Bandpassed \DeltaF/F');
        if saveOutputs
            saveas(fig2, fullfile(outFolder, sprintf('pooled_%s_activity_psth_like.png', ev)));
        end
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
    title('Sensory-triggered phase concentration');
    xlim([min(tSec) max(tSec)]);
    box off;

    nexttile;
    imagesc(tSec, 1:size(X,1), X);
    set(gca, 'YDir', 'normal');
    xline(0, '--k', 'LineWidth', 1);

    xlabel('Time from event onset (s)');
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
    title('Sensory-triggered neural activity');
    xlim([min(tSec) max(tSec)]);
    box off;

    nexttile;
    imagesc(tSec, 1:size(X,1), X);
    set(gca, 'YDir', 'normal');
    xline(0, '--k', 'LineWidth', 1);

    xlabel('Time from event onset (s)');
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