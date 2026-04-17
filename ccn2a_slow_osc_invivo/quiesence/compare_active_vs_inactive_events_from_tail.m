%% =========================================================================
% compare_active_vs_inactive_events_from_tail.m
%
% Tail-defined least-active vs most-active windows:
% compare neural strong-event statistics from OASIS-deconvolved activity
%
% Requires:
%   - deconvolveCa.m on MATLAB path
%   - tail_quick.mat
%   - dffs_repact_respcells.mat
%
% Optional:
%   - existing spks_oasis.mat in each suite2p folder
%
% Output:
%   - per-recording summary table
%   - paired signrank statistics
%   - one summary figure
% =========================================================================

clear; clc; close all;

%% =========================
% 0) USER SETTINGS
% =========================

rootFolder   = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low';
tailFileName = 'tail_quick.mat';
dffFileName  = 'dffs_repact_respcells.mat';
spkFileName  = 'spks_oasis.mat';

% Tail-defined windows
winMin  = 20;   % minutes
stepMin = 5;    % minutes

% NaN handling
maxNanFracPerNeuron = 0.20;

% OASIS options
opt = struct;
opt.type                   = 'ar1';
opt.method                 = 'constrained';
opt.sn                     = [];
opt.b                      = 0;
opt.optimize_b             = true;
opt.optimize_pars          = true;
opt.smin                   = -2;
opt.thresh_factor          = 1.0;
opt.remove_large_residuals = false;

% Strong-event detection parameters
base_factor   = 4;
min_dist_sec  = 3.5;
prom_factor   = 2.5;
peak_quantile = 0.9;
min_norm_amp  = 10;

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_tail_events_results');

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

recordings = struct('folder', {}, 'tailPath', {}, 'dffPath', {}, 'spkPath', {});
for i = 1:numel(tailFiles)
    thisFolder = tailFiles(i).folder;
    if isKey(dffMap, thisFolder)
        recordings(end+1).folder  = thisFolder; %#ok<SAGROW>
        recordings(end).tailPath  = fullfile(tailFiles(i).folder, tailFiles(i).name);
        recordings(end).dffPath   = dffMap(thisFolder);
        recordings(end).spkPath   = fullfile(thisFolder, spkFileName);
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

% Event metrics
R.low_eventRate_mean      = nan(nR,1);   % mean over neurons (events/s)
R.high_eventRate_mean     = nan(nR,1);

R.low_eventRate_median    = nan(nR,1);
R.high_eventRate_median   = nan(nR,1);

R.low_totalEventCount     = nan(nR,1);
R.high_totalEventCount    = nan(nR,1);

R.low_fracActiveNeurons   = nan(nR,1);
R.high_fracActiveNeurons  = nan(nR,1);

R.low_eventAmp_mean       = nan(nR,1);   % mean z-amplitude of strong events
R.high_eventAmp_mean      = nan(nR,1);

allLowRates               = cell(nR,1);
allHighRates              = cell(nR,1);
allLowEventAmp            = cell(nR,1);
allHighEventAmp           = cell(nR,1);

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
    % 3.4 Select usable neurons
    % ---------------------------------------------------------------------
    dff_low  = dff(:, low_startFrame:low_endFrame);
    dff_high = dff(:, high_startFrame:high_endFrame);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;

    keepUse = keepLow & keepHigh & ~all(isnan(dff),2);

    dff_use = dff(keepUse,:);
    R.nNeurons_used_low(r)  = sum(keepLow);
    R.nNeurons_used_high(r) = sum(keepHigh);

    if size(dff_use,1) < 1
        warning('No shared usable neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    dff_use = fill_nans_rowwise(dff_use);
    [Nuse, Tuse] = size(dff_use); %#ok<NASGU>

    % 3.5 Load or compute OASIS spikes
    spks = [];
    
    % dff_use is N x T
    Ncur = size(dff_use, 1);
    Tcur = size(dff_use, 2);
    
    % ---------------------------------------------------------------------
    % Try loading existing spikes
    % Expected convention for this script: spks = T x N
    % ---------------------------------------------------------------------
    if isfile(recordings(r).spkPath)
        L = load(recordings(r).spkPath);
    
        if isfield(L, 'spks')
            spks = L.spks;
    
            fprintf('  Loaded existing spks_oasis.mat with size [%d x %d]\n', ...
                size(spks,1), size(spks,2));
    
            % Good case: already T x N
            if size(spks,1) == Tcur && size(spks,2) == Ncur
                % keep as is
    
            % Old file saved as N x T -> transpose it
            elseif size(spks,1) == Ncur && size(spks,2) == Tcur
                spks = spks.';
                fprintf('  Transposed existing spks to [%d x %d] (T x N)\n', ...
                    size(spks,1), size(spks,2));
    
            % Vectorized output -> reshape back to T x N
            elseif isvector(spks) && numel(spks) == Tcur * Ncur
                spks = reshape(spks, Tcur, Ncur);
                fprintf('  Reshaped vectorized existing spks to [%d x %d] (T x N)\n', ...
                    size(spks,1), size(spks,2));
    
            % Anything else is unusable, recompute
            else
                warning('Existing spks_oasis size mismatch in %s. Recomputing.', recordings(r).spkPath);
                spks = [];
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % Recompute with OASIS if needed
    % Input to deconvolveCa is T x N
    % ---------------------------------------------------------------------
    if isempty(spks)
        fprintf('  Running OASIS deconvolution...\n');
    
        try
            [~, s_est, opts_used] = deconvolveCa(dff_use.', opt);   % input is T x N
    
            % Force output back to T x N
            if isvector(s_est)
                if numel(s_est) ~= Tcur * Ncur
                    error('OASIS returned vector of length %d, expected %d (= T*N).', ...
                        numel(s_est), Tcur * Ncur);
                end
                spks = reshape(s_est, Tcur, Ncur);
    
            elseif ismatrix(s_est)
                if size(s_est,1) == Tcur && size(s_est,2) == Ncur
                    spks = s_est;
    
                elseif size(s_est,1) == Ncur && size(s_est,2) == Tcur
                    spks = s_est.';
    
                elseif numel(s_est) == Tcur * Ncur
                    spks = reshape(s_est, Tcur, Ncur);
    
                else
                    error('Unexpected OASIS output size: [%d x %d], expected [%d x %d] or [%d x %d].', ...
                        size(s_est,1), size(s_est,2), Tcur, Ncur, Ncur, Tcur);
                end
    
            else
                error('Unexpected OASIS output type for s_est.');
            end
    
            fprintf('  OASIS output reshaped to spks = [%d x %d] (T x N)\n', ...
                size(spks,1), size(spks,2));
    
            if saveOutputs
                save(recordings(r).spkPath, 'spks', 'opts_used', '-v7.3');
            end
    
        catch ME
            warning('OASIS failed in %s\n%s', recordings(r).folder, ME.message);
            continue;
        end
    end
    
    % ---------------------------------------------------------------------
    % Final safety check
    % ---------------------------------------------------------------------
    if size(spks,1) ~= Tcur || size(spks,2) ~= Ncur
        warning('Final spks size mismatch in %s: got [%d x %d], expected [%d x %d]. Skipping.', ...
            recordings(r).folder, size(spks,1), size(spks,2), Tcur, Ncur);
        continue;
    end

    % ---------------------------------------------------------------------
    % 3.6 Strong-event detection
    % ---------------------------------------------------------------------
    fprintf('  Detecting strong events...\n');
    [eventsMat, eventsStruct, z_amp] = strong_events_mad( ...
        spks, fps, base_factor, min_dist_sec, prom_factor, peak_quantile, min_norm_amp); %#ok<ASGLU>
    
    fprintf('  size(dff_use)  = [%d x %d]\n', size(dff_use,1), size(dff_use,2));
    fprintf('  size(spks)     = [%d x %d]\n', size(spks,1), size(spks,2));
    fprintf('  size(eventsMat)= [%d x %d]\n', size(eventsMat,1), size(eventsMat,2));
    
    % If eventsMat came out as T x N, transpose it
    if size(eventsMat,1) == size(dff_use,2) && size(eventsMat,2) == size(dff_use,1)
        eventsMat = eventsMat.';
        fprintf('  Transposed eventsMat to match neurons x frames.\n');
    end
    
    % z_amp is expected T x N in summarize_events_in_window
    % keep it that way unless your helper returned N x T
    if size(z_amp,1) == size(dff_use,1) && size(z_amp,2) == size(dff_use,2)
        z_amp = z_amp.';   % convert N x T -> T x N
        fprintf('  Transposed z_amp to match frames x neurons.\n');
    end
    
    if size(eventsMat,1) ~= size(dff_use,1) || size(eventsMat,2) ~= size(dff_use,2)
        warning('Event matrix size mismatch in %s. Skipping.', recordings(r).folder);
        continue;
    end
    
    if size(z_amp,1) ~= size(dff_use,2) || size(z_amp,2) ~= size(dff_use,1)
        warning('z_amp size mismatch in %s. Skipping.', recordings(r).folder);
        continue;
    end
        
    % z_amp is expected T x N in summarize_events_in_window
    % keep it that way unless your helper returned N x T
    if size(z_amp,1) == size(dff_use,1) && size(z_amp,2) == size(dff_use,2)
        z_amp = z_amp.';   % convert N x T -> T x N
        fprintf('  Transposed z_amp to match frames x neurons.\n');
    end
    
    if size(eventsMat,1) ~= size(dff_use,1) || size(eventsMat,2) ~= size(dff_use,2)
        warning('Event matrix size mismatch in %s. Skipping.', recordings(r).folder);
        continue;
    end
    
    if size(z_amp,1) ~= size(dff_use,2) || size(z_amp,2) ~= size(dff_use,1)
        warning('z_amp size mismatch in %s. Skipping.', recordings(r).folder);
        continue;
    end

    % ---------------------------------------------------------------------
    % 3.7 Quantify events in low/high windows
    % ---------------------------------------------------------------------
    low_idx  = low_startFrame:low_endFrame;
    high_idx = high_startFrame:high_endFrame;

    [lowMetrics, lowAmp]   = summarize_events_in_window(eventsMat, z_amp, low_idx, fps);
    [highMetrics, highAmp] = summarize_events_in_window(eventsMat, z_amp, high_idx, fps);

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

    allLowRates{r}      = lowMetrics.rates_per_neuron;
    allHighRates{r}     = highMetrics.rates_per_neuron;
    allLowEventAmp{r}   = lowAmp;
    allHighEventAmp{r}  = highAmp;

    fprintf('  Low  events: total=%d, mean rate=%.4g, frac active=%.3f\n', ...
        lowMetrics.total_count, lowMetrics.rate_mean, lowMetrics.frac_active_neurons);
    fprintf('  High events: total=%d, mean rate=%.4g, frac active=%.3f\n', ...
        highMetrics.total_count, highMetrics.rate_mean, highMetrics.frac_active_neurons);
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used_low, R.nNeurons_used_high, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.low_startFrame, R.low_endFrame, R.high_startFrame, R.high_endFrame, ...
    R.low_boutCount, R.high_boutCount, R.low_boutRate_perMin, R.high_boutRate_perMin, ...
    R.low_eventRate_mean, R.high_eventRate_mean, ...
    R.low_eventRate_median, R.high_eventRate_median, ...
    R.low_totalEventCount, R.high_totalEventCount, ...
    R.low_fracActiveNeurons, R.high_fracActiveNeurons, ...
    R.low_eventAmp_mean, R.high_eventAmp_mean, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used_low','nNeurons_used_high', ...
    'fps','nFrames','durationMin', ...
    'low_startFrame','low_endFrame','high_startFrame','high_endFrame', ...
    'low_boutCount','high_boutCount','low_boutRate_perMin','high_boutRate_perMin', ...
    'low_eventRate_mean','high_eventRate_mean', ...
    'low_eventRate_median','high_eventRate_median', ...
    'low_totalEventCount','high_totalEventCount', ...
    'low_fracActiveNeurons','high_fracActiveNeurons', ...
    'low_eventAmp_mean','high_eventAmp_mean'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_events_active_vs_inactive.csv'));
    save(fullfile(outFolder, 'summary_events_active_vs_inactive.mat'), ...
        'R', 'T', 'allLowRates', 'allHighRates', 'allLowEventAmp', 'allHighEventAmp');
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

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_events.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURE
% =========================

validRec = isfinite(R.low_eventRate_mean) & isfinite(R.high_eventRate_mean) & ...
           isfinite(R.low_totalEventCount) & isfinite(R.high_totalEventCount) & ...
           isfinite(R.low_fracActiveNeurons) & isfinite(R.high_fracActiveNeurons) & ...
           isfinite(R.low_eventAmp_mean) & isfinite(R.high_eventAmp_mean);

if any(validRec)

    xPair = [1 2];
    idxV  = find(validRec(:))';

    figE = figure('Name','Active vs inactive event comparison', ...
        'Position',[100 100 1400 700]);
    tiledlayout(1,4, 'Padding','compact', 'TileSpacing','compact');

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
    ylabel('Total strong-event count');
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
    ylabel('Mean strong-event amplitude (z)');
    title(sprintf('Event amplitude\nn=%d, p=%.3g %s', ...
        statsSummary.ampMean.n, statsSummary.ampMean.p, p_to_stars(statsSummary.ampMean.p)));
    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.ampMean.p);

    sgtitle('Strong-event comparison: tail-defined inactive vs active windows');

    if saveOutputs
        saveas(figE, fullfile(outFolder, 'summary_events_active_vs_inactive.png'));
    end
else
    warning('No valid recordings available for event summary figure.');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% ================== HELPERS ==================
function [eventsMat, eventsStruct, z_amp] = strong_events_mad(spks, fs, base_factor, min_dist_sec, prom_factor, peak_quantile, min_norm_amp)
% spks: T x N
% returns:
%   eventsMat   : N x T logical
%   eventsStruct: 1 x N cell, each with fields idx, amp
%   z_amp       : T x N z-scored spike amplitudes

    [T, N] = size(spks);
    eventsMat = false(N, T);
    eventsStruct = cell(N,1);
    z_amp = nan(T,N);

    minDist = max(1, round(min_dist_sec * fs));

    for n = 1:N
        x = double(spks(:,n));
        x(~isfinite(x)) = 0;

        medx = median(x, 'omitnan');
        sigx = mad(x,1);
        if ~isfinite(sigx) || sigx <= 0
            sigx = std(x, 0, 'omitnan');
        end
        if ~isfinite(sigx) || sigx <= 0
            sigx = 1;
        end

        z = (x - medx) ./ sigx;
        z_amp(:,n) = z;

        minH = base_factor;
        minP = prom_factor;

        [pks, locs, ~, prom] = findpeaks(z, ...
            'MinPeakHeight', minH, ...
            'MinPeakDistance', minDist, ...
            'MinPeakProminence', minP);

        if isempty(locs)
            eventsStruct{n} = struct('idx', [], 'amp', [], 'prom', []);
            continue;
        end

        % Keep only stronger peaks for each neuron
        qThr = quantile(pks, peak_quantile);
        keep = pks >= max(qThr, min_norm_amp);

        locs = locs(keep);
        pks  = pks(keep);
        prom = prom(keep);

        eventsMat(n, locs) = true;
        eventsStruct{n} = struct('idx', locs(:), 'amp', pks(:), 'prom', prom(:));
    end
end

function [M, ampVals] = summarize_events_in_window(eventsMat, z_amp, idx, fps)
% eventsMat: N x T logical
% z_amp    : T x N
% idx      : frame indices for one window

    idx = idx(:)';
    idx = idx(idx >= 1 & idx <= size(eventsMat,2));

    evWin = eventsMat(:, idx);   % N x Tw
    durSec = numel(idx) / fps;

    counts = sum(evWin, 2);      % N x 1
    rates  = counts / durSec;    % events/s

    fracActive = mean(counts > 0, 'omitnan');
    totalCount = sum(counts, 'omitnan');

    % collect amplitudes of events in this window
    ampVals = [];
    [rows, cols] = find(evWin);
    if ~isempty(rows)
        linIdx = sub2ind(size(z_amp), idx(cols(:))', rows(:)); % wrong orientation if used directly
        % safer:
        ampVals = arrayfun(@(k) z_amp(idx(cols(k)), rows(k)), 1:numel(rows)).';
    end

    M = struct();
    M.counts_per_neuron   = counts;
    M.rates_per_neuron    = rates;
    M.rate_mean           = mean(rates, 'omitnan');
    M.rate_median         = median(rates, 'omitnan');
    M.total_count         = totalCount;
    M.frac_active_neurons = fracActive;
    M.amp_mean            = mean(ampVals, 'omitnan');
end

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