%% =========================================================================
% compare_active_vs_inactive_corr_vs_distance_from_tail.m
%
% Tail-defined least-active vs most-active windows:
% compare neural correlation-vs-distance structure in the least active and
% most active periods of each recording.
%
% Requires:
%   - tail_quick.mat
%   - dffs_repact_respcells.mat
%   - corr_vs_distance_AO.m on MATLAB path
%   - Statistics Toolbox
%
% Output:
%   - per-recording summary table
%   - paired signrank statistics
%   - summary figures
%
% Notes:
%   - corr_vs_distance_AO is applied separately to the least-active and
%     most-active windows.
%   - hemispheres are averaged exactly like in your snippet:
%         band_curve = mean(cd_band, 1, 'omitnan')
%   - distance bins are defined per recording from XY positions of the
%     usable neurons.
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

% corr_vs_distance settings
please_split = 1;
side2split   = 2;     % split by Y
nDistBins    = 50;    % same as your snippet
minMaxDistFallback = 100;

% Optional randomized control
doRandomControl = true;
nRand = 100;

% Summary-bin settings for scalar metrics
nNearBins = 5;        % average first few bins as "short-range correlation"

% Save outputs
saveOutputs = true;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_corr_vs_distance_results');

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

% Scalar corr-distance summaries
R.low_corr_near_mean      = nan(nR,1);
R.high_corr_near_mean     = nan(nR,1);

R.low_corr_all_mean       = nan(nR,1);
R.high_corr_all_mean      = nan(nR,1);

R.low_corr_auc            = nan(nR,1);
R.high_corr_auc           = nan(nR,1);

R.maxDist                 = nan(nR,1);

% Curves
allBinCenters             = cell(nR,1);
allLowCurve               = cell(nR,1);
allHighCurve              = cell(nR,1);
allLowCD                  = cell(nR,1);
allHighCD                 = cell(nR,1);

% Optional randomized control
allCtrlLow                = cell(nR,1);
allCtrlHigh               = cell(nR,1);

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
    % 3.1 Load DFF, FPS, positions
    % ---------------------------------------------------------------------
    S_dff = load(recordings(r).dffPath);
    [dff, fps] = extract_dff_and_fps(S_dff);
    positions = extract_positions(S_dff);

    if isempty(dff)
        warning('Could not find DFF matrix in %s. Skipping.', recordings(r).dffPath);
        continue;
    end
    if isempty(fps) || ~isscalar(fps) || isnan(fps) || fps <= 0
        warning('Could not find valid fps in %s. Skipping.', recordings(r).dffPath);
        continue;
    end
    if isempty(positions) || size(positions,2) < 2
        warning('Could not find valid positions in %s. Skipping.', recordings(r).dffPath);
        continue;
    end

    dff = double(dff);
    if size(dff,1) > size(dff,2)
        dff = dff.';
    end

    nNeurons = size(dff,1);
    nFrames  = size(dff,2);

    % Try to match position rows to neuron rows
    positions = double(positions);
    if size(positions,1) ~= nNeurons
        if size(positions,2) == nNeurons && size(positions,1) >= 2
            positions = positions.';
        end
    end

    if size(positions,1) ~= nNeurons
        warning('Positions size mismatch in %s. dff has %d neurons, positions has %d rows. Skipping.', ...
            recordings(r).folder, nNeurons, size(positions,1));
        continue;
    end

    R.nNeurons_total(r) = nNeurons;
    R.fps(r)            = fps;
    R.nFrames(r)        = nFrames;
    R.durationMin(r)    = nFrames / fps / 60;

    fprintf('  DFF: %d neurons x %d frames | fps = %.4f | duration = %.2f min\n', ...
        nNeurons, nFrames, fps, R.durationMin(r));
    fprintf('  Positions: %d neurons x %d dims\n', size(positions,1), size(positions,2));

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
    % 3.4 Select usable neurons shared across both windows
    % ---------------------------------------------------------------------
    dff_low  = dff(:, low_startFrame:low_endFrame);
    dff_high = dff(:, high_startFrame:high_endFrame);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;

    keepUse = keepLow & keepHigh & ~all(isnan(dff),2) & ...
              all(isfinite(positions(:,1:2)),2);

    dff_use = dff(keepUse,:);
    pos_use = positions(keepUse,:);

    R.nNeurons_used_low(r)    = sum(keepLow);
    R.nNeurons_used_high(r)   = sum(keepHigh);
    R.nNeurons_used_shared(r) = size(dff_use,1);

    if size(dff_use,1) < 5
        warning('Too few shared usable neurons in %s. Skipping.', recordings(r).folder);
        continue;
    end

    dff_use = fill_nans_rowwise(dff_use);

    % ---------------------------------------------------------------------
    % 3.5 Build shared distance bins for this recording
    % ---------------------------------------------------------------------
    posXY = pos_use(:,1:2);
    pd = squareform(pdist(posXY));
    maxDist = max(pd(:), [], 'omitnan');

    if ~isfinite(maxDist) || maxDist < 50
        maxDist = minMaxDistFallback;
    end

    steps_used = floor(linspace(0, maxDist, nDistBins+1));
    steps_used = unique(steps_used);

    if numel(steps_used) < 3
        steps_used = round(linspace(0, maxDist, nDistBins+1));
        steps_used = unique(steps_used);
    end

    if numel(steps_used) < 3
        warning('Could not build valid distance bins in %s. Skipping.', recordings(r).folder);
        continue;
    end

    binCenters = steps_used(2:end) - diff(steps_used)/2;
    R.maxDist(r) = maxDist;

    % ---------------------------------------------------------------------
    % 3.6 Compute corr vs distance in low and high windows
    % ---------------------------------------------------------------------
    dff_low_use  = dff_use(:, low_startFrame:low_endFrame);
    dff_high_use = dff_use(:, high_startFrame:high_endFrame);

    try
        [cd_low, ~, ~, ~] = corr_vs_distance_AO(dff_low_use, pos_use, ...
            please_split, side2split, [], steps_used);

        [cd_high, ~, ~, ~] = corr_vs_distance_AO(dff_high_use, pos_use, ...
            please_split, side2split, [], steps_used);
    catch ME
        warning('corr_vs_distance_AO failed in %s\n%s', recordings(r).folder, ME.message);
        continue;
    end

    low_curve  = mean(cd_low,  1, 'omitnan');
    high_curve = mean(cd_high, 1, 'omitnan');

    if numel(low_curve) ~= numel(binCenters) || numel(high_curve) ~= numel(binCenters)
        warning('corr_vs_distance curve length mismatch in %s. Skipping.', recordings(r).folder);
        continue;
    end

    allBinCenters{r} = binCenters;
    allLowCurve{r}   = low_curve;
    allHighCurve{r}  = high_curve;
    allLowCD{r}      = cd_low;
    allHighCD{r}     = cd_high;

    % ---------------------------------------------------------------------
    % 3.7 Optional randomized control
    % ---------------------------------------------------------------------
    if doRandomControl
        Ncells = size(dff_use,1);
        nBand  = Ncells;

        ctrl_low = nan(nRand, numel(low_curve));
        ctrl_high = nan(nRand, numel(high_curve));

        for rr = 1:nRand
            ridx = randperm(Ncells, nBand);

            [cd_r_low, ~, ~, ~] = corr_vs_distance_AO(dff_low_use(ridx,:), pos_use, ...
                please_split, side2split, [], steps_used);

            [cd_r_high, ~, ~, ~] = corr_vs_distance_AO(dff_high_use(ridx,:), pos_use, ...
                please_split, side2split, [], steps_used);

            ctrl_low(rr,:)  = mean(cd_r_low,  1, 'omitnan');
            ctrl_high(rr,:) = mean(cd_r_high, 1, 'omitnan');
        end

        allCtrlLow{r}  = ctrl_low;
        allCtrlHigh{r} = ctrl_high;
    end

    % ---------------------------------------------------------------------
    % 3.8 Scalar summaries (distance-range based)
    % ---------------------------------------------------------------------
    nearRange = [19 40];
    allRange  = [19 120];

    idxNear = isfinite(binCenters) & binCenters >= nearRange(1) & binCenters <= nearRange(2);
    idxAll  = isfinite(binCenters) & binCenters >= allRange(1)  & binCenters <= allRange(2);

    % Mean correlation in 19-40 um
    if any(idxNear)
        R.low_corr_near_mean(r)  = mean(low_curve(idxNear),  'omitnan');
        R.high_corr_near_mean(r) = mean(high_curve(idxNear), 'omitnan');
    else
        R.low_corr_near_mean(r)  = NaN;
        R.high_corr_near_mean(r) = NaN;
    end

    % Mean correlation in 19-120 um
    if any(idxAll)
        R.low_corr_all_mean(r)   = mean(low_curve(idxAll),  'omitnan');
        R.high_corr_all_mean(r)  = mean(high_curve(idxAll), 'omitnan');
    else
        R.low_corr_all_mean(r)   = NaN;
        R.high_corr_all_mean(r)  = NaN;
    end

    % AUC in 19-120 um
    idxLowAUC  = idxAll & isfinite(low_curve);
    idxHighAUC = idxAll & isfinite(high_curve);

    if sum(idxLowAUC) >= 2
        R.low_corr_auc(r) = trapz(binCenters(idxLowAUC), low_curve(idxLowAUC));
    else
        R.low_corr_auc(r) = NaN;
    end

    if sum(idxHighAUC) >= 2
        R.high_corr_auc(r) = trapz(binCenters(idxHighAUC), high_curve(idxHighAUC));
    else
        R.high_corr_auc(r) = NaN;
    end

    fprintf('  Low  curve: near-mean=%.4g | all-mean=%.4g | AUC=%.4g\n', ...
        R.low_corr_near_mean(r), R.low_corr_all_mean(r), R.low_corr_auc(r));
    fprintf('  High curve: near-mean=%.4g | all-mean=%.4g | AUC=%.4g\n', ...
        R.high_corr_near_mean(r), R.high_corr_all_mean(r), R.high_corr_auc(r));

    % ---------------------------------------------------------------------
    % 3.9 Per-recording QC figure
    % ---------------------------------------------------------------------
    try
        figQC = figure('Name', ['corr_vs_distance_QC_' R.recordingLabel{r}], ...
            'Position', [100 100 1200 500], 'Color', 'w');

        tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

        % left: window traces
        nexttile; hold on;
        plot((low_startFrame:low_endFrame)/fps, mean(dff_low_use,1,'omitnan'), 'b-', 'LineWidth', 1);
        plot((high_startFrame:high_endFrame)/fps, mean(dff_high_use,1,'omitnan'), 'r-', 'LineWidth', 1);
        xlabel('Time (s)');
        ylabel('Mean dF/F');
        title(sprintf('%s\nLeast vs most active windows', R.recordingLabel{r}), 'Interpreter','none');
        legend({'Least active','Most active'}, 'Location','best');
        box off;

        % right: corr vs distance
        nexttile; hold on;
        plot(binCenters, low_curve,  'b-', 'LineWidth', 2);
        plot(binCenters, high_curve, 'r-', 'LineWidth', 2);
        xlabel('Distance');
        ylabel('Mean correlation');
        title('corr vs distance');
        legend({'Least active','Most active'}, 'Location','best');
        box off;

        if saveOutputs
            saveas(figQC, fullfile(outFolder, ['QC_corr_vs_distance_' sanitize_filename(R.recordingLabel{r}) '.png']));
        end
        close(figQC);
    catch
    end
end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, R.folder, ...
    R.nNeurons_total, R.nNeurons_used_low, R.nNeurons_used_high, R.nNeurons_used_shared, ...
    R.fps, R.nFrames, R.durationMin, ...
    R.low_startFrame, R.low_endFrame, R.high_startFrame, R.high_endFrame, ...
    R.low_boutCount, R.high_boutCount, R.low_boutRate_perMin, R.high_boutRate_perMin, ...
    R.maxDist, ...
    R.low_corr_near_mean, R.high_corr_near_mean, ...
    R.low_corr_all_mean, R.high_corr_all_mean, ...
    R.low_corr_auc, R.high_corr_auc, ...
    'VariableNames', { ...
    'recordingLabel','folder', ...
    'nNeurons_total','nNeurons_used_low','nNeurons_used_high','nNeurons_used_shared', ...
    'fps','nFrames','durationMin', ...
    'low_startFrame','low_endFrame','high_startFrame','high_endFrame', ...
    'low_boutCount','high_boutCount','low_boutRate_perMin','high_boutRate_perMin', ...
    'maxDist', ...
    'low_corr_near_mean','high_corr_near_mean', ...
    'low_corr_all_mean','high_corr_all_mean', ...
    'low_corr_auc','high_corr_auc'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_corr_vs_distance_active_vs_inactive.csv'));
    save(fullfile(outFolder, 'summary_corr_vs_distance_active_vs_inactive.mat'), ...
        'R', 'T', 'allBinCenters', 'allLowCurve', 'allHighCurve', 'allLowCD', 'allHighCD', ...
        'allCtrlLow', 'allCtrlHigh', '-v7.3');
end

%% =========================
% 5) PAIRED STATISTICS
% =========================

statsSummary = struct();
statsSummary.nearMean = run_paired_signrank(R.low_corr_near_mean, R.high_corr_near_mean);
statsSummary.allMean  = run_paired_signrank(R.low_corr_all_mean,  R.high_corr_all_mean);
statsSummary.auc      = run_paired_signrank(R.low_corr_auc,       R.high_corr_auc);

disp('================ PAIRED SIGNRANK TESTS ================')
fprintf('Near-bin mean corr: n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.nearMean.n, statsSummary.nearMean.p, statsSummary.nearMean.deltaMedian);
fprintf('All-bin mean corr:  n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.allMean.n, statsSummary.allMean.p, statsSummary.allMean.deltaMedian);
fprintf('Curve AUC:          n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.auc.n, statsSummary.auc.p, statsSummary.auc.deltaMedian);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary_corr_vs_distance.mat'), 'statsSummary');
end

%% =========================
% 6) SUMMARY FIGURES
% =========================

validRecCurves = ~cellfun(@isempty, allLowCurve) & ~cellfun(@isempty, allHighCurve);

validRecScalars = validRecCurves & ...
                  isfinite(R.low_corr_near_mean) & isfinite(R.high_corr_near_mean) & ...
                  isfinite(R.low_corr_all_mean)  & isfinite(R.high_corr_all_mean);

validRecAUC = validRecScalars & ...
              isfinite(R.low_corr_auc) & isfinite(R.high_corr_auc);

idxCurve  = find(validRecCurves(:))';
idxScalar = find(validRecScalars(:))';
idxAUC    = find(validRecAUC(:))';

% -------------------------------------------------------------------------
% 6.1 Shaded summary curves (real distance axis, truncated to 150)
% -------------------------------------------------------------------------
if any(validRecCurves)

    % common real distance axis
    xMaxPlot = 150;
    xCommon = linspace(0, xMaxPlot, 151);

    lowInterp = nan(numel(idxCurve), numel(xCommon));
    highInterp = nan(numel(idxCurve), numel(xCommon));

    if doRandomControl
        ctrlLowInterp = nan(numel(idxCurve)*nRand, numel(xCommon));
        ctrlHighInterp = nan(numel(idxCurve)*nRand, numel(xCommon));
    else
        ctrlLowInterp = [];
        ctrlHighInterp = [];
    end

    cLow = 1;
    cHigh = 1;

    for k = 1:numel(idxCurve)
        i = idxCurve(k);

        xL = allBinCenters{i};
        yL = allLowCurve{i};
        xH = allBinCenters{i};
        yH = allHighCurve{i};

        okL = isfinite(xL) & isfinite(yL) & xL <= xMaxPlot;
        okH = isfinite(xH) & isfinite(yH) & xH <= xMaxPlot;

        if sum(okL) >= 2
            lowInterp(k,:) = interp1(xL(okL), yL(okL), xCommon, 'linear', NaN);
        end
        if sum(okH) >= 2
            highInterp(k,:) = interp1(xH(okH), yH(okH), xCommon, 'linear', NaN);
        end

        if doRandomControl
            ctrlL = allCtrlLow{i};
            ctrlH = allCtrlHigh{i};

            if ~isempty(ctrlL)
                for rr = 1:size(ctrlL,1)
                    y = ctrlL(rr,:);
                    ok = isfinite(xL) & isfinite(y) & xL <= xMaxPlot;
                    if sum(ok) >= 2
                        ctrlLowInterp(cLow,:) = interp1(xL(ok), y(ok), xCommon, 'linear', NaN);
                        cLow = cLow + 1;
                    end
                end
            end

            if ~isempty(ctrlH)
                for rr = 1:size(ctrlH,1)
                    y = ctrlH(rr,:);
                    ok = isfinite(xH) & isfinite(y) & xH <= xMaxPlot;
                    if sum(ok) >= 2
                        ctrlHighInterp(cHigh,:) = interp1(xH(ok), y(ok), xCommon, 'linear', NaN);
                        cHigh = cHigh + 1;
                    end
                end
            end
        end
    end

    if doRandomControl
        ctrlLowInterp = ctrlLowInterp(1:cLow-1,:);
        ctrlHighInterp = ctrlHighInterp(1:cHigh-1,:);
    end

    % summary stats
    mLow  = mean(lowInterp,  1, 'omitnan');
    sLow  = std(lowInterp,   0, 1, 'omitnan') ./ sqrt(sum(~isnan(lowInterp),1));

    mHigh = mean(highInterp, 1, 'omitnan');
    sHigh = std(highInterp,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(highInterp),1));

    if doRandomControl && ~isempty(ctrlLowInterp)
        mCtrlLow = mean(ctrlLowInterp, 1, 'omitnan');
        sCtrlLow = std(ctrlLowInterp,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(ctrlLowInterp),1));
    else
        mCtrlLow = nan(size(xCommon));
        sCtrlLow = nan(size(xCommon));
    end

    if doRandomControl && ~isempty(ctrlHighInterp)
        mCtrlHigh = mean(ctrlHighInterp, 1, 'omitnan');
        sCtrlHigh = std(ctrlHighInterp,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(ctrlHighInterp),1));
    else
        mCtrlHigh = nan(size(xCommon));
        sCtrlHigh = nan(size(xCommon));
    end

    fig1 = figure('Name','corr_vs_distance shaded summaries', ...
        'Position',[100 100 900 900], 'Color','w');
    tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

    % ---------------- top-left: low vs shuffled low
    nexttile; hold on;
    plot_shaded_sem(xCommon, mLow, sLow, [0.30 0.30 1.00], 0.20);
    if doRandomControl && any(isfinite(mCtrlLow))
        plot_shaded_sem(xCommon, mCtrlLow, sCtrlLow, [0.70 0.50 0.80], 0.25);
    end
    xlim([0 150]);
    xlabel('distance');
    ylabel('Pairwise correlation');
    title('least active');
    box off;

    % ---------------- top-right: high vs shuffled high
    nexttile; hold on;
    plot_shaded_sem(xCommon, mHigh, sHigh, [0.20 0.90 0.90], 0.20);
    if doRandomControl && any(isfinite(mCtrlHigh))
        plot_shaded_sem(xCommon, mCtrlHigh, sCtrlHigh, [0.70 0.50 0.80], 0.25);
    end
    xlim([0 150]);
    xlabel('distance');
    ylabel('Pairwise correlation');
    title('most active');
    box off;

    % ---------------- bottom: low vs high
    nexttile([1 2]); hold on;
    plot_shaded_sem(xCommon, mLow, sLow, [0.30 0.30 1.00], 0.20);
    plot_shaded_sem(xCommon, mHigh, sHigh, [0.20 0.90 0.90], 0.20);
    xlim([0 150]);
    xlabel('distance');
    ylabel('Pairwise correlation');
    title('least active vs most active');
    legend({'least active','most active'}, 'Location','northeast');
    box off;

    sgtitle('corr vs distance: tail-defined inactive vs active windows');

    if saveOutputs
        saveas(fig1, fullfile(outFolder, 'summary_corr_vs_distance_shaded.png'));
    end
end

    % ---------------------------------------------------------------------
    % 6.2 Scalar summaries
    % ---------------------------------------------------------------------
    if any(validRecScalars)

        fig2 = figure('Name','corr_vs_distance scalar summaries', ...
            'Position',[100 100 1400 500], 'Color','w');
        tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

        xPair = [1 2];

        % ---------------- near mean
        nexttile; hold on;
        yLow = R.low_corr_near_mean(validRecScalars);
        yHigh = R.high_corr_near_mean(validRecScalars);
        boxchart(ones(numel(yLow),1)*1, yLow);
        boxchart(ones(numel(yHigh),1)*2, yHigh);
        for i = idxScalar
            plot(xPair, [R.low_corr_near_mean(i), R.high_corr_near_mean(i)], '-o', ...
                'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
        end
        xlim([0.5 2.5]);
        set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
        ylabel(sprintf('Mean corr (first %d bins)', nNearBins));
        title(sprintf('Short-range correlation\nn=%d, p=%.3g %s', ...
            statsSummary.nearMean.n, statsSummary.nearMean.p, p_to_stars(statsSummary.nearMean.p)));
        add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.nearMean.p);

        % ---------------- all mean
        nexttile; hold on;
        yLow = R.low_corr_all_mean(validRecScalars);
        yHigh = R.high_corr_all_mean(validRecScalars);
        boxchart(ones(numel(yLow),1)*1, yLow);
        boxchart(ones(numel(yHigh),1)*2, yHigh);
        for i = idxScalar
            plot(xPair, [R.low_corr_all_mean(i), R.high_corr_all_mean(i)], '-o', ...
                'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
        end
        xlim([0.5 2.5]);
        set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
        ylabel('Mean corr (all bins)');
        title(sprintf('All-bin mean correlation\nn=%d, p=%.3g %s', ...
            statsSummary.allMean.n, statsSummary.allMean.p, p_to_stars(statsSummary.allMean.p)));
        add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.allMean.p);

        % ---------------- auc
        nexttile; hold on;
        if any(validRecAUC)
            yLow = R.low_corr_auc(validRecAUC);
            yHigh = R.high_corr_auc(validRecAUC);
            boxchart(ones(numel(yLow),1)*1, yLow);
            boxchart(ones(numel(yHigh),1)*2, yHigh);
            for i = idxAUC
                plot(xPair, [R.low_corr_auc(i), R.high_corr_auc(i)], '-o', ...
                    'Color',[0.4 0.4 0.4], 'MarkerFaceColor','w');
            end
            xlim([0.5 2.5]);
            set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
            ylabel('Curve AUC');
            title(sprintf('corr-distance AUC\nn=%d, p=%.3g %s', ...
                statsSummary.auc.n, statsSummary.auc.p, p_to_stars(statsSummary.auc.p)));
            add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.auc.p);
        else
            text(0.5, 0.5, 'No valid AUC data', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize', 12);
            axis off;
            title('corr-distance AUC');
        end

        sgtitle('corr vs distance: tail-defined inactive vs active windows');

        if saveOutputs
            saveas(fig2, fullfile(outFolder, 'summary_corr_vs_distance_scalars.png'));
        end
    else
        warning('No valid recordings available for corr-vs-distance scalar summary figure.');
    end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% ================== HELPERS ==================

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

function positions = extract_positions(S)
% Robust extraction of neuron positions

    positions = [];

    vars = fieldnames(S);

    preferredNames = { ...
        'position', 'pos', 'allPos', 'xyz', 'XYZ', ...
        'centroids', 'centroid', 'coords', 'coordinates', ...
        'statPos', 'cellPositions', 'roi_positions'};

    for i = 1:numel(preferredNames)
        nm = preferredNames{i};
        if isfield(S, nm) && isnumeric(S.(nm)) && ismatrix(S.(nm))
            P = double(S.(nm));
            if min(size(P)) >= 2
                positions = P;
                return;
            end
        end
    end

    % Search nested structs
    for i = 1:numel(vars)
        v = S.(vars{i});
        if isstruct(v)
            fns = fieldnames(v);
            for j = 1:numel(fns)
                val = v.(fns{j});
                if isnumeric(val) && ismatrix(val) && min(size(val)) >= 2
                    nm = lower(fns{j});
                    if contains(nm, 'pos') || contains(nm, 'coord') || contains(nm, 'cent')
                        positions = double(val);
                        return;
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

function plot_shaded_sem(x, m, s, lineColor, faceAlpha)
% x, m, s are row vectors

    x = x(:)';
    m = m(:)';
    s = s(:)';

    ok = isfinite(x) & isfinite(m) & isfinite(s);
    if sum(ok) < 2
        return;
    end

    xx = x(ok);
    mm = m(ok);
    ss = s(ok);

    fill([xx fliplr(xx)], [mm-ss fliplr(mm+ss)], lineColor, ...
        'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
    plot(xx, mm, '-', 'Color', lineColor, 'LineWidth', 2.5);
end