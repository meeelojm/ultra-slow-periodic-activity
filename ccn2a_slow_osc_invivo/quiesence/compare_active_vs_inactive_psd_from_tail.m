%% =========================================================================
% compare_active_vs_inactive_psd_from_tail.m
%
% SELF-CONTAINED SCRIPT
%
% What it does:
%   1) Recursively finds recordings that contain:
%        - tail_quick.mat
%        - dffs_repact_respcells.mat
%   2) Loads tail data and neural dff data
%   3) Extracts tail bout onsets
%   4) Uses sliding windows to rank periods from least active to most active
%      based on tail bout onsets
%   5) Selects:
%        - least active window
%        - most active window
%   6) Runs calcsort_autocorr_freq_analysis_v3 on both windows
%   7) Extracts per-recording summaries:
%        - bout count
%        - bout rate
%        - median dominant frequency
%        - median band power in frequencyRange
%   8) Makes summary plots comparing most active vs least active periods
%
% IMPORTANT:
%   This script assumes you already have:
%       calcsort_autocorr_freq_analysis_v3.m
%   on your MATLAB path.
%
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
winMin     = 20;   % window length in minutes
stepMin    = 5;    % overlap step in minutes

% PSD/autocorr settings for calcsort_autocorr_freq_analysis_v3
%numLags        = floor(fps*(winMin/2));            % adjust if needed
windowSize     = 1024;            % adjust if needed
maxPeakCount   = 15;             % adjust if needed
minPeakHeight  = 2;              % adjust if needed
frequencyRange = [0.01 0.05];    % ultraslow band of interest (Hz)

% Whether to exclude neurons with too many NaNs inside selected windows
maxNanFracPerNeuron = 0.20;  % remove neurons with >20% NaNs in selected window

% Whether to save outputs
saveOutputs = false;
outFolder   = fullfile(rootFolder, 'active_vs_inactive_tail_psd_results');

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
        recordings(end+1).folder   = thisFolder; %#ok<SAGROW>
        recordings(end).tailPath   = fullfile(tailFiles(i).folder, tailFiles(i).name);
        recordings(end).dffPath    = dffMap(thisFolder);
    end
end

if isempty(recordings)
    error('No folders were found that contain BOTH %s and %s.', tailFileName, dffFileName);
end

fprintf('Found %d matched recordings.\n', numel(recordings));

%% =========================
% 2) PREALLOCATE RESULTS
% =========================

R = struct();
R.recordingLabel          = cell(numel(recordings),1);
R.folder                  = cell(numel(recordings),1);

R.nNeurons_total          = nan(numel(recordings),1);
R.nNeurons_used_low       = nan(numel(recordings),1);
R.nNeurons_used_high      = nan(numel(recordings),1);

R.fps                     = nan(numel(recordings),1);
R.nFrames                 = nan(numel(recordings),1);
R.durationMin             = nan(numel(recordings),1);

R.low_startFrame          = nan(numel(recordings),1);
R.low_endFrame            = nan(numel(recordings),1);
R.high_startFrame         = nan(numel(recordings),1);
R.high_endFrame           = nan(numel(recordings),1);

R.low_boutCount           = nan(numel(recordings),1);
R.high_boutCount          = nan(numel(recordings),1);
R.low_boutRate_perMin     = nan(numel(recordings),1);
R.high_boutRate_perMin    = nan(numel(recordings),1);

R.low_domFreq_median      = nan(numel(recordings),1);
R.high_domFreq_median     = nan(numel(recordings),1);

R.low_bandPower_median    = nan(numel(recordings),1);
R.high_bandPower_median   = nan(numel(recordings),1);

R.low_bandPower_mean      = nan(numel(recordings),1);
R.high_bandPower_mean     = nan(numel(recordings),1);

R.low_totalPSD_mean       = [];
R.high_totalPSD_mean      = [];
R.freqBins                = [];

allLowDomFreqs   = cell(numel(recordings),1);
allHighDomFreqs  = cell(numel(recordings),1);
allLowBandPower  = cell(numel(recordings),1);
allHighBandPower = cell(numel(recordings),1);

%% =========================
% 3) PROCESS EACH RECORDING
% =========================

for r = 1:numel(recordings)

    fprintf('\n============================================================\n');
    fprintf('Recording %d / %d\n', r, numel(recordings));
    fprintf('%s\n', recordings(r).folder);

    R.folder{r} = recordings(r).folder;
    R.recordingLabel{r} = make_recording_label(recordings(r).folder);

    %% 3.1 Load DFF and FPS
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

    % Make sure neurons x frames
    if size(dff,1) > size(dff,2)
        % usually neurons << frames, so if opposite, transpose
        fprintf('Transposing dff because rows > columns.\n');
        dff = dff';
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

    %% 3.2 Load tail data and extract bout onsets
    S_tail = load(recordings(r).tailPath);
    [boutOnsetsFrames, tailActivity] = extract_bout_onsets_from_tail(S_tail, nFrames, fps);

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

    %% 3.3 Rank sliding windows by activity
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

    %% 3.4 Select DFF segments and remove neurons with too many NaNs
    dff_low  = dff(:, low_startFrame:low_endFrame);
    dff_high = dff(:, high_startFrame:high_endFrame);

    keepLow  = mean(isnan(dff_low),  2) <= maxNanFracPerNeuron;
    keepHigh = mean(isnan(dff_high), 2) <= maxNanFracPerNeuron;

    dff_low_use  = dff_low(keepLow, :);
    dff_high_use = dff_high(keepHigh, :);

    % Also remove rows that are all NaN
    dff_low_use  = dff_low_use(~all(isnan(dff_low_use),2), :);
    dff_high_use = dff_high_use(~all(isnan(dff_high_use),2), :);

    R.nNeurons_used_low(r)  = size(dff_low_use,1);
    R.nNeurons_used_high(r) = size(dff_high_use,1);
    numLags = round(0.5 * winMin * 60 * fps);

    if isempty(dff_low_use) || isempty(dff_high_use)
        warning('No usable neurons after NaN filtering in %s. Skipping.', recordings(r).folder);
        continue;
    end

    %% 3.5 Run autocorr / PSD analysis on least active window
    try
        [autocorrs_low, domFreqs_low, powerSpecs_low, freqBins_low, isoIdx_low] = ...
            calcsort_autocorr_freq_analysis_v3( ...
                dff_low_use, ...
                numLags, fps, ...
                windowSize, maxPeakCount, minPeakHeight, frequencyRange); close
    catch ME
        warning('Failed low-activity PSD analysis in %s\n%s', recordings(r).folder, ME.message);
        continue;
    end

    %% 3.6 Run autocorr / PSD analysis on most active window
    try
        [autocorrs_high, domFreqs_high, powerSpecs_high, freqBins_high, isoIdx_high] = ...
            calcsort_autocorr_freq_analysis_v3( ...
                dff_high_use, ...
                numLags, fps, ...
                windowSize, maxPeakCount, minPeakHeight, frequencyRange);  close
    catch ME
        warning('Failed high-activity PSD analysis in %s\n%s', recordings(r).folder, ME.message);
        continue;
    end

    % Basic consistency check
    if isempty(freqBins_low) || isempty(freqBins_high)
        warning('Empty freq bins for %s. Skipping summary extraction.', recordings(r).folder);
        continue;
    end

    if isempty(R.freqBins)
        R.freqBins = freqBins_low(:);
    end

%% 3.7 Extract band power and dominant frequency summaries

% ----- frequency vectors -----
if ismatrix(freqBins_low) && size(freqBins_low,2) > 1
    freqVec_low = freqBins_low(:,1);
else
    freqVec_low = freqBins_low(:);
end

if ismatrix(freqBins_high) && size(freqBins_high,2) > 1
    freqVec_high = freqBins_high(:,1);
else
    freqVec_high = freqBins_high(:);
end

% ----- ensure PSD arrays are [nFreq x nNeurons] -----
if size(powerSpecs_low,1) ~= numel(freqVec_low) && size(powerSpecs_low,2) == numel(freqVec_low)
    powerSpecs_low = powerSpecs_low.';
end

if size(powerSpecs_high,1) ~= numel(freqVec_high) && size(powerSpecs_high,2) == numel(freqVec_high)
    powerSpecs_high = powerSpecs_high.';
end

% ----- safety checks -----
if size(powerSpecs_low,1) ~= numel(freqVec_low)
    warning('LOW mismatch in %s: PSD rows = %d, freq bins = %d. Skipping summary extraction.', ...
        recordings(r).folder, size(powerSpecs_low,1), numel(freqVec_low));
    continue;
end

if size(powerSpecs_high,1) ~= numel(freqVec_high)
    warning('HIGH mismatch in %s: PSD rows = %d, freq bins = %d. Skipping summary extraction.', ...
        recordings(r).folder, size(powerSpecs_high,1), numel(freqVec_high));
    continue;
end

% ----- frequency masks -----
bandMask_low  = freqVec_low  >= frequencyRange(1) & freqVec_low  <= frequencyRange(2);
bandMask_high = freqVec_high >= frequencyRange(1) & freqVec_high <= frequencyRange(2);

fprintf('  low band bins:  %d / %d\n', sum(bandMask_low),  numel(bandMask_low));
fprintf('  high band bins: %d / %d\n', sum(bandMask_high), numel(bandMask_high));

if ~any(bandMask_low) || ~any(bandMask_high)
    warning('No PSD bins fell inside frequencyRange for %s. Skipping summary extraction.', ...
        recordings(r).folder);
    continue;
end

% ----- band power per neuron -----
bandPower_low  = nansum(powerSpecs_low(bandMask_low, :), 1).';
bandPower_high = nansum(powerSpecs_high(bandMask_high, :), 1).';

% ----- dominant frequency within band per neuron -----
Pband_low  = powerSpecs_low(bandMask_low, :);
Pband_high = powerSpecs_high(bandMask_high, :);

fband_low  = freqVec_low(bandMask_low);
fband_high = freqVec_high(bandMask_high);

[~, idx_low]  = max(Pband_low,  [], 1);
[~, idx_high] = max(Pband_high, [], 1);

domFreqs_low_band  = fband_low(idx_low(:));
domFreqs_high_band = fband_high(idx_high(:));

% ----- store per-neuron outputs -----
allLowDomFreqs{r}   = domFreqs_low_band;
allHighDomFreqs{r}  = domFreqs_high_band;
allLowBandPower{r}  = bandPower_low;
allHighBandPower{r} = bandPower_high;

% ----- recording-level summaries -----
R.low_domFreq_median(r)   = nanmedian(domFreqs_low_band);
R.high_domFreq_median(r)  = nanmedian(domFreqs_high_band);

R.low_bandPower_median(r)  = nanmedian(bandPower_low);
R.high_bandPower_median(r) = nanmedian(bandPower_high);

R.low_bandPower_mean(r)  = nanmean(bandPower_low);
R.high_bandPower_mean(r) = nanmean(bandPower_high);

% ----- mean PSD across neurons for this recording -----
R.low_totalPSD_mean(:,r)  = nanmean(powerSpecs_low,  2);
R.high_totalPSD_mean(:,r) = nanmean(powerSpecs_high, 2);

if isempty(R.freqBins)
    R.freqBins = freqVec_low;
end

fprintf('  low domFreq median  = %.4f Hz\n', R.low_domFreq_median(r));
fprintf('  high domFreq median = %.4f Hz\n', R.high_domFreq_median(r));
fprintf('  low bandPower mean  = %.6f\n', R.low_bandPower_mean(r));
fprintf('  high bandPower mean = %.6f\n', R.high_bandPower_mean(r));

  %% 3.8 Optional per-recording diagnostic figure
    fig1 = figure('Name', ['Activity windows - ' R.recordingLabel{r}], ...
        'Position', [100 100 1400 800]);
    
    tiledlayout(3,2, 'Padding','compact', 'TileSpacing','compact');
    
    % -------------------------------------------------------------------------
    % Prep vectors for plotting
    % -------------------------------------------------------------------------
    if ismatrix(freqBins_low) && size(freqBins_low,2) > 1
        freqVec_low = freqBins_low(:,1);
    else
        freqVec_low = freqBins_low(:);
    end
    
    if ismatrix(freqBins_high) && size(freqBins_high,2) > 1
        freqVec_high = freqBins_high(:,1);
    else
        freqVec_high = freqBins_high(:);
    end
    
    % Ensure PSD is [nFreq x nNeurons]
    if size(powerSpecs_low,1) ~= numel(freqVec_low) && size(powerSpecs_low,2) == numel(freqVec_low)
        powerSpecs_low = powerSpecs_low.';
    end
    if size(powerSpecs_high,1) ~= numel(freqVec_high) && size(powerSpecs_high,2) == numel(freqVec_high)
        powerSpecs_high = powerSpecs_high.';
    end
    
    % Mean PSD across neurons
    meanPSD_low  = nanmean(powerSpecs_low,  2);
    meanPSD_high = nanmean(powerSpecs_high, 2);
    
    % Limit PSD display to first 30 bins or less
    nShow = min([30, numel(freqVec_low), numel(freqVec_high), numel(meanPSD_low), numel(meanPSD_high)]);
    fShow_low   = freqVec_low(1:nShow);
    fShow_high  = freqVec_high(1:nShow);
    psdShow_low  = meanPSD_low(1:nShow);
    psdShow_high = meanPSD_high(1:nShow);
    
    % Lag axis for mean autocorrelation
    lagSec_low  = (0:size(autocorrs_low,2)-1)  ./ fps;
    lagSec_high = (0:size(autocorrs_high,2)-1) ./ fps;
    
    % -------------------------------------------------------------------------
    % Tail activity trace
    % -------------------------------------------------------------------------
    nexttile;
    if ~isempty(tailActivity)
        tmin = (0:numel(tailActivity)-1) / fps / 60;
        plot(tmin, tailActivity, 'k');
        hold on;
        xline(low_startFrame/fps/60,  '--b', 'Low start');
        xline(low_endFrame/fps/60,    '--b', 'Low end');
        xline(high_startFrame/fps/60, '--r', 'High start');
        xline(high_endFrame/fps/60,   '--r', 'High end');
        xlabel('Time (min)');
        ylabel('Tail activity');
        title('Tail activity with selected windows');
    else
        text(0.1,0.5,'No continuous tail activity trace found','FontSize',12);
        axis off;
    end
    
    % -------------------------------------------------------------------------
    % Sliding bout count
    % -------------------------------------------------------------------------
    nexttile;
    plot(winTable.centerFrame/fps/60, winTable.boutCount, '-k', 'LineWidth', 1.2);
    hold on;
    scatter(winTable.centerFrame(lowIdx)/fps/60,  winTable.boutCount(lowIdx),  80, 'b', 'filled');
    scatter(winTable.centerFrame(highIdx)/fps/60, winTable.boutCount(highIdx), 80, 'r', 'filled');
    xlabel('Window center (min)');
    ylabel('Bout count');
    title('Sliding-window bout count');
    
    % -------------------------------------------------------------------------
    % Mean PSD
    % -------------------------------------------------------------------------
    nexttile;
    plot(fShow_low,  psdShow_low,  'b', 'LineWidth', 1.5); hold on;
    plot(fShow_high, psdShow_high, 'r', 'LineWidth', 1.5);
    
    xt = 1:3:nShow;
    xticks(fShow_low(xt));
    xticklabels(compose('%.3f', fShow_low(xt)));
    
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title('Mean PSD');
    legend({'Least active','Most active'}, 'Location','best');
    xline(frequencyRange(1), '--k');
    xline(frequencyRange(2), '--k');
    
    % -------------------------------------------------------------------------
    % Dominant frequency distribution in 4 categorical bins
    % -------------------------------------------------------------------------
    nexttile;
    hold on;
    
    edgesHist = [0, 0.018, 0.055, 0.091, 0.10];
    binLabels = {'0–0.018','0.018–0.055','0.055–0.091','0.091–0.10'};
    nBinsHist = numel(edgesHist) - 1;
    xCats     = 1:nBinsHist;
    
    % Low window
    domLow = domFreqs_low_band(:);
    domLow = domLow(isfinite(domLow) & domLow >= edgesHist(1) & domLow <= edgesHist(end));
    propLow = nan(1,nBinsHist);
    if ~isempty(domLow)
        bLow = discretize(domLow, edgesHist);
        cLow = accumarray(bLow(~isnan(bLow)), 1, [nBinsHist 1], @sum, 0);
        propLow = (cLow ./ sum(cLow)).';
    end
    
    % High window
    domHigh = domFreqs_high_band(:);
    domHigh = domHigh(isfinite(domHigh) & domHigh >= edgesHist(1) & domHigh <= edgesHist(end));
    propHigh = nan(1,nBinsHist);
    if ~isempty(domHigh)
        bHigh = discretize(domHigh, edgesHist);
        cHigh = accumarray(bHigh(~isnan(bHigh)), 1, [nBinsHist 1], @sum, 0);
        propHigh = (cHigh ./ sum(cHigh)).';
    end
    
    % small jitter so blue/red don't fully overlap
    % jitterWidth = 0.12;
    % xLow  = xCats - 0.10 + (rand(size(xCats)) - 0.5) * jitterWidth;
    % xHigh = xCats + 0.10 + (rand(size(xCats)) - 0.5) * jitterWidth;
    % 
    % scatter(xLow(isfinite(propLow)),   propLow(isfinite(propLow)),   60, 'b', 'filled', ...
    %     'MarkerFaceAlpha',0.65, 'MarkerEdgeAlpha',0.65);
    % scatter(xHigh(isfinite(propHigh)), propHigh(isfinite(propHigh)), 60, 'r', 'filled', ...
    %     'MarkerFaceAlpha',0.65, 'MarkerEdgeAlpha',0.65);
    
    plot(xCats, propLow,  '-o', 'Color', 'b', 'LineWidth', 1.8, 'MarkerFaceColor', 'b');
    plot(xCats, propHigh, '-o', 'Color', 'r', 'LineWidth', 1.8, 'MarkerFaceColor', 'r');
    
    xticks(xCats);
    xticklabels(binLabels);
    xlim([0.5, nBinsHist+0.5]);
    xlabel('Dominant frequency bin (Hz)');
    ylabel('Proportion of neurons');
    title('Dominant frequency distribution');
    legend({'Least active','Most active'}, 'Location','best');
    
    % -------------------------------------------------------------------------
    % Band power distributions
    % -------------------------------------------------------------------------
    nexttile;
    hold on;
    
    allBP = [bandPower_low(:); bandPower_high(:)];
    allBP = allBP(isfinite(allBP));
    
    if isempty(allBP)
        text(0.1,0.5,'No finite band-power values','FontSize',12);
        axis off;
    else
        nBins = 20;  % adjust if you want fewer/more bins
        edgesBP = linspace(min(allBP), max(allBP), nBins+1);
    
        histogram(bandPower_low,  edgesBP, ...
            'Normalization','probability', 'FaceAlpha',0.5);
    
        histogram(bandPower_high, edgesBP, ...
            'Normalization','probability', 'FaceAlpha',0.5);
    
        xlabel(sprintf('Band power [%.3f %.3f] Hz', frequencyRange(1), frequencyRange(2)));
        ylabel('Probability');
        title('Band power distribution');
        legend({'Least active','Most active'}, 'Location','best');
    end
    
    % -------------------------------------------------------------------------
    % Mean autocorrelation
    % -------------------------------------------------------------------------
    nexttile;
    if ~isempty(autocorrs_low) && ~isempty(autocorrs_high)
        plot(lagSec_low,  nanmean(autocorrs_low,1),  'b', 'LineWidth',1.5); hold on;
        plot(lagSec_high, nanmean(autocorrs_high,1), 'r', 'LineWidth',1.5);
        xlabel('Lag (s)'); xlim([0 100])
        ylabel('Autocorrelation');
        title('Mean autocorrelation');
        legend({'Least active','Most active'}, 'Location','best');
    else
        text(0.1,0.5,'Autocorr output empty','FontSize',12);
        axis off;
    end
    
    sgtitle(R.recordingLabel{r}, 'Interpreter','none');
    
    if saveOutputs
        saveas(fig1, fullfile(outFolder, ...
            sprintf('diagnostic_%02d_%s.png', r, sanitize_filename(R.recordingLabel{r}))));
    end

end

%% =========================
% 4) BUILD SUMMARY TABLE
% =========================

T = table( ...
    R.recordingLabel, ...
    R.folder, ...
    R.nNeurons_total, ...
    R.nNeurons_used_low, ...
    R.nNeurons_used_high, ...
    R.fps, ...
    R.nFrames, ...
    R.durationMin, ...
    R.low_startFrame, ...
    R.low_endFrame, ...
    R.high_startFrame, ...
    R.high_endFrame, ...
    R.low_boutCount, ...
    R.high_boutCount, ...
    R.low_boutRate_perMin, ...
    R.high_boutRate_perMin, ...
    R.low_domFreq_median, ...
    R.high_domFreq_median, ...
    R.low_bandPower_median, ...
    R.high_bandPower_median, ...
    R.low_bandPower_mean, ...
    R.high_bandPower_mean, ...
    'VariableNames', { ...
    'recordingLabel', 'folder', ...
    'nNeurons_total', 'nNeurons_used_low', 'nNeurons_used_high', ...
    'fps', 'nFrames', 'durationMin', ...
    'low_startFrame', 'low_endFrame', ...
    'high_startFrame', 'high_endFrame', ...
    'low_boutCount', 'high_boutCount', ...
    'low_boutRate_perMin', 'high_boutRate_perMin', ...
    'low_domFreq_median', 'high_domFreq_median', ...
    'low_bandPower_median', 'high_bandPower_median', ...
    'low_bandPower_mean', 'high_bandPower_mean'});

disp(T);

if saveOutputs
    writetable(T, fullfile(outFolder, 'summary_active_vs_inactive.csv'));
    save(fullfile(outFolder, 'summary_active_vs_inactive.mat'), ...
        'R', 'T', 'allLowDomFreqs', 'allHighDomFreqs', 'allLowBandPower', 'allHighBandPower');
end

%% =========================
% 5) SUMMARY FIGURES
% =========================

validRec = ~isnan(R.low_boutCount) & ~isnan(R.high_boutCount) & ...
           ~isnan(R.low_domFreq_median) & ~isnan(R.high_domFreq_median) & ...
           ~isnan(R.low_bandPower_median) & ~isnan(R.high_bandPower_median);

%% =========================
% 6) PAIRED STATISTICS + SUMMARY FIGURES
% =========================

% -------------------------------------------------------------------------
% helper for paired signrank stats
% -------------------------------------------------------------------------
statsSummary = struct();

statsSummary.boutCount       = run_paired_signrank(R.low_boutCount,        R.high_boutCount);
statsSummary.boutRate        = run_paired_signrank(R.low_boutRate_perMin,  R.high_boutRate_perMin);
statsSummary.domFreq         = run_paired_signrank(R.low_domFreq_median,   R.high_domFreq_median);
statsSummary.bandPowerMedian = run_paired_signrank(R.low_bandPower_median, R.high_bandPower_median);
statsSummary.bandPowerMean   = run_paired_signrank(R.low_bandPower_mean,   R.high_bandPower_mean);

disp('================ PAIRED SIGNRANK TESTS ================')
fprintf('Bout count:            n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.boutCount.n, statsSummary.boutCount.p, statsSummary.boutCount.deltaMedian);
fprintf('Bout rate:             n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.boutRate.n, statsSummary.boutRate.p, statsSummary.boutRate.deltaMedian);
fprintf('Median dom freq:       n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.domFreq.n, statsSummary.domFreq.p, statsSummary.domFreq.deltaMedian);
fprintf('Median band power:     n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.bandPowerMedian.n, statsSummary.bandPowerMedian.p, statsSummary.bandPowerMedian.deltaMedian);
fprintf('Mean band power:       n = %d, p = %.4g, median Δ = %.4g\n', ...
    statsSummary.bandPowerMean.n, statsSummary.bandPowerMean.p, statsSummary.bandPowerMean.deltaMedian);

if saveOutputs
    save(fullfile(outFolder, 'stats_summary.mat'), 'statsSummary');
end

% -------------------------------------------------------------------------
% Figure A: movement / bouts
% -------------------------------------------------------------------------
validRecA = isfinite(R.low_boutCount) & isfinite(R.high_boutCount) & ...
            isfinite(R.low_boutRate_perMin) & isfinite(R.high_boutRate_perMin);

if any(validRecA)

    xPair = [1 2];

    figA = figure('Name','Active vs inactive bout comparison', ...
        'Position',[100 100 900 700]);
    tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

    % ----- bout count
    nexttile; hold on;
    yLow  = R.low_boutCount(validRecA);
    yHigh = R.high_boutCount(validRecA);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    idxA = find(validRecA(:))';
    for i = idxA
        plot(xPair, [R.low_boutCount(i), R.high_boutCount(i)], '-o', ...
            'Color', [0.4 0.4 0.4], 'MarkerFaceColor', 'w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Bout count');
    title(sprintf('Bout count\nn=%d, p=%.3g %s', ...
        statsSummary.boutCount.n, statsSummary.boutCount.p, p_to_stars(statsSummary.boutCount.p)));

    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.boutCount.p);

    % ----- bout rate
    nexttile; hold on;
    yLow  = R.low_boutRate_perMin(validRecA);
    yHigh = R.high_boutRate_perMin(validRecA);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxA
        plot(xPair, [R.low_boutRate_perMin(i), R.high_boutRate_perMin(i)], '-o', ...
            'Color', [0.4 0.4 0.4], 'MarkerFaceColor', 'w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Bouts / min');
    title(sprintf('Bout rate\nn=%d, p=%.3g %s', ...
        statsSummary.boutRate.n, statsSummary.boutRate.p, p_to_stars(statsSummary.boutRate.p)));

    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.boutRate.p);

    sgtitle('Tail-defined inactive vs active periods');

    if saveOutputs
        saveas(figA, fullfile(outFolder, 'summary_bouts_active_vs_inactive.png'));
    end
else
    warning('No valid recordings available for movement summary figure.');
end

% -------------------------------------------------------------------------
% Figure B: spectral summaries
% -------------------------------------------------------------------------
validRecB = isfinite(R.low_domFreq_median) & isfinite(R.high_domFreq_median) & ...
            isfinite(R.low_bandPower_median) & isfinite(R.high_bandPower_median);

if any(validRecB)

    xPair = [1 2];

    figB = figure('Name','Active vs inactive PSD comparison', ...
        'Position',[100 100 1200 700]);
    tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

    % ----- dominant frequency
    nexttile; hold on;
   yLow  = R.low_domFreq_median(validRecB);
yHigh = R.high_domFreq_median(validRecB);
idxThis = find(validRecB(:))';

boxchart(ones(numel(yLow),1)*1, yLow);
boxchart(ones(numel(yHigh),1)*2, yHigh);

jit = 0.04;   % small horizontal jitter

for k = 1:numel(idxThis)
    i = idxThis(k);

    x1 = 1 + (rand-0.5)*2*jit;
    x2 = 2 + (rand-0.5)*2*jit;

    plot([x1 x2], [R.low_domFreq_median(i), R.high_domFreq_median(i)], ...
        '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

    scatter(x1, R.low_domFreq_median(i), 36, 'w', 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3]);
    scatter(x2, R.high_domFreq_median(i), 36, 'w', 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3]);
end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel('Median dominant frequency (Hz)');
    title(sprintf('Dominant frequency\nn=%d, p=%.3g %s', ...
        statsSummary.domFreq.n, statsSummary.domFreq.p, p_to_stars(statsSummary.domFreq.p)));

    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.domFreq.p);

    % ----- band power
    nexttile; hold on;
    yLow  = R.low_bandPower_median(validRecB);
    yHigh = R.high_bandPower_median(validRecB);

    boxchart(ones(numel(yLow),1)*1, yLow);
    boxchart(ones(numel(yHigh),1)*2, yHigh);

    for i = idxB
        plot(xPair, [R.low_bandPower_median(i), R.high_bandPower_median(i)], '-o', ...
            'Color', [0.4 0.4 0.4], 'MarkerFaceColor', 'w');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',xPair,'XTickLabel',{'Least active','Most active'});
    ylabel(sprintf('Median band power [%.2f %.2f] Hz', frequencyRange(1), frequencyRange(2)));
    title(sprintf('Band power\nn=%d, p=%.3g %s', ...
        statsSummary.bandPowerMedian.n, statsSummary.bandPowerMedian.p, p_to_stars(statsSummary.bandPowerMedian.p)));

    add_sig_bar(gca, xPair, [yLow; yHigh], statsSummary.bandPowerMedian.p);

    % ----- average PSD with shaded SEM
    nexttile; hold on;

    freqVec = R.freqBins(:);
    lowPSD  = R.low_totalPSD_mean(:,validRecB);
    highPSD = R.high_totalPSD_mean(:,validRecB);

    nPlot = min([numel(freqVec), size(lowPSD,1), size(highPSD,1)]);
    freqVec = freqVec(1:nPlot);
    lowPSD  = lowPSD(1:nPlot,:);
    highPSD = highPSD(1:nPlot,:);

    lowMean  = nanmean(lowPSD, 2);
    highMean = nanmean(highPSD, 2);

    lowN  = sum(~isnan(lowPSD), 2);
    highN = sum(~isnan(highPSD), 2);

% mean and SD across recordings
lowMean  = nanmean(lowPSD, 2);
highMean = nanmean(highPSD, 2);

lowSD  = nanstd(lowPSD, 0, 2);
highSD = nanstd(highPSD, 0, 2);

lowSEM  = nanstd(lowPSD, 0, 2) ./ sqrt(max(lowN,1));
highSEM = nanstd(highPSD, 0, 2) ./ sqrt(max(highN,1));
% shaded error regions = SD
fill([freqVec; flipud(freqVec)], ...
     [lowMean-lowSEM; flipud(lowMean+lowSEM)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([freqVec; flipud(freqVec)], ...
     [highMean-highSEM; flipud(highMean+highSEM)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% mean lines
plot(freqVec, lowMean,  'b', 'LineWidth', 1.8);
plot(freqVec, highMean, 'r', 'LineWidth', 1.8);

xlabel('Frequency (Hz)');
ylabel('Mean PSD across recordings');
title('Average PSD');
legend({'Least active SD','Most active SD','Least active mean','Most active mean'}, ...
    'Location','best');

xline(frequencyRange(1), '--k');
xline(frequencyRange(2), '--k');
xlim([0 0.2]);

    sgtitle('Neural spectral comparison: tail-defined inactive vs active windows');

    if saveOutputs
        saveas(figB, fullfile(outFolder, 'summary_psd_active_vs_inactive.png'));
    end
else
    warning('No valid recordings available for spectral summary figure.');
end

fprintf('\nDone.\n');
if saveOutputs
    fprintf('Outputs saved to:\n%s\n', outFolder);
end

%% =========================
% LOCAL HELPER FUNCTIONS
% =========================

function out = run_paired_signrank(xLow, xHigh)
    idx = isfinite(xLow) & isfinite(xHigh);
    xLow = xLow(idx);
    xHigh = xHigh(idx);

    out = struct('p',NaN,'h',NaN,'signedrank',NaN,'n',numel(xLow), ...
                 'deltaMedian',NaN,'xLow',xLow,'xHigh',xHigh);

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

    axes(ax); hold(ax, 'on');
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
%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

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

function label = make_recording_label(folderPath)
    parts = regexp(folderPath, '[\\/]', 'split');
    if numel(parts) >= 3
        label = strjoin(parts(end-2:end), '_');
    else
        label = folderPath;
    end
end

function s = sanitize_filename(s)
    s = regexprep(s, '[^\w\-]', '_');
end