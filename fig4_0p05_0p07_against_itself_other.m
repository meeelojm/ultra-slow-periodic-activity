%% ============================================================
%  FIGURE 4C — Band(0.05–0.07 Hz) correlations using
%  calcsort_autocorr_freq_analysis_v3 for domFreqs
%  ============================================================

%% ---------- CONFIG ----------
baseFolder     = '\\forskning.it.ntnu.no\ntnu\mh-kin\yaksi5\anna\Data\Processed 2P data\mGluR multimodal\cell_detect_data\';
folderPattern  = '^[A-Za-z0-9]{4}_Wt$';   % or 'Het'

origStartFrame = 200;    % canonical start in raw movie domain
origFps        = 2.5;    % canonical fps used for origStartFrame

% Frequency-analysis params (used by your calcsort_* function)
windowSize     = 256;
maxPeakCount   = 15;
minPeakHeight  = 2;
frequencyRange = [0.01 0.10];     % search band for domFreqs

% Band of interest for the test
lfs  = 0.05;                      % Hz
lffs = 0.07;                      % Hz

% Output folder
outRoot = fullfile(baseFolder, '_perfish_stability_outputs');
if ~exist(outRoot,'dir'), mkdir(outRoot); end

%% ---------- FIND FISH FOLDERS ----------
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    [~, nm] = fileparts(allFolders{i});
    if ~isempty(regexp(nm, folderPattern, 'once'))
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
    end
end
fprintf('Found %d matching folders.\n', numel(validFolders));

resultsFiles = {};
fishLabel    = {};
for i = 1:numel(validFolders)
    D = dir(fullfile(validFolders{i}, '*_Results_dff.mat'));
    if ~isempty(D)
        resultsFiles{end+1}  = fullfile(validFolders{i}, D(1).name); %#ok<SAGROW>
        [~, fishLabel{end+1}] = fileparts(validFolders{i});         %#ok<SAGROW>
    else
        warning('No *_Results_dff.mat in %s', validFolders{i});
    end
end
nF = numel(resultsFiles);
fprintf('Will analyze %d fish.\n', nF);

%% ---------- ACCUMULATORS ----------
BandIdx_perFish    = cell(nF,1);
domFreqs_perFish   = cell(nF,1);
CorrMats_perFish   = cell(nF,1);

perFish_bb_mean = nan(nF,1);   % Band–Band mean corr (per fish)
perFish_br_mean = nan(nF,1);   % Band–Rest mean corr (per fish)
perFish_bb_med  = nan(nF,1);
perFish_br_med  = nan(nF,1);
perFish_nBand   = nan(nF,1);

AllBandCorrs     = [];         % pooled vectors for hist/stats
AllBandRestCorrs = [];

%% ---------- PER-FISH LOOP ----------
for f = 1:nF
    fprintf('\n[%d/%d] %s\n', f, nF, fishLabel{f});
    S = load(resultsFiles{f});

    % DFF (neurons x frames)
    if isfield(S,'results') && isfield(S.results,'DV_DFFmovwindow')
        dff = double(S.results.DV_DFFmovwindow);
    else
        error('Expected S.results.DV_DFFmovwindow in %s', resultsFiles{f});
    end
    if size(dff,1) > size(dff,2), dff = dff.'; end

    % fps (fixed)
    fps = 2.5;

    % Four fixed 2-min periods starting from a canonical start
    start_time_sec = origStartFrame / origFps;
    start_frame    = max(1, round(start_time_sec * fps));
    framesPerPer   = round(120 * fps);
    P = struct( ...
        't0', start_frame + [0,1,2,3]*framesPerPer, ...
        't1', start_frame + [1,2,3,4]*framesPerPer );
    P.t0 = min(P.t0, size(dff,2));
    P.t1 = min(P.t1, size(dff,2));

    % === Frequency window for calcsort_* ===
    startFrameFreq = P.t0(1);        % use P1–P4 concatenated window
    endFrameFreq   = P.t1(4);
    startFrameFreq = max(1, startFrameFreq);
    endFrameFreq   = min(size(dff,2), endFrameFreq);

    % Number of lags for autocorr-based analysis:
    % choose ~1.5 cycles of the slowest period to give the peak room
    % slowest period ~ 1/fmin
    numLagsFreq = ceil(1.5 * fps / max(frequencyRange(1), eps));

    % --- Recalculate domFreqs with your function ---
    % NOTE: calcsort_* expects dff_red; we map our dff -> dff_red
    dff_red = dff;
    [autocorrs, domFreqs, powerSpecs, freqBins, isoIdx] = ...
        calcsort_autocorr_freq_analysis_v3( ...
            dff_red(:, startFrameFreq:endFrameFreq), ...
            numLagsFreq, fps, ...
            windowSize, maxPeakCount, minPeakHeight, frequencyRange); %#ok<ASGLU>

    domF = double(domFreqs(:));
    if numel(domF) ~= size(dff,1)
        warning('domFreqs length (%d) != neurons (%d) in %s — trimming/padding with NaN.', ...
                numel(domF), size(dff,1), fishLabel{f});
        domF = nan(size(dff,1),1);
        domF(1:min(end, numel(domFreqs))) = domFreqs(1:min(end, numel(domFreqs)));
    end
    domFreqs_perFish{f} = domF;

    % Select band neurons (0.05–0.07 Hz)
    inBand = isfinite(domF) & domF >= lfs & domF <= lffs;
    nBand  = sum(inBand);
    BandIdx_perFish{f} = inBand;
    perFish_nBand(f)   = nBand;
    fprintf('  Found %d band cells out of %d total.\n', nBand, numel(domF));

    % Full Pearson correlation (neurons × neurons), pairwise to tolerate NaNs
    C_inBand = corr(dff(inBand,:).', 'Rows','pairwise');
    C_outBand = corr(dff(inBand,:).', dff(~inBand,:).', 'Rows','pairwise');
    C_inBand(1:size(C_inBand,1)+1:end) = NaN;          % ignore self-corr
    C_outBand(1:size(C_outBand,1)+1:end) = NaN;          % ignore self-corr
    CorrMats_perFish{f,1} = C_inBand;
    CorrMats_perFish{f,2} = C_outBand;

    % sizes
    nBand = sum(inBand);
    nRest = size(dff,1) - nBand;
    
    % optional: remove self-corr from the square in-band matrix
    if nBand > 0
        C_inBand(1:nBand+1:end) = NaN;   % diagonal only exists for the square matrix
    end
    % NOTE: C_outBand is rectangular (nBand × nRest) — it has no self-diagonal to remove.
    
    % ---- Band–Band (upper triangle only) ----
    if nBand >= 2
        ut    = triu(true(nBand), 1);
        bbVec = C_inBand(ut);
    else
        bbVec = [];
    end
    
    % ---- Band–Rest (all entries; no duplicates, no diagonal) ----
    if nBand >= 1 && nRest >= 1
        brVec = C_outBand(:);
    else
        brVec = [];
    end
    
    % ---- save per-fish matrices if you want both blocks later ----
    CorrMats_perFish{f,1} = C_inBand;      % band × band
    CorrMats_perFish{f,2} = C_outBand;     % band × rest
    
    % ---- per-fish summaries ----
    perFish_bb_mean(f) = mean(bbVec,'omitnan');
    perFish_br_mean(f) = mean(brVec,'omitnan');
    perFish_bb_med(f)  = median(bbVec,'omitnan');
    perFish_br_med(f)  = median(brVec,'omitnan');
    
    % ---- pooled vectors for group stats ----
    AllBandCorrs     = [AllBandCorrs;     bbVec(:)];
    AllBandRestCorrs = [AllBandRestCorrs; brVec(:)];

end

%% ---------- SAVE ----------
save(fullfile(outRoot,'fig4C_corrBlocks_perFish_CALCSORT.mat'), ...
    'CorrMats_perFish','BandIdx_perFish','domFreqs_perFish','fishLabel', ...
    'perFish_bb_mean','perFish_br_mean','perFish_bb_med','perFish_br_med','perFish_nBand', ...
    'AllBandCorrs','AllBandRestCorrs', ...
    'lfs','lffs','frequencyRange','windowSize','maxPeakCount','minPeakHeight');

%% ---------- GROUP STATS ----------
keep = isfinite(perFish_bb_mean) & isfinite(perFish_br_mean) & perFish_nBand>0;
fprintf('\n=== Band (0.05–0.07 Hz) correlation advantage (calcsort_*) ===\n');
fprintf('Fish with band cells: %d / %d\n', sum(keep), nF);
delta = perFish_bb_mean(keep) - perFish_br_mean(keep);
fprintf('Δ(mean) BB–BR: median = %.3f  [IQR %.3f–%.3f]\n', ...
    median(delta,'omitnan'), prctile(delta,25), prctile(delta,75));
[p,~,st] = signrank(perFish_bb_mean(keep), perFish_br_mean(keep));

if isfield(st,'zval')
    fprintf('Wilcoxon signed-rank (per-fish means): p=%.3g (z=%.2f)\n', p, st.zval);
elseif isfield(st,'signedrank')
    fprintf('Wilcoxon signed-rank (per-fish means): p=%.3g (signedrank=%d)\n', p, st.signedrank);
else
    fprintf('Wilcoxon signed-rank (per-fish means): p=%.3g\n', p);
end

A = [perFish_bb_mean(keep), perFish_br_mean(keep)];
conds = {'Band–Band','Band–Rest'};
N = size(A,1);

figure('Color','w'); hold on; grid on; box on;

% prepare grouping vector for boxchart
groups = repelem(1:2, N)';   % column of 1 1 … 2 2 …
values = A(:);               % all data stacked

% boxplots with whiskers
boxchart(groups, values, 'BoxFaceAlpha',0.2);

% paired lines
for i=1:N
    plot(1:2, A(i,:), '-', 'Color',[.7 .7 .7]);
end

% scatter points (on top)
scatter(groups, values, 40, 'filled','k');

xlim([0.5 2.5]);
set(gca,'XTick',1:2,'XTickLabel',conds);
ylabel('Mean correlation');

