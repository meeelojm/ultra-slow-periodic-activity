function OUT = bout_onset_manifold_windows_analysis(rootFolder, varargin)
% BOUT_ONSET_MANIFOLD_WINDOWS_ANALYSIS
% -------------------------------------------------------------------------
% Population-level manifold analysis comparing least-active vs most-active
% tail-bout windows, while fitting the latent manifold on the FULL recording.
%
% For each recording the function:
%   1) finds matching dff + tail files
%   2) extracts bout onsets from tail
%   3) ranks sliding windows by bout rate
%   4) fits PCA / latent manifold on the full recording
%   5) evaluates low-activity and high-activity windows using the SAME latent axes
%   6) computes per-window summaries:
%        - explained variance retained by PCs
%        - latent PSDs and dominant frequencies
%        - PLV matrix
%        - LDS eigenvalues / oscillator candidate pairs
%        - window trajectory slices on full manifold
%   7) saves per-recording figures and population summary figures
%
% IMPORTANT DESIGN CHOICE
% -----------------------
% The latent manifold is computed on the full recording only once. Low/high
% windows are then projected into that fixed latent space. This means window
% differences reflect changes in trajectory occupancy / local dynamics, not
% changes caused by rotating the PCA basis separately in each window.
%
% OUTPUT
% ------
% OUT.summaryTable        : per-recording summary table
% OUT.population          : pooled population metrics
% OUT.recordings          : per-recording structs
% OUT.cfg                 : configuration used
%
% -------------------------------------------------------------------------

%% ========================= PARSE INPUTS ================================
ip = inputParser;
addRequired(ip, 'rootFolder', @(x) ischar(x) || isstring(x));

addParameter(ip, 'TailFileName', 'tail_quick.mat', @(x) ischar(x) || isstring(x));
addParameter(ip, 'DffFileName',  'dffs_repact_respcells.mat', @(x) ischar(x) || isstring(x));

addParameter(ip, 'WinMin', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'StepMin', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'MaxNanFracPerNeuron', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);

addParameter(ip, 'NumPCs', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2);
addParameter(ip, 'BandHz', [0.01 0.10], @(x) isnumeric(x) && numel(x)==2 && x(1) > 0 && x(2) > x(1));
addParameter(ip, 'WelchWinSec', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'WelchOverlapFrac', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
addParameter(ip, 'PhaseAmpQuantileMin', 0.20, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
addParameter(ip, 'PLVThreshold', 0.70, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(ip, 'LDSReg', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(ip, 'MinEigMag', 0.90, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'MaxEigMag', 1.05, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'UseZScore', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'DoPlots', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'SaveFigures', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'OutDirName', 'bout_manifold_analysis', @(x) ischar(x) || isstring(x));

parse(ip, rootFolder, varargin{:});
P = ip.Results;

rootFolder = char(rootFolder);
tailFileName = char(P.TailFileName);
dffFileName  = char(P.DffFileName);

outDir = fullfile(rootFolder, char(P.OutDirName));
figDir = fullfile(outDir, 'figures_per_recording');
popDir = fullfile(outDir, 'population_figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(popDir, 'dir'), mkdir(popDir); end

%% ========================= 1) FIND MATCHING RECORDINGS ================
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
        recordings(end+1).folder  = thisFolder; %#ok<AGROW>
        recordings(end).tailPath  = fullfile(tailFiles(i).folder, tailFiles(i).name);
        recordings(end).dffPath   = dffMap(thisFolder);
    end
end

if isempty(recordings)
    error('No folders contain both %s and %s.', tailFileName, dffFileName);
end

fprintf('Found %d matched recordings.\n', numel(recordings));

%% ========================= 2) PREALLOCATE =============================
nR = numel(recordings);
Rec = repmat(struct(), nR, 1);
rows = repmat(struct(), nR, 1);

modeColors = turbo(max(P.NumPCs, 2));

%% ========================= 3) PROCESS EACH RECORDING ==================
for r = 1:nR
    fprintf('\n============================================================\n');
    fprintf('Recording %d / %d\n', r, nR);
    fprintf('%s\n', recordings(r).folder);

    try
        recLabel = make_recording_label(recordings(r).folder);

        % ----- load dff / fps -----
        S_dff = load(recordings(r).dffPath);
        [dff, fps] = extract_dff_and_fps(S_dff);
        if isempty(dff), error('Could not extract DFF.'); end
        if isempty(fps) || ~isscalar(fps) || ~isfinite(fps) || fps <= 0
            error('Could not extract valid fps.');
        end

        dff = double(dff);
        if size(dff,1) > size(dff,2)
            dff = dff.';
        end

        nNeurons = size(dff,1);
        nFrames  = size(dff,2);
        durMin   = nFrames / fps / 60;

        if P.Verbose
            fprintf('  DFF: %d neurons x %d frames | fps = %.4f | duration = %.2f min\n', ...
                nNeurons, nFrames, fps, durMin);
        end

        % ----- bout onsets from tail -----
        S_tail = load(recordings(r).tailPath);
        [boutOnsetsFrames, tailTrace] = extract_bout_onsets_from_tail(S_tail, nFrames, fps);
        boutOnsetsFrames = unique(round(boutOnsetsFrames(:)));
        boutOnsetsFrames = boutOnsetsFrames(isfinite(boutOnsetsFrames));
        boutOnsetsFrames = boutOnsetsFrames(boutOnsetsFrames >= 1 & boutOnsetsFrames <= nFrames);
        if isempty(boutOnsetsFrames)
            error('No bout onsets extracted.');
        end

        % ----- rank windows -----
        winFrames  = max(1, round(P.WinMin  * 60 * fps));
        stepFrames = max(1, round(P.StepMin * 60 * fps));
        if nFrames < winFrames
            error('Recording shorter than one analysis window.');
        end

        [winTable, lowIdx, highIdx] = rank_windows_by_bout_activity( ...
            boutOnsetsFrames, nFrames, winFrames, stepFrames, fps);
        if isempty(winTable)
            error('Could not rank windows.');
        end

        lowF  = [winTable.startFrame(lowIdx),  winTable.endFrame(lowIdx)];
        highF = [winTable.startFrame(highIdx), winTable.endFrame(highIdx)];

        if P.Verbose
            fprintf('  Least active window:  %d-%d | bouts = %d | %.3f bouts/min\n', ...
                lowF(1), lowF(2), winTable.boutCount(lowIdx), winTable.boutRatePerMin(lowIdx));
            fprintf('  Most active window:   %d-%d | bouts = %d | %.3f bouts/min\n', ...
                highF(1), highF(2), winTable.boutCount(highIdx), winTable.boutRatePerMin(highIdx));
        end

        % ----- neuron QC on BOTH windows, manifold fitted on full recording -----
        dff_low  = dff(:, lowF(1):lowF(2));
        dff_high = dff(:, highF(1):highF(2));
        keepLow  = mean(isnan(dff_low),  2) <= P.MaxNanFracPerNeuron;
        keepHigh = mean(isnan(dff_high), 2) <= P.MaxNanFracPerNeuron;
        keepUse  = keepLow & keepHigh & ~all(isnan(dff),2);
        if nnz(keepUse) < P.NumPCs
            error('Too few neurons survive shared QC.');
        end

        dff_use = dff(keepUse, :);
        dff_use = fill_nans_rowwise(dff_use);

        if P.UseZScore
            Xfull = zscore_fill_rows(dff_use);
        else
            Xfull = dff_use;
            Xfull(~isfinite(Xfull)) = 0;
        end

        % ----- PCA ON FULL RECORDING -----
        maxPC = min([P.NumPCs, rank(Xfull'), size(Xfull,1), size(Xfull,2)]);
        [coeff, score, latent, ~, explained, mu] = pca(Xfull', 'NumComponents', maxPC);
        latentTS = score(:,1:maxPC);           % time x PC
        coeffUse = coeff(:,1:maxPC);           % neuron x PC
        numPCs   = maxPC;

        % neuron mode assignment from full recording
        [modeLoadingAbs, modeID] = max(abs(coeffUse), [], 2);
        [~, sortOrd] = sortrows([modeID(:), -modeLoadingAbs(:)], [1 2]);
        modeMean = compute_mode_mean(Xfull, modeID, numPCs);
        [~, domMode_t] = max(abs(modeMean), [], 1);

        % ----- evaluate FULL / LOW / HIGH using same latent manifold -----
        EvalFull = eval_window_in_fixed_latent(latentTS, 1:size(latentTS,1), fps, P);
        EvalLow  = eval_window_in_fixed_latent(latentTS, lowF(1):lowF(2), fps, P);
        EvalHigh = eval_window_in_fixed_latent(latentTS, highF(1):highF(2), fps, P);

        % explained variance within window, still using global PCs
        winExpLow  = window_pc_variance_fraction(latentTS(lowF(1):lowF(2),:));
        winExpHigh = window_pc_variance_fraction(latentTS(highF(1):highF(2),:));
        winExpFull = window_pc_variance_fraction(latentTS);

        % ----- store recording -----
        Rec(r).label            = recLabel;
        Rec(r).folder           = recordings(r).folder;
        Rec(r).fps              = fps;
        Rec(r).nNeurons         = nNeurons;
        Rec(r).nFrames          = nFrames;
        Rec(r).durationMin      = durMin;
        Rec(r).nKeptShared      = nnz(keepUse);
        Rec(r).boutOnsetsFrames = boutOnsetsFrames;
        Rec(r).tailTrace        = tailTrace;
        Rec(r).winTable         = winTable;
        Rec(r).lowFrames        = lowF;
        Rec(r).highFrames       = highF;
        Rec(r).fullPCA.coeff    = coeffUse;
        Rec(r).fullPCA.latentTS = latentTS;
        Rec(r).fullPCA.explained = explained(1:numPCs);
        Rec(r).fullPCA.latent    = latent(1:numPCs);
        Rec(r).fullPCA.mu        = mu;
        Rec(r).modeID            = modeID;
        Rec(r).modeMean          = modeMean;
        Rec(r).domMode_t         = domMode_t;
        Rec(r).sortOrd           = sortOrd;
        Rec(r).EvalFull          = EvalFull;
        Rec(r).EvalLow           = EvalLow;
        Rec(r).EvalHigh          = EvalHigh;
        Rec(r).WinExpFull        = winExpFull;
        Rec(r).WinExpLow         = winExpLow;
        Rec(r).WinExpHigh        = winExpHigh;

        % ----- summary row -----
        rows(r).recordingLabel       = string(recLabel);
        rows(r).folder               = string(recordings(r).folder);
        rows(r).fps                  = fps;
        rows(r).nNeurons             = nNeurons;
        rows(r).nKeptShared          = nnz(keepUse);
        rows(r).durationMin          = durMin;
        rows(r).lowBoutRatePerMin    = winTable.boutRatePerMin(lowIdx);
        rows(r).highBoutRatePerMin   = winTable.boutRatePerMin(highIdx);
        rows(r).full_meanPLV         = EvalFull.meanPLV;
        rows(r).low_meanPLV          = EvalLow.meanPLV;
        rows(r).high_meanPLV         = EvalHigh.meanPLV;
        rows(r).full_nPairsLDS       = EvalFull.nOscillatorPairsLDS;
        rows(r).low_nPairsLDS        = EvalLow.nOscillatorPairsLDS;
        rows(r).high_nPairsLDS       = EvalHigh.nOscillatorPairsLDS;
        rows(r).full_pcFreqCV        = EvalFull.pcFreqCV;
        rows(r).low_pcFreqCV         = EvalLow.pcFreqCV;
        rows(r).high_pcFreqCV        = EvalHigh.pcFreqCV;
        rows(r).delta_highMinusLow_PLV      = EvalHigh.meanPLV - EvalLow.meanPLV;
        rows(r).delta_highMinusLow_pcFreqCV = EvalHigh.pcFreqCV - EvalLow.pcFreqCV;
        rows(r).delta_highMinusLow_pairs    = EvalHigh.nOscillatorPairsLDS - EvalLow.nOscillatorPairsLDS;
        rows(r).full_pc1_exp         = winExpFull(1);
        rows(r).low_pc1_exp          = winExpLow(1);
        rows(r).high_pc1_exp         = winExpHigh(1);

        % ----- plots per fish -----
        if P.DoPlots
            F = make_recording_figures(Rec(r), modeColors, numPCs);
            if P.SaveFigures
                exportgraphics(F.winCompare, fullfile(figDir, [recLabel '_windowCompare.png']), 'Resolution', 200);
                exportgraphics(F.modeDynamics, fullfile(figDir, [recLabel '_modeDynamics.png']), 'Resolution', 200);
            end
            close_if_valid(F);
        end

    catch ME
        warning('Failed on folder:\n%s\n\n%s', recordings(r).folder, getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

% remove empty rows
ok = ~arrayfun(@isempty, {Rec.label});
Rec = Rec(ok);
rows = rows(ok);
if isempty(Rec)
    error('No recordings completed successfully.');
end

%% ========================= 4) BUILD TABLE =============================
summaryTable = struct2table(rows);

%% ========================= 5) POPULATION SUMMARIES ====================
Pop = build_population_summary(Rec, P.NumPCs);

if P.DoPlots
    Fpop = make_population_figures(Pop, summaryTable);
    if P.SaveFigures
        exportgraphics(Fpop.main, fullfile(popDir, 'population_window_comparison.png'), 'Resolution', 200);
        exportgraphics(Fpop.psd, fullfile(popDir, 'population_psd_comparison.png'), 'Resolution', 200);
    end
    close_if_valid(Fpop);
end

writetable(summaryTable, fullfile(outDir, 'summary_table.csv'));
save(fullfile(outDir, 'summary_struct.mat'), 'Rec', 'Pop', 'summaryTable', 'P', '-v7.3');

OUT = struct();
OUT.summaryTable = summaryTable;
OUT.population   = Pop;
OUT.recordings   = Rec;
OUT.cfg          = P;

end

%% =====================================================================
function [dff, fps] = extract_dff_and_fps(S)
dff = [];
fps = [];

vars = fieldnames(S);

% ---------- 1) direct candidate names ----------
directNames = {'dff_new','dff','DFF','dffMat','DV_DFFmovwindow'};

for i = 1:numel(directNames)
    vn = directNames{i};
    if isfield(S, vn)
        cand = S.(vn);
        if isnumeric(cand) && ismatrix(cand) && all(size(cand) > 10)
            dff = cand;
            break;
        end
    end
end

% ---------- 2) nested structs ----------
if isempty(dff)
    for i = 1:numel(vars)
        x = S.(vars{i});
        if ~isstruct(x), continue; end

        for j = 1:numel(directNames)
            vn = directNames{j};
            if isfield(x, vn)
                cand = x.(vn);
                if isnumeric(cand) && ismatrix(cand) && all(size(cand) > 10)
                    dff = cand;
                    break;
                end
            end
        end

        if ~isempty(dff), break; end
    end
end

% ---------- 3) fallback: choose the most plausible matrix ----------
if isempty(dff)
    bestScore = -inf;
    bestCand = [];

    for i = 1:numel(vars)
        cand = S.(vars{i});
        if ~isnumeric(cand) || ~ismatrix(cand)
            continue;
        end

        sz = size(cand);
        if min(sz) < 10 || max(sz) < 100
            continue;
        end

        % prefer long time-series matrices
        longDim = max(sz);
        shortDim = min(sz);

        % neurons x time should usually have many more frames than neurons
        ratioScore = longDim / max(shortDim,1);

        % penalize very square matrices
        score = ratioScore + 0.001*longDim;

        if score > bestScore
            bestScore = score;
            bestCand = cand;
        end
    end

    dff = bestCand;
end

if isempty(dff)
    error('Could not find a plausible DFF matrix.');
end

% ---------- orient as neurons x time ----------
if size(dff,1) > size(dff,2)
    dff = dff.';
end

% sanity check: time dimension should be reasonably long
if size(dff,2) < 1000
    warning('Extracted DFF has only %d frames. This is suspicious.', size(dff,2));
end

% ---------- fps ----------
candFps = {'fps','fs','Fs','origFps'};

for i = 1:numel(candFps)
    vn = candFps{i};
    if isfield(S, vn) && isnumeric(S.(vn)) && isscalar(S.(vn))
        fps = double(S.(vn));
        break;
    end
end

if isempty(fps)
    for i = 1:numel(vars)
        x = S.(vars{i});
        if ~isstruct(x), continue; end

        for j = 1:numel(candFps)
            vn = candFps{j};
            if isfield(x, vn) && isnumeric(x.(vn)) && isscalar(x.(vn))
                fps = double(x.(vn));
                break;
            end
        end

        if ~isempty(fps), break; end
    end
end

if isempty(fps)
    error('Could not find fps.');
end
end
%% =====================================================================
function [boutOnsetsFrames, tailTrace] = extract_bout_onsets_from_tail(S, nFrames, fps)
boutOnsetsFrames = [];
tailTrace = [];

if isfield(S, 'bou_ons_fra')
    boutOnsetsFrames = S.bou_ons_fra(:);
end
if isempty(boutOnsetsFrames)
    cand = {'bou_ons_fra','bout_onsets_frames','bout_frames','eventFrames'};
    for i = 1:numel(cand)
        if isfield(S, cand{i}) && isnumeric(S.(cand{i}))
            boutOnsetsFrames = S.(cand{i})(:);
            break;
        end
    end
end
if isfield(S, 'tai_ang_uni_frames')
    tailTrace = S.tai_ang_uni_frames(:);
elseif isfield(S, 'tail_trace')
    tailTrace = S.tail_trace(:);
end
if isempty(boutOnsetsFrames) && ~isempty(tailTrace)
    thr = median(abs(tailTrace), 'omitnan') + 3*mad(tailTrace,1);
    [~, locs] = findpeaks(abs(tailTrace), 'MinPeakHeight', thr, 'MinPeakDistance', max(1, round(0.2*fps)));
    boutOnsetsFrames = locs(:);
end
boutOnsetsFrames = unique(round(boutOnsetsFrames(:)));
boutOnsetsFrames = boutOnsetsFrames(isfinite(boutOnsetsFrames));
boutOnsetsFrames = boutOnsetsFrames(boutOnsetsFrames >= 1 & boutOnsetsFrames <= nFrames);
end

%% =====================================================================
function [winTable, lowIdx, highIdx] = rank_windows_by_bout_activity(boutFrames, nFrames, winFrames, stepFrames, fps)
starts = 1:stepFrames:(nFrames-winFrames+1);
if isempty(starts)
    winTable = table(); lowIdx = []; highIdx = [];
    return;
end
nW = numel(starts);
boutCount = zeros(nW,1);
for i = 1:nW
    s = starts(i);
    e = s + winFrames - 1;
    boutCount(i) = sum(boutFrames >= s & boutFrames <= e);
end
boutRatePerMin = boutCount ./ (winFrames/fps/60);
winTable = table(starts(:), starts(:)+winFrames-1, boutCount, boutRatePerMin, ...
    'VariableNames', {'startFrame','endFrame','boutCount','boutRatePerMin'});
[~, lowIdx]  = min(winTable.boutRatePerMin);
[~, highIdx] = max(winTable.boutRatePerMin);
end

%% =====================================================================
function X = fill_nans_rowwise(X)
for i = 1:size(X,1)
    xi = X(i,:);
    if any(~isfinite(xi))
        xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
        xi(~isfinite(xi)) = 0;
        X(i,:) = xi;
    end
end
end

%% =====================================================================
function Xz = zscore_fill_rows(X)
mu = mean(X,2,'omitnan');
sd = std(X,0,2,'omitnan');
sd(sd == 0 | ~isfinite(sd)) = 1;
Xz = (X - mu) ./ sd;
Xz(~isfinite(Xz)) = 0;
end

%% =====================================================================
function modeMean = compute_mode_mean(X, modeID, numPCs)
modeMean = nan(numPCs, size(X,2));
for k = 1:numPCs
    idx = modeID == k;
    if any(idx)
        modeMean(k,:) = mean(X(idx,:), 1, 'omitnan');
    end
end
end

%% =====================================================================
function Eval = eval_window_in_fixed_latent(latentTS_full, frameIdx, fps, P)
Xi = latentTS_full(frameIdx, :);
numPCs = size(Xi,2);

[pcDomFreq, pcPeakPower, fAxis, PxxPC] = compute_pc_psd(Xi, fps, P.WelchWinSec, P.WelchOverlapFrac, P.BandHz);

latentBand = bandpass_latent(Xi, fps, P.BandHz);
latentHilb = hilbert(latentBand);
latentPhase = angle(latentHilb);
latentAmp   = abs(latentHilb);

phaseMask = true(size(latentAmp));
for i = 1:numPCs
    ai = latentAmp(:,i);
    if all(~isfinite(ai))
        phaseMask(:,i) = false;
    else
        thr = quantile(ai(isfinite(ai)), P.PhaseAmpQuantileMin);
        phaseMask(:,i) = ai >= thr;
    end
end

PLV = nan(numPCs, numPCs);
for i = 1:numPCs
    for j = 1:numPCs
        good = isfinite(latentPhase(:,i)) & isfinite(latentPhase(:,j)) & phaseMask(:,i) & phaseMask(:,j);
        if nnz(good) < 20
            PLV(i,j) = NaN;
        else
            dphi = latentPhase(good,i) - latentPhase(good,j);
            PLV(i,j) = abs(mean(exp(1i*dphi)));
        end
    end
end

maskOff = triu(true(numPCs),1);
plvVals = PLV(maskOff);
plvVals = plvVals(isfinite(plvVals));
if isempty(plvVals)
    meanPLV = NaN;
else
    meanPLV = mean(plvVals);
end

[A, lambda, eigMag, eigAng, eigFreqHz, isCandidate, pairs, pairFreqHz, pairMag] = fit_latent_lds(Xi, fps, P);
pcFreqCV = std(pcDomFreq, 'omitnan') / mean(pcDomFreq, 'omitnan');

Eval = struct();
Eval.frameIdx             = frameIdx(:)';
Eval.latentTS             = Xi;
Eval.fAxis                = fAxis;
Eval.PxxPC                = PxxPC;
Eval.pcDomFreqHz          = pcDomFreq;
Eval.pcPeakPower          = pcPeakPower;
Eval.pcFreqCV             = pcFreqCV;
Eval.latentBand           = latentBand;
Eval.latentPhase          = latentPhase;
Eval.latentAmp            = latentAmp;
Eval.phaseMask            = phaseMask;
Eval.PLV                  = PLV;
Eval.meanPLV              = meanPLV;
Eval.A                    = A;
Eval.lambda               = lambda;
Eval.eigMag               = eigMag;
Eval.eigAng               = eigAng;
Eval.eigFreqHz            = eigFreqHz;
Eval.isCandidate          = isCandidate;
Eval.pairs                = pairs;
Eval.pairFreqHz           = pairFreqHz;
Eval.pairMag              = pairMag;
Eval.nOscillatorPairsLDS  = size(pairs,1);
end

%% =====================================================================
function frac = window_pc_variance_fraction(latentWin)
vars = var(latentWin, 0, 1, 'omitnan');
frac = vars ./ sum(vars, 'omitnan');
frac = frac(:)';
end

%% =====================================================================
function [A, lambda, eigMag, eigAng, eigFreqHz, isCandidate, pairs, pairFreqHz, pairMag] = fit_latent_lds(latentTS, fps, P)
X1 = latentTS(1:end-1, :)';
X2 = latentTS(2:end,   :)';
A = (X2 * X1') / (X1 * X1' + P.LDSReg * eye(size(X1,1)));
[~,D] = eig(A);
lambda = diag(D);
eigMag = abs(lambda);
eigAng = angle(lambda);
eigFreqHz = abs(eigAng) * fps / (2*pi);
isComplex = abs(imag(lambda)) > 1e-10;
isCandidate = isComplex & eigMag >= P.MinEigMag & eigMag <= P.MaxEigMag & ...
    eigFreqHz >= P.BandHz(1) & eigFreqHz <= P.BandHz(2);
oscIdx = find(isCandidate);
used = false(numel(oscIdx),1);
pairs = [];
for ii = 1:numel(oscIdx)
    if used(ii), continue; end
    lam = lambda(oscIdx(ii));
    rem = find(~used);
    remVals = lambda(oscIdx(rem));
    [~,jjrel] = min(abs(remVals - conj(lam)));
    jj = rem(jjrel);
    if jj ~= ii
        used(ii) = true; used(jj) = true;
        pairs = [pairs; oscIdx(ii) oscIdx(jj)]; %#ok<AGROW>
    else
        used(ii) = true;
    end
end
pairFreqHz = nan(size(pairs,1),1);
pairMag    = nan(size(pairs,1),1);
for i = 1:size(pairs,1)
    pairFreqHz(i) = mean(eigFreqHz(pairs(i,:)), 'omitnan');
    pairMag(i)    = mean(eigMag(pairs(i,:)), 'omitnan');
end
end

%% =====================================================================
function [pcDomFreq, pcPeakPower, fAxis, PxxPC] = compute_pc_psd(latentTS, fps, winSec, overlapFrac, bandHz)
nPC = size(latentTS,2);
nwin = max(16, round(winSec * fps));
nwin = min(nwin, size(latentTS,1));
nover = round(nwin * overlapFrac);
nfft = max(256, 2^nextpow2(nwin));
pcDomFreq = nan(nPC,1);
pcPeakPower = nan(nPC,1);
PxxPC = [];
fAxis = [];
for i = 1:nPC
    xi = latentTS(:,i);
    xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
    [Pxx, f] = pwelch(xi, hann(nwin), nover, nfft, fps);
    if isempty(PxxPC)
        PxxPC = nan(nPC, numel(Pxx));
        fAxis = f;
    end
    PxxPC(i,:) = Pxx(:)';
    mask = f >= bandHz(1) & f <= bandHz(2);
    if any(mask)
        [mx, idx] = max(Pxx(mask));
        fsub = f(mask);
        pcDomFreq(i) = fsub(idx);
        pcPeakPower(i) = mx;
    end
end
end

%% =====================================================================
function Xbp = bandpass_latent(X, fps, bandHz)
Wn = bandHz / (fps/2);
[b,a] = butter(3, Wn, 'bandpass');
Xbp = nan(size(X));
for i = 1:size(X,2)
    xi = fillmissing(X(:,i), 'linear', 'EndValues', 'nearest');
    Xbp(:,i) = filtfilt(b,a,xi);
end
end

%% =====================================================================
function F = make_recording_figures(R, modeColors, numPCs)
F = struct();

% ---------------- figure 1: FULL manifold + low/high window comparison ---
F.winCompare = figure('Color','w','Position',[60 60 1750 980], 'Name', ['Window compare: ' R.label]);
tlo = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

% trajectory panel with low/high slices on full manifold
ax1 = nexttile(tlo,[2 1]); hold(ax1,'on');
LT = R.fullPCA.latentTS;
plot3(ax1, LT(:,1), LT(:,2), LT(:,3), '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
idxL = R.lowFrames(1):R.lowFrames(2);
idxH = R.highFrames(1):R.highFrames(2);
plot3(ax1, LT(idxL,1), LT(idxL,2), LT(idxL,3), '-', 'Color', [0 0.45 0.85], 'LineWidth', 1.8);
plot3(ax1, LT(idxH,1), LT(idxH,2), LT(idxH,3), '-', 'Color', [0.85 0.25 0.1], 'LineWidth', 1.8);
scatter3(ax1, LT(R.boutOnsetsFrames,1), LT(R.boutOnsetsFrames,2), LT(R.boutOnsetsFrames,3), 10, [0.2 0.2 0.2], 'filled');
view(ax1, 40, 25); grid(ax1,'on');
xlabel(ax1,'PC1'); ylabel(ax1,'PC2'); zlabel(ax1,'PC3');
title(ax1, sprintf('Full manifold with low/high window slices\nPC1-3 = %.1f%%', sum(R.fullPCA.explained(1:min(3,end)))));
legend(ax1, {'Full','Low activity','High activity','Bout onsets'}, 'Location','best');

% explained variance
ax2 = nexttile(tlo); hold(ax2,'on');
plot(ax2, 1:numPCs, R.WinExpFull(1:numPCs)*100, '-o', 'Color', [0.3 0.3 0.3], 'LineWidth',1.4);
plot(ax2, 1:numPCs, R.WinExpLow(1:numPCs)*100, '-o', 'Color', [0 0.45 0.85], 'LineWidth',1.6);
plot(ax2, 1:numPCs, R.WinExpHigh(1:numPCs)*100, '-o', 'Color', [0.85 0.25 0.1], 'LineWidth',1.6);
xlabel(ax2,'PC'); ylabel(ax2,'% variance in window'); title(ax2,'Explained variance in fixed latent space');
legend(ax2, {'Full','Low','High'}, 'Location','best');

% PSD low vs high
ax3 = nexttile(tlo); hold(ax3,'on');
plot(ax3, R.EvalLow.fAxis, R.EvalLow.PxxPC(1,:), '-', 'Color', [0 0.45 0.85], 'LineWidth', 1.6);
plot(ax3, R.EvalHigh.fAxis, R.EvalHigh.PxxPC(1,:), '-', 'Color', [0.85 0.25 0.1], 'LineWidth', 1.6);
for k = 2:min(3,numPCs)
    plot(ax3, R.EvalLow.fAxis, R.EvalLow.PxxPC(k,:), '-', 'Color', [0 0.45 0.85 0.35], 'LineWidth', 0.8);
    plot(ax3, R.EvalHigh.fAxis, R.EvalHigh.PxxPC(k,:), '-', 'Color', [0.85 0.25 0.1 0.35], 'LineWidth', 0.8);
end
xlabel(ax3,'Frequency (Hz)'); ylabel(ax3,'Power'); title(ax3,'Latent PSDs: low vs high');
legend(ax3, {'Low PC1','High PC1'}, 'Location','best'); xlim(ax3,[0 0.15]);

% LDS eigenvalues low/high
ax4 = nexttile(tlo); hold(ax4,'on');
theta = linspace(0, 2*pi, 400);
plot(ax4, cos(theta), sin(theta), 'k--');
scatter(ax4, real(R.EvalLow.lambda),  imag(R.EvalLow.lambda), 30, [0 0.45 0.85], 'filled', 'MarkerFaceAlpha',0.55);
scatter(ax4, real(R.EvalHigh.lambda), imag(R.EvalHigh.lambda), 30, [0.85 0.25 0.1], 'filled', 'MarkerFaceAlpha',0.55);
axis(ax4,'equal'); grid(ax4,'on');
xlabel(ax4,'Real(\lambda)'); ylabel(ax4,'Imag(\lambda)');
title(ax4, sprintf('LDS eigenvalues\nlow pairs=%d | high pairs=%d', R.EvalLow.nOscillatorPairsLDS, R.EvalHigh.nOscillatorPairsLDS));

% PLV low
ax5 = nexttile(tlo);
imagesc(ax5, R.EvalLow.PLV, [0 1]); axis(ax5,'image'); colorbar(ax5);
xlabel(ax5,'PC'); ylabel(ax5,'PC'); title(ax5, sprintf('PLV | low | mean=%.3f', R.EvalLow.meanPLV));

% PLV high
ax6 = nexttile(tlo);
imagesc(ax6, R.EvalHigh.PLV, [0 1]); axis(ax6,'image'); colorbar(ax6);
xlabel(ax6,'PC'); ylabel(ax6,'PC'); title(ax6, sprintf('PLV | high | mean=%.3f', R.EvalHigh.meanPLV));

% text summary
ax7 = nexttile(tlo,[1 2]); axis(ax7,'off');
text(ax7,0.01,0.98,{ ...
    ['Recording: ' R.label], ...
    sprintf('fps = %.3f | neurons kept = %d', R.fps, R.nKeptShared), ...
    sprintf('Low window:  %d-%d | %.3f bouts/min', R.lowFrames(1), R.lowFrames(2), R.winTable.boutRatePerMin(find(R.winTable.startFrame==R.lowFrames(1),1))), ...
    sprintf('High window: %d-%d | %.3f bouts/min', R.highFrames(1), R.highFrames(2), R.winTable.boutRatePerMin(find(R.winTable.startFrame==R.highFrames(1),1))), ...
    sprintf('Low dom freqs (Hz):  %s', mat2str(R.EvalLow.pcDomFreqHz,4)), ...
    sprintf('High dom freqs (Hz): %s', mat2str(R.EvalHigh.pcDomFreqHz,4)), ...
    sprintf('Low mean PLV = %.3f | High mean PLV = %.3f', R.EvalLow.meanPLV, R.EvalHigh.meanPLV), ...
    sprintf('Low pcFreqCV = %.3f | High pcFreqCV = %.3f', R.EvalLow.pcFreqCV, R.EvalHigh.pcFreqCV), ...
    sprintf('Low LDS pairs = %d | High LDS pairs = %d', R.EvalLow.nOscillatorPairsLDS, R.EvalHigh.nOscillatorPairsLDS), ...
    sprintf('PC1 variance full/low/high = %.1f / %.1f / %.1f %%', 100*R.WinExpFull(1), 100*R.WinExpLow(1), 100*R.WinExpHigh(1)) ...
    }, 'VerticalAlignment','top','FontName','Courier');

% ---------------- figure 2: mode dynamics on full manifold ------------
F.modeDynamics = figure('Color','w','Position',[100 100 1700 260*numPCs], 'Name', ['Mode dynamics: ' R.label]);
tlo2 = tiledlayout(numPCs, 2, 'TileSpacing','compact','Padding','compact');
tSec = (0:size(R.modeMean,2)-1) / R.fps;
for k = 1:numPCs
    axL = nexttile(tlo2); hold(axL,'on');
    plot(axL, tSec, R.modeMean(k,:), 'Color', modeColors(k,:), 'LineWidth', 1.1);
    xline(axL, (R.lowFrames(1)-1)/R.fps, '--', 'Color', [0 0.45 0.85]);
    xline(axL, (R.lowFrames(2)-1)/R.fps, '--', 'Color', [0 0.45 0.85]);
    xline(axL, (R.highFrames(1)-1)/R.fps, '--', 'Color', [0.85 0.25 0.1]);
    xline(axL, (R.highFrames(2)-1)/R.fps, '--', 'Color', [0.85 0.25 0.1]);
    ylabel(axL, sprintf('mode %d', k), 'Color', modeColors(k,:));
    if k == 1, title(axL, 'Mode activity over time'); end
    if k == numPCs, xlabel(axL, 'time (s)'); end
    grid(axL,'on'); axis(axL,'tight');

    axR = nexttile(tlo2); hold(axR,'on');
    plot3(axR, LT(:,1), LT(:,2), LT(:,3), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.6);
    idx = R.domMode_t == k;
    plot3_colored_mask(axR, LT(:,1), LT(:,2), LT(:,3), idx, modeColors(k,:));
    if k == 1, title(axR, 'Latent trajectory colored by dominant mode'); end
    if k == numPCs
        xlabel(axR,'PC1'); ylabel(axR,'PC2'); zlabel(axR,'PC3');
    end
    grid(axR,'on'); view(axR,40,25);
end
end

%% =====================================================================
function plot3_colored_mask(ax, x, y, z, mask, col)
idx = find(mask(:));
if isempty(idx), return; end
breaks = [1; find(diff(idx) > 1)+1; numel(idx)+1];
for b = 1:numel(breaks)-1
    ii = idx(breaks(b):breaks(b+1)-1);
    if numel(ii) >= 2
        plot3(ax, x(ii), y(ii), z(ii), '-', 'Color', col, 'LineWidth', 1.2);
    else
        plot3(ax, x(ii), y(ii), z(ii), '.', 'Color', col, 'MarkerSize', 8);
    end
end
end

%% =====================================================================
function Pop = build_population_summary(Rec, maxNumPCs)
nR = numel(Rec);
Pop.nRecordings = nR;
Pop.lowPLV  = nan(nR,1); Pop.highPLV = nan(nR,1);
Pop.lowPairs = nan(nR,1); Pop.highPairs = nan(nR,1);
Pop.lowFreqCV = nan(nR,1); Pop.highFreqCV = nan(nR,1);
Pop.lowBoutRate = nan(nR,1); Pop.highBoutRate = nan(nR,1);
Pop.fullExp = nan(nR, maxNumPCs); Pop.lowExp = nan(nR, maxNumPCs); Pop.highExp = nan(nR, maxNumPCs);
Pop.lowDomFreq = nan(nR, maxNumPCs); Pop.highDomFreq = nan(nR, maxNumPCs);
for r = 1:nR
    Pop.lowPLV(r)   = Rec(r).EvalLow.meanPLV;
    Pop.highPLV(r)  = Rec(r).EvalHigh.meanPLV;
    Pop.lowPairs(r) = Rec(r).EvalLow.nOscillatorPairsLDS;
    Pop.highPairs(r)= Rec(r).EvalHigh.nOscillatorPairsLDS;
    Pop.lowFreqCV(r)= Rec(r).EvalLow.pcFreqCV;
    Pop.highFreqCV(r)= Rec(r).EvalHigh.pcFreqCV;
    iLow  = find(Rec(r).winTable.startFrame == Rec(r).lowFrames(1),1);
    iHigh = find(Rec(r).winTable.startFrame == Rec(r).highFrames(1),1);
    Pop.lowBoutRate(r)  = Rec(r).winTable.boutRatePerMin(iLow);
    Pop.highBoutRate(r) = Rec(r).winTable.boutRatePerMin(iHigh);
    k1 = min(maxNumPCs, numel(Rec(r).WinExpFull));
    Pop.fullExp(r,1:k1) = Rec(r).WinExpFull(1:k1);
    k2 = min(maxNumPCs, numel(Rec(r).WinExpLow));
    Pop.lowExp(r,1:k2) = Rec(r).WinExpLow(1:k2);
    k3 = min(maxNumPCs, numel(Rec(r).WinExpHigh));
    Pop.highExp(r,1:k3) = Rec(r).WinExpHigh(1:k3);
    kdL = min(maxNumPCs, numel(Rec(r).EvalLow.pcDomFreqHz));
    kdH = min(maxNumPCs, numel(Rec(r).EvalHigh.pcDomFreqHz));
    Pop.lowDomFreq(r,1:kdL) = Rec(r).EvalLow.pcDomFreqHz(1:kdL);
    Pop.highDomFreq(r,1:kdH)= Rec(r).EvalHigh.pcDomFreqHz(1:kdH);
end
end

%% =====================================================================
function F = make_population_figures(Pop, summaryTable)
F = struct();

F.main = figure('Color','w','Position',[80 80 1600 900], 'Name','Population low vs high activity');
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tlo); hold(ax1,'on');
boxchart(ax1, [ones(size(Pop.lowPLV)); 2*ones(size(Pop.highPLV))], [Pop.lowPLV; Pop.highPLV]);
plot(ax1, [ones(size(Pop.lowPLV)) 2*ones(size(Pop.highPLV))]', [Pop.lowPLV Pop.highPLV]', '-k', 'Color', [0.7 0.7 0.7]);
scatter(ax1, ones(size(Pop.lowPLV)), Pop.lowPLV, 22, [0 0.45 0.85], 'filled');
scatter(ax1, 2*ones(size(Pop.highPLV)), Pop.highPLV, 22, [0.85 0.25 0.1], 'filled');
set(ax1,'XTick',[1 2],'XTickLabel',{'Low','High'}); ylabel(ax1,'Mean PLV'); title(ax1,'Phase-locking');

ax2 = nexttile(tlo); hold(ax2,'on');
boxchart(ax2, [ones(size(Pop.lowPairs)); 2*ones(size(Pop.highPairs))], [Pop.lowPairs; Pop.highPairs]);
plot(ax2, [ones(size(Pop.lowPairs)) 2*ones(size(Pop.highPairs))]', [Pop.lowPairs Pop.highPairs]', '-k', 'Color', [0.7 0.7 0.7]);
scatter(ax2, ones(size(Pop.lowPairs)), Pop.lowPairs, 22, [0 0.45 0.85], 'filled');
scatter(ax2, 2*ones(size(Pop.highPairs)), Pop.highPairs, 22, [0.85 0.25 0.1], 'filled');
set(ax2,'XTick',[1 2],'XTickLabel',{'Low','High'}); ylabel(ax2,'# LDS oscillator pairs'); title(ax2,'LDS candidate pairs');

ax3 = nexttile(tlo); hold(ax3,'on');
boxchart(ax3, [ones(size(Pop.lowFreqCV)); 2*ones(size(Pop.highFreqCV))], [Pop.lowFreqCV; Pop.highFreqCV]);
plot(ax3, [ones(size(Pop.lowFreqCV)) 2*ones(size(Pop.highFreqCV))]', [Pop.lowFreqCV Pop.highFreqCV]', '-k', 'Color', [0.7 0.7 0.7]);
scatter(ax3, ones(size(Pop.lowFreqCV)), Pop.lowFreqCV, 22, [0 0.45 0.85], 'filled');
scatter(ax3, 2*ones(size(Pop.highFreqCV)), Pop.highFreqCV, 22, [0.85 0.25 0.1], 'filled');
set(ax3,'XTick',[1 2],'XTickLabel',{'Low','High'}); ylabel(ax3,'PC frequency CV'); title(ax3,'Timescale spread');

ax4 = nexttile(tlo); hold(ax4,'on');
meanLow = mean(Pop.lowExp,1,'omitnan')*100;
meanHigh = mean(Pop.highExp,1,'omitnan')*100;
plot(ax4, 1:numel(meanLow), meanLow, '-o', 'Color', [0 0.45 0.85], 'LineWidth', 1.6);
plot(ax4, 1:numel(meanHigh), meanHigh, '-o', 'Color', [0.85 0.25 0.1], 'LineWidth', 1.6);
xlabel(ax4,'PC'); ylabel(ax4,'Mean % variance'); title(ax4,'Windowed explained variance'); legend(ax4, {'Low','High'}, 'Location','best');

ax5 = nexttile(tlo); hold(ax5,'on');
imagesc(ax5, meanLow(:)'); colormap(ax5, parula); colorbar(ax5); title(ax5,'Mean low explained variance'); xlabel(ax5,'PC'); set(ax5,'YTick',[]);

ax6 = nexttile(tlo); hold(ax6,'on');
imagesc(ax6, meanHigh(:)'); colormap(ax6, parula); colorbar(ax6); title(ax6,'Mean high explained variance'); xlabel(ax6,'PC'); set(ax6,'YTick',[]);

F.psd = figure('Color','w','Position',[120 120 1400 800], 'Name','Population PSD / frequencies');
tlo2 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax7 = nexttile(tlo2); hold(ax7,'on');
boxchart(ax7, [repmat((1:size(Pop.lowDomFreq,2))', size(Pop.lowDomFreq,1), 1); repmat((1:size(Pop.highDomFreq,2))'+0.35, size(Pop.highDomFreq,1),1)], ...
    [Pop.lowDomFreq(:); Pop.highDomFreq(:)]);
scatter(ax7, repmat((1:size(Pop.lowDomFreq,2))', size(Pop.lowDomFreq,1), 1), Pop.lowDomFreq(:), 8, [0 0.45 0.85], 'filled');
scatter(ax7, repmat((1:size(Pop.highDomFreq,2))'+0.35, size(Pop.highDomFreq,1), 1), Pop.highDomFreq(:), 8, [0.85 0.25 0.1], 'filled');
xlabel(ax7,'PC'); ylabel(ax7,'Dominant frequency (Hz)'); title(ax7,'Low vs high dominant frequencies');

ax8 = nexttile(tlo2); hold(ax8,'on');
scatter(ax8, Pop.lowBoutRate, Pop.lowPLV, 35, [0 0.45 0.85], 'filled');
scatter(ax8, Pop.highBoutRate, Pop.highPLV, 35, [0.85 0.25 0.1], 'filled');
xlabel(ax8,'Bout rate (bouts/min)'); ylabel(ax8,'Mean PLV'); title(ax8,'Activity vs phase-locking'); legend(ax8, {'Low','High'}, 'Location','best');

ax9 = nexttile(tlo2); hold(ax9,'on');
scatter(ax9, Pop.lowBoutRate, Pop.lowPairs, 35, [0 0.45 0.85], 'filled');
scatter(ax9, Pop.highBoutRate, Pop.highPairs, 35, [0.85 0.25 0.1], 'filled');
xlabel(ax9,'Bout rate (bouts/min)'); ylabel(ax9,'# LDS pairs'); title(ax9,'Activity vs LDS pairs');

ax10 = nexttile(tlo2); axis(ax10,'off');
text(ax10, 0.01, 0.98, { ...
    sprintf('n recordings = %d', height(summaryTable)), ...
    sprintf('Mean low bout rate  = %.3f', mean(Pop.lowBoutRate,'omitnan')), ...
    sprintf('Mean high bout rate = %.3f', mean(Pop.highBoutRate,'omitnan')), ...
    sprintf('Mean low PLV  = %.3f', mean(Pop.lowPLV,'omitnan')), ...
    sprintf('Mean high PLV = %.3f', mean(Pop.highPLV,'omitnan')), ...
    sprintf('Mean low PCfreqCV  = %.3f', mean(Pop.lowFreqCV,'omitnan')), ...
    sprintf('Mean high PCfreqCV = %.3f', mean(Pop.highFreqCV,'omitnan')), ...
    sprintf('Mean low pairs  = %.3f', mean(Pop.lowPairs,'omitnan')), ...
    sprintf('Mean high pairs = %.3f', mean(Pop.highPairs,'omitnan')) ...
    }, 'VerticalAlignment','top','FontName','Courier');
end

%% =====================================================================
function label = make_recording_label(folderPath)
[~, parent] = fileparts(folderPath);
label = parent;
end

%% =====================================================================
function close_if_valid(S)
fn = fieldnames(S);
for i = 1:numel(fn)
    try
        if isgraphics(S.(fn{i}))
            close(S.(fn{i}));
        end
    catch
    end
end
end
