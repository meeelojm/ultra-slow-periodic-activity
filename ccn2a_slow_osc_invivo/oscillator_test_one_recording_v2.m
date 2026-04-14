function OUT = oscillator_test_one_recording_v2(dff_full, fps, varargin)
% OSCILLATOR_TEST_ONE_RECORDING_V2
%
% Purpose
% -------
% Analyze one recording at the population level to test whether the latent
% dynamics are more consistent with:
%   - one dominant oscillator,
%   - multiple oscillatory modes,
%   - or weak / overlapping / nonstationary slow modes.
%
% This version also:
%   - optionally filters neurons using ACF peak prominence / peak count
%   - assigns neurons to latent modes using PCA loadings
%   - produces cleaner mode figures:
%       1) neuron membership heatmap + color strip
%       2) stacked per-mode activity traces + one colored manifold plot
%       3) optional anatomical scatter3 plots if neuron positions are provided
%
% Inputs
% ------
% dff_full : [neurons x time]
% fps      : scalar
%
% Name-value parameters
% ---------------------
% 'EpochFrames'            : [start end] frames, default = full recording
% 'EventFrames'            : event frames in absolute frame index
% 'FiniteFracThresh'       : minimum fraction of finite samples per neuron
% 'UseZScore'              : z-score neurons before PCA
% 'NumPCs'                 : number of latent PCs to use
% 'BandHz'                 : band for phase analyses
% 'WelchWinSec'            : PSD Welch window length in sec
% 'WelchOverlapFrac'       : PSD overlap fraction
% 'LDSReg'                 : ridge regularization for LDS fit
% 'MinEigMag'              : minimum eigenvalue magnitude for oscillator candidate
% 'MaxEigMag'              : maximum eigenvalue magnitude for oscillator candidate
% 'PhaseAmpQuantileMin'    : mask lowest-amplitude fraction for phase usage
% 'PLVThreshold'           : threshold for phase-lock graph
% 'WinSec'                 : sliding window length in sec
% 'StepSec'                : sliding window step in sec
%
% ACF filter parameters
% ---------------------
% 'EnableAcfPeakFilter'    : true/false
% 'AcfPeakLagRangeSec'     : lag range where peaks are counted
% 'AcfMinProminence'       : MinPeakProminence for ACF peaks
% 'MinNumAcPeaks'          : minimum number of ACF peaks to keep neuron
% 'AcfMaxLagSec'           : maximum lag for ACF calculation
%
% Plotting
% --------
% 'MakeModeMembershipPlot' : true/false
% 'NeuronPositions'        : optional [neurons x 3] or [neurons x 2]
% 'DoPlots'                : true/false
% 'Verbose'                : true/false
%
% Output
% ------
% OUT struct with:
%   .data
%   .filter
%   .pca
%   .mode
%   .psd
%   .phase
%   .lds
%   .sliding
%   .summary
%   .fig
%
% ---------------------------------------------------------------------

%% ============================= PARSE ===================================
ip = inputParser;
addRequired(ip, 'dff_full', @(x) isnumeric(x) && ndims(x)==2);
addRequired(ip, 'fps', @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'EpochFrames', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'EventFrames', [], @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'FiniteFracThresh', 0.95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(ip, 'UseZScore', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'NumPCs', 6, @(x) isnumeric(x) && isscalar(x) && x >= 2);

addParameter(ip, 'BandHz', [0.01 0.10], @(x) isnumeric(x) && numel(x)==2 && x(1)>0 && x(2)>x(1));
addParameter(ip, 'WelchWinSec', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'WelchOverlapFrac', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);

addParameter(ip, 'LDSReg', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(ip, 'MinEigMag', 0.90, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'MaxEigMag', 1.05, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'PhaseAmpQuantileMin', 0.20, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
addParameter(ip, 'PLVThreshold', 0.70, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

addParameter(ip, 'WinSec', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'StepSec', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);

% ACF peak-prominence filter
addParameter(ip, 'EnableAcfPeakFilter', false, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'AcfPeakLagRangeSec', [10 120], @(x) isnumeric(x) && numel(x)==2 && x(1)>=0 && x(2)>x(1));
addParameter(ip, 'AcfMinProminence', 0.08, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(ip, 'MinNumAcPeaks', 2, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(ip, 'AcfMaxLagSec', 180, @(x) isnumeric(x) && isscalar(x) && x > 0);

% mode plots / anatomy
addParameter(ip, 'MakeModeMembershipPlot', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'NeuronPositions', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x) && (size(x,2)==2 || size(x,2)==3)));

addParameter(ip, 'DoPlots', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(ip, dff_full, fps, varargin{:});
P = ip.Results;

[N0, T0] = size(dff_full);

%% ============================= EPOCH ===================================
if isempty(P.EpochFrames)
    epochFrames = [1 T0];
else
    epochFrames = round(P.EpochFrames(:)');
    epochFrames(1) = max(1, epochFrames(1));
    epochFrames(2) = min(T0, epochFrames(2));
end

tIdx = epochFrames(1):epochFrames(2);
dff_epoch = dff_full(:, tIdx);
Te = numel(tIdx);

%% ============================= EVENTS ==================================
eventFrames = round(P.EventFrames(:));
eventFrames = eventFrames(isfinite(eventFrames));
eventFrames = eventFrames(eventFrames >= epochFrames(1) & eventFrames <= epochFrames(2));
eventFramesLocal = eventFrames - epochFrames(1) + 1;
eventTimesSec = (eventFramesLocal - 1) / fps;

%% ============================= BASE QC =================================
fracFinite = mean(isfinite(dff_epoch), 2);
keepBase = fracFinite >= P.FiniteFracThresh;

%% ======================== OPTIONAL ACF FILTER ===========================
nAcPeaks_all = nan(N0,1);
keepAcf = true(N0,1);

if P.EnableAcfPeakFilter
    maxLag = round(P.AcfMaxLagSec * fps);
    lagAxis = (-maxLag:maxLag) / fps;
    lagMask = lagAxis >= P.AcfPeakLagRangeSec(1) & lagAxis <= P.AcfPeakLagRangeSec(2);

    for i = 1:N0
        if ~keepBase(i)
            keepAcf(i) = false;
            nAcPeaks_all(i) = NaN;
            continue;
        end

        xi = dff_epoch(i,:);
        xi = xi(isfinite(xi));

        if numel(xi) < maxLag + 5
            keepAcf(i) = false;
            nAcPeaks_all(i) = 0;
            continue;
        end

        xi = xi - mean(xi);
        aci = xcorr(xi, maxLag, 'coeff');
        aciWin = aci(lagMask);

        if all(~isfinite(aciWin))
            nAcPeaks_all(i) = 0;
            keepAcf(i) = false;
            continue;
        end

        pks = findpeaks(aciWin, 'MinPeakProminence', P.AcfMinProminence);
        nAcPeaks_all(i) = numel(pks);
        keepAcf(i) = nAcPeaks_all(i) >= P.MinNumAcPeaks;
    end
else
    nAcPeaks_all(:) = NaN;
end

keep = keepBase & keepAcf;
keepIdx = find(keep);
dropIdx = find(~keep);

if isempty(keepIdx)
    error('No neurons survived filtering.');
end

Xraw = dff_epoch(keepIdx, :);

if P.UseZScore
    X = zscore_fill_rows(Xraw);
else
    X = Xraw;
    X(~isfinite(X)) = 0;
end

% positions (optional)
positionsKept = [];
if ~isempty(P.NeuronPositions)
    pos = P.NeuronPositions;
    if size(pos,1) ~= N0
        warning('NeuronPositions has %d rows but dff_full has %d neurons. Ignoring positions.', size(pos,1), N0);
    else
        if size(pos,2) == 2
            pos = [pos, zeros(size(pos,1),1)];
        end
        positionsKept = pos(keepIdx, :);
    end
end

if P.Verbose
    fprintf('\n=== oscillator_test_one_recording_v2 ===\n');
    fprintf('Input neurons               : %d\n', N0);
    fprintf('Epoch frames                : %d to %d (%d frames)\n', epochFrames(1), epochFrames(2), Te);
    fprintf('Kept neurons                : %d\n', numel(keepIdx));
    fprintf('Band for phase analyses     : [%.3f %.3f] Hz\n', P.BandHz(1), P.BandHz(2));
    fprintf('Stim events in epoch        : %d\n', numel(eventFrames));
    if P.EnableAcfPeakFilter
        fprintf('ACF peak filter ON          : lag=[%.1f %.1f] s | prom=%.3f | min peaks=%d\n', ...
            P.AcfPeakLagRangeSec(1), P.AcfPeakLagRangeSec(2), P.AcfMinProminence, P.MinNumAcPeaks);
        fprintf('Neurons kept after ACF filt : %d / %d\n', nnz(keep), N0);
    end
end

%% ============================== PCA ====================================
maxPC = min([P.NumPCs, rank(X'), size(X,1), size(X,2)]);
[coeff, score, latent, ~, explained, mu] = pca(X', 'NumComponents', maxPC);
numPCs = maxPC;
latentTS = score(:, 1:numPCs);      % time x PCs
coeffUse = coeff(:, 1:numPCs);      % kept neurons x PCs

%% ======================= MODE ASSIGNMENT ===============================
% Assign each neuron to the PC on which it loads most strongly
[modeLoadingAbs, modeID] = max(abs(coeffUse), [], 2);
modeSign = sign(coeffUse(sub2ind(size(coeffUse), (1:size(coeffUse,1))', modeID)));

% Sort neurons by mode, then by loading strength within mode
sortTable = [modeID(:), -modeLoadingAbs(:)];
[~, sortOrd] = sortrows(sortTable, [1 2]);
modeID_sorted = modeID(sortOrd);

% Mean activity per mode
modeMean = nan(numPCs, Te);
modeCount = zeros(numPCs,1);

for k = 1:numPCs
    idx = (modeID == k);
    modeCount(k) = nnz(idx);
    if any(idx)
        modeMean(k,:) = mean(X(idx,:), 1, 'omitnan');
    end
end

% Dominant mode over time from mode mean activity
[~, domMode_t] = max(abs(modeMean), [], 1);

%% =========================== LATENT PSD ================================
[pcDomFreq, pcPeakPower, fAxis, PxxPC] = compute_pc_psd( ...
    latentTS, fps, P.WelchWinSec, P.WelchOverlapFrac, P.BandHz);

%% ======================= BANDPASS + HILBERT ============================
latentBand = bandpass_latent(latentTS, fps, P.BandHz);
latentHilb = hilbert(latentBand);
latentPhase = angle(latentHilb);
latentAmp = abs(latentHilb);

% Mask low-amplitude phase samples
phaseMask = true(size(latentAmp));
for i = 1:numPCs
    ai = latentAmp(:,i);
    goodAi = ai(isfinite(ai));
    if isempty(goodAi)
        phaseMask(:,i) = false(size(ai));
        continue;
    end
    thr = quantile(goodAi, P.PhaseAmpQuantileMin);
    phaseMask(:,i) = ai >= thr;
end

PLV = nan(numPCs, numPCs);
for i = 1:numPCs
    for j = 1:numPCs
        good = isfinite(latentPhase(:,i)) & isfinite(latentPhase(:,j)) & ...
               phaseMask(:,i) & phaseMask(:,j);
        if nnz(good) < 20
            PLV(i,j) = NaN;
        else
            dphi = latentPhase(good,i) - latentPhase(good,j);
            PLV(i,j) = abs(mean(exp(1i * dphi)));
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

adj = false(numPCs, numPCs);
goodPLV = isfinite(PLV);
adj(goodPLV) = PLV(goodPLV) >= P.PLVThreshold;
adj(1:numPCs+1:end) = false;
nPhaseComponents = count_connected_components(adj);

%% =============================== LDS ===================================
[A, lambda, eigMag, eigAng, eigFreqHz, isCandidate, pairs, pairFreqHz, pairMag] = ...
    fit_latent_lds(latentTS, fps, P);

pcFreqCV = std(pcDomFreq, 'omitnan') / mean(pcDomFreq, 'omitnan');

%% ============================= SLIDING =================================
sliding = run_sliding_oscillator_analysis(latentTS, fps, eventFramesLocal, P);

%% ============================= SUMMARY =================================
summary = struct();
summary.nKeptNeurons        = numel(keepIdx);
summary.numPCs              = numPCs;
summary.pcDomFreqHz         = pcDomFreq;
summary.pcFreqCV            = pcFreqCV;
summary.meanPLV             = meanPLV;
summary.nPhaseComponents    = nPhaseComponents;
summary.nOscillatorPairsLDS = size(pairs,1);
summary.ldsPairFreqHz       = pairFreqHz;
summary.ldsPairMagnitude    = pairMag;
summary.meanSlidingPairs    = mean(sliding.nOscPairs, 'omitnan');
summary.maxSlidingPairs     = max(sliding.nOscPairs, [], 'omitnan');
summary.modeCount           = modeCount;

if size(pairs,1) <= 1 && isfinite(meanPLV) && meanPLV > 0.8
    summary.interpretation = "More consistent with one dominant oscillator + noise.";
elseif size(pairs,1) >= 2 && (pcFreqCV > 0.15 || (~isnan(meanPLV) && meanPLV < 0.7) || nPhaseComponents >= 2)
    summary.interpretation = "More consistent with multiple oscillatory modes.";
elseif summary.maxSlidingPairs >= 2
    summary.interpretation = "Global fit is mixed, but sliding-window dynamics suggest transient multiple oscillatory modes.";
else
    summary.interpretation = "Mixed evidence; possible weakly coupled or partially overlapping oscillators.";
end

if P.Verbose
    fprintf('Latent PC dominant freqs (Hz): %s\n', mat2str(pcDomFreq, 4));
    fprintf('PC dominant frequency CV     : %.4f\n', pcFreqCV);
    fprintf('Mean off-diagonal PLV        : %.4f\n', meanPLV);
    fprintf('Phase-lock components        : %d (threshold %.2f)\n', nPhaseComponents, P.PLVThreshold);
    fprintf('LDS oscillatory eigen-pairs  : %d\n', size(pairs,1));
    fprintf('Mean sliding-window pairs    : %.3f\n', summary.meanSlidingPairs);
    fprintf('Max sliding-window pairs     : %.3f\n', summary.maxSlidingPairs);
    if ~isempty(pairFreqHz)
        fprintf('LDS pair freqs (Hz)          : %s\n', mat2str(pairFreqHz, 4));
    end
    fprintf('Interpretation               : %s\n', summary.interpretation);
end

%% ============================== OUTPUT =================================
OUT = struct();
OUT.params = P;

OUT.data = struct();
OUT.data.keepIdx          = keepIdx;
OUT.data.dropIdx          = dropIdx;
OUT.data.epochFrames      = epochFrames;
OUT.data.eventFramesAbs   = eventFrames;
OUT.data.eventFramesLocal = eventFramesLocal;
OUT.data.eventTimesSec    = eventTimesSec;
OUT.data.Xraw             = Xraw;
OUT.data.Xz               = X;
OUT.data.positions3D      = positionsKept;

OUT.filter = struct();
OUT.filter.keepBase       = keepBase;
OUT.filter.keepAcf        = keepAcf;
OUT.filter.keepFinal      = keep;
OUT.filter.nAcPeaks_all   = nAcPeaks_all;

OUT.pca = struct();
OUT.pca.coeff       = coeff;
OUT.pca.score       = score;
OUT.pca.latentTS    = latentTS;
OUT.pca.latent      = latent;
OUT.pca.explained   = explained;
OUT.pca.mu          = mu;

OUT.mode = struct();
OUT.mode.modeID          = modeID;
OUT.mode.modeLoadingAbs  = modeLoadingAbs;
OUT.mode.modeSign        = modeSign;
OUT.mode.sortOrd         = sortOrd;
OUT.mode.modeID_sorted   = modeID_sorted;
OUT.mode.modeMean        = modeMean;
OUT.mode.domMode_t       = domMode_t;
OUT.mode.modeCount       = modeCount;

OUT.psd = struct();
OUT.psd.fAxis        = fAxis;
OUT.psd.PxxPC        = PxxPC;
OUT.psd.pcDomFreqHz  = pcDomFreq;
OUT.psd.pcPeakPower  = pcPeakPower;

OUT.phase = struct();
OUT.phase.latentBand  = latentBand;
OUT.phase.latentHilb  = latentHilb;
OUT.phase.latentPhase = latentPhase;
OUT.phase.latentAmp   = latentAmp;
OUT.phase.phaseMask   = phaseMask;
OUT.phase.PLV         = PLV;

OUT.lds = struct();
OUT.lds.A           = A;
OUT.lds.lambda      = lambda;
OUT.lds.eigMag      = eigMag;
OUT.lds.eigAng      = eigAng;
OUT.lds.eigFreqHz   = eigFreqHz;
OUT.lds.isCandidate = isCandidate;
OUT.lds.pairs       = pairs;
OUT.lds.pairFreqHz  = pairFreqHz;
OUT.lds.pairMag     = pairMag;

OUT.sliding = sliding;
OUT.summary = summary;
OUT.fig = struct();

if P.DoPlots
    OUT.fig = make_all_plots(OUT, fps);
end

end

%% =======================================================================
function Xz = zscore_fill_rows(X)
mu = mean(X, 2, 'omitnan');
sd = std(X, 0, 2, 'omitnan');
sd(sd == 0 | ~isfinite(sd)) = 1;
Xz = (X - mu) ./ sd;
Xz(~isfinite(Xz)) = 0;
end

%% =======================================================================
function [A, lambda, eigMag, eigAng, eigFreqHz, isCandidate, pairs, pairFreqHz, pairMag] = fit_latent_lds(latentTS, fps, P)
X1 = latentTS(1:end-1, :)';
X2 = latentTS(2:end,   :)';
reg = P.LDSReg;

A = (X2 * X1') / (X1 * X1' + reg * eye(size(X1,1)));

[~, D] = eig(A);
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
    k = oscIdx(ii);
    lam = lambda(k);

    rem = find(~used);
    remVals = lambda(oscIdx(rem));
    [~, jjrel] = min(abs(remVals - conj(lam)));
    jj = rem(jjrel);

    if jj ~= ii
        used(ii) = true;
        used(jj) = true;
        pairs = [pairs; oscIdx(ii) oscIdx(jj)]; %#ok<AGROW>
    else
        used(ii) = true;
    end
end

pairFreqHz = nan(size(pairs,1),1);
pairMag = nan(size(pairs,1),1);
for i = 1:size(pairs,1)
    pairFreqHz(i) = mean(eigFreqHz(pairs(i,:)), 'omitnan');
    pairMag(i) = mean(eigMag(pairs(i,:)), 'omitnan');
end
end

%% =======================================================================
function sliding = run_sliding_oscillator_analysis(latentTS, fps, eventFramesLocal, P)
T = size(latentTS,1);
numPCs = size(latentTS,2);

winF = max(20, round(P.WinSec * fps));
stepF = max(1, round(P.StepSec * fps));

startFrames = 1:stepF:max(1, T-winF+1);
nWin = numel(startFrames);

sliding = struct();
sliding.startFrames   = startFrames(:);
sliding.centerFrames  = startFrames(:) + floor(winF/2);
sliding.centerSec     = (sliding.centerFrames - 1) / fps;
sliding.nOscPairs     = nan(nWin,1);
sliding.pcDomFreqHz   = nan(nWin, numPCs);
sliding.nEvents       = zeros(nWin,1);

for w = 1:nWin
    s = startFrames(w);
    e = min(T, s + winF - 1);
    Xi = latentTS(s:e, :);

    if ~isempty(eventFramesLocal)
        sliding.nEvents(w) = sum(eventFramesLocal >= s & eventFramesLocal <= e);
    end

    [pcDomFreq, ~, ~, ~] = compute_pc_psd(Xi, fps, min(P.WelchWinSec, (e-s+1)/fps), P.WelchOverlapFrac, P.BandHz);
    sliding.pcDomFreqHz(w, 1:numel(pcDomFreq)) = pcDomFreq(:)';

    if size(Xi,1) >= 20
        [~, ~, ~, ~, ~, ~, pairs, ~, ~] = fit_latent_lds(Xi, fps, P);
        sliding.nOscPairs(w) = size(pairs,1);
    end
end
end

%% =======================================================================
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

%% =======================================================================
function Xbp = bandpass_latent(X, fps, bandHz)
Wn = bandHz / (fps/2);
if Wn(2) >= 1
    error('Upper band edge must be below Nyquist.');
end

[b, a] = butter(3, Wn, 'bandpass');
Xbp = nan(size(X));

for i = 1:size(X,2)
    xi = X(:,i);
    xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
    Xbp(:,i) = filtfilt(b, a, xi);
end
end

%% =======================================================================
function nComp = count_connected_components(adj)
n = size(adj,1);
visited = false(n,1);
nComp = 0;

for i = 1:n
    if visited(i), continue; end
    nComp = nComp + 1;
    stack = i;
    while ~isempty(stack)
        v = stack(end);
        stack(end) = [];
        if visited(v), continue; end
        visited(v) = true;
        nbrs = find(adj(v,:) | adj(:,v)');
        nbrs = nbrs(~visited(nbrs));
        stack = [stack, nbrs]; %#ok<AGROW>
    end
end
end

%% =======================================================================
function figs = make_all_plots(OUT, fps)
figs = struct();

sumFigs = make_summary_plots(OUT, fps);
figs.summary = sumFigs.summary;
figs.phase   = sumFigs.phase;
figs.sliding = sumFigs.sliding;

if OUT.params.MakeModeMembershipPlot
    modeFigs = make_mode_plots(OUT, fps);
    figs.modeMembership = modeFigs.modeMembership;
    figs.modeDynamics   = modeFigs.modeDynamics;
    if isfield(modeFigs, 'modeBrain')
        figs.modeBrain = modeFigs.modeBrain;
    end
end
end

%% =======================================================================
function figs = make_summary_plots(OUT, fps)
figs = struct();

latentTS     = OUT.pca.latentTS;
explained    = OUT.pca.explained;
fAxis        = OUT.psd.fAxis;
PxxPC        = OUT.psd.PxxPC;
pcDomFreq    = OUT.psd.pcDomFreqHz;
PLV          = OUT.phase.PLV;
lambda       = OUT.lds.lambda;
isCand       = OUT.lds.isCandidate;
pairs        = OUT.lds.pairs;
eventTimesSec = OUT.data.eventTimesSec;

numPCs = size(latentTS,2);
tSec   = (0:size(latentTS,1)-1) / fps;

% ===================== SUMMARY FIGURE =====================
figs.summary = figure('Color','w','Position',[80 60 1600 950], 'Name','Oscillator test summary');
tlo = tiledlayout(2,4, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile(tlo, [2 1]); hold(ax1,'on');
if size(latentTS,2) >= 3
    plot3(ax1, latentTS(:,1), latentTS(:,2), latentTS(:,3), '-', ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    if ~isempty(eventTimesSec)
        evIdx = round(eventTimesSec * fps) + 1;
        evIdx = evIdx(evIdx>=1 & evIdx<=size(latentTS,1));
        scatter3(ax1, latentTS(evIdx,1), latentTS(evIdx,2), latentTS(evIdx,3), ...
            30, 'r', 'filled');
    end
    xlabel(ax1,'PC1'); ylabel(ax1,'PC2'); zlabel(ax1,'PC3');
    title(ax1, sprintf('PC trajectory (PC1-3 = %.1f%%)', sum(explained(1:min(3,end)))));
    view(ax1,45,30); grid(ax1,'on'); box(ax1,'off');
end

ax2 = nexttile(tlo);
bar(ax2, explained(1:min(10,numel(explained))));
xlabel(ax2,'PC');
ylabel(ax2,'% variance explained');
title(ax2,'Explained variance');

ax3 = nexttile(tlo); hold(ax3,'on');
for i = 1:numPCs
    plot(ax3, fAxis, PxxPC(i,:), 'LineWidth', 1.2);
end
xlim(ax3, [0 max(OUT.params.BandHz(2)*1.5, 0.12)]);
xlabel(ax3,'Frequency (Hz)');
ylabel(ax3,'Power');
title(ax3,'PSD of latent PCs');
legtxt = arrayfun(@(i) sprintf('PC%d (%.3f Hz)', i, pcDomFreq(i)), ...
    1:numPCs, 'UniformOutput', false);
legend(ax3, legtxt, 'Location','best');

ax4 = nexttile(tlo); hold(ax4,'on');
theta = linspace(0, 2*pi, 400);
plot(ax4, cos(theta), sin(theta), 'k--');
scatter(ax4, real(lambda(~isCand)), imag(lambda(~isCand)), 35, [0.3 0.3 0.3], 'filled');
scatter(ax4, real(lambda(isCand)), imag(lambda(isCand)), 55, 'r', 'filled');
axis(ax4,'equal');
xlabel(ax4,'Real(\lambda)');
ylabel(ax4,'Imag(\lambda)');
title(ax4, sprintf('LDS eigenvalues | candidate pairs = %d', size(pairs,1)));
grid(ax4,'on');

ax5 = nexttile(tlo);
imagesc(ax5, PLV, [0 1]);
axis(ax5,'image');
xlabel(ax5,'PC');
ylabel(ax5,'PC');
title(ax5,'Phase-locking value (PLV)');
colorbar(ax5);

ax6 = nexttile(tlo, [1 2]);
axis(ax6,'off');
txt = {
    sprintf('Kept neurons: %d', OUT.summary.nKeptNeurons)
    sprintf('Latent PCs: %d', OUT.summary.numPCs)
    sprintf('PC dom freqs (Hz): %s', mat2str(OUT.summary.pcDomFreqHz, 4))
    sprintf('PC freq CV: %.3f', OUT.summary.pcFreqCV)
    sprintf('Mean PLV: %.3f', OUT.summary.meanPLV)
    sprintf('Phase components: %d', OUT.summary.nPhaseComponents)
    sprintf('LDS oscillator pairs: %d', OUT.summary.nOscillatorPairsLDS)
    sprintf('Mean sliding pairs: %.3f', OUT.summary.meanSlidingPairs)
    sprintf('Max sliding pairs: %.3f', OUT.summary.maxSlidingPairs)
    sprintf('Mode counts: %s', mat2str(OUT.summary.modeCount(:)'))
    sprintf('LDS pair freqs (Hz): %s', mat2str(OUT.summary.ldsPairFreqHz, 4))
    ['Interpretation: ' char(OUT.summary.interpretation)]
    };
text(ax6, 0.01, 0.98, txt, 'VerticalAlignment','top', 'FontName','Courier');

% ===================== PHASE FIGURE =====================
figs.phase = figure('Color','w','Position',[100 100 1450 850], 'Name','Phase coupling details');
tlo2 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

axA = nexttile(tlo2); hold(axA,'on');
for i = 1:numPCs
    plot(axA, tSec, OUT.phase.latentBand(:,i), 'LineWidth', 1);
end
if ~isempty(eventTimesSec)
    yl = ylim(axA);
    for i = 1:numel(eventTimesSec)
        plot(axA, [eventTimesSec(i) eventTimesSec(i)], yl, 'r:');
    end
end
xlabel(axA,'Time (s)');
ylabel(axA,'Bandpassed latent signal');
title(axA,'Bandpassed PCs');
legend(axA, compose('PC%d', 1:numPCs), 'Location','best');
axis(axA,'tight');

axB = nexttile(tlo2); hold(axB,'on');
for i = 1:min(numPCs,4)
    ph = unwrap(OUT.phase.latentPhase(:,i));
    good = OUT.phase.phaseMask(:,i);
    ph(~good) = NaN;
    plot(axB, tSec, ph, 'LineWidth', 1);
end
if ~isempty(eventTimesSec)
    yl = ylim(axB);
    for i = 1:numel(eventTimesSec)
        plot(axB, [eventTimesSec(i) eventTimesSec(i)], yl, 'r:');
    end
end
xlabel(axB,'Time (s)');
ylabel(axB,'Unwrapped phase (rad)');
title(axB,'Unwrapped phases (masked low amp)');
legend(axB, compose('PC%d', 1:min(numPCs,4)), 'Location','best');
axis(axB,'tight');

axC = nexttile(tlo2);
imagesc(axC, tSec, 1:numPCs, OUT.phase.latentAmp');
xlabel(axC,'Time (s)');
ylabel(axC,'PC');
title(axC,'Instantaneous amplitude');
colorbar(axC);
axis(axC,'tight');

axD = nexttile(tlo2);
imagesc(axD, OUT.phase.PLV, [0 1]);
axis(axD,'image');
xlabel(axD,'PC');
ylabel(axD,'PC');
title(axD,'PLV matrix');
colorbar(axD);

% ===================== SLIDING FIGURE =====================
figs.sliding = figure('Color','w','Position',[110 110 1500 900], 'Name','Sliding-window oscillator summary');
tlo3 = tiledlayout(3,2, 'TileSpacing','compact', 'Padding','compact');

% first row: manifold spanning both columns
axS3 = nexttile(tlo3, [1 2]); hold(axS3,'on');
if size(latentTS,2) >= 3
    plot3(axS3, latentTS(:,1), latentTS(:,2), latentTS(:,3), '-', 'Color', [0.8 0.8 0.8]);
    cFrames = round(OUT.sliding.centerSec * fps) + 1;
    cFrames = cFrames(cFrames >= 1 & cFrames <= size(latentTS,1));
    scatter3(axS3, latentTS(cFrames,1), latentTS(cFrames,2), latentTS(cFrames,3), ...
        45, OUT.sliding.nOscPairs(1:numel(cFrames)), 'filled');
    xlabel(axS3,'PC1'); ylabel(axS3,'PC2'); zlabel(axS3,'PC3');
    title(axS3,'Trajectory colored by sliding-window oscillator count');
    colorbar(axS3);
    view(axS3,45,30); grid(axS3,'on'); box(axS3,'off');
end

axS1 = nexttile(tlo3); hold(axS1,'on');
plot(axS1, OUT.sliding.centerSec/60, OUT.sliding.nOscPairs, '-o', 'LineWidth', 1);
xlabel(axS1,'Window center (min)');
ylabel(axS1,'# oscillator pairs');
title(axS1,'Sliding-window LDS oscillator count');

axS2 = nexttile(tlo3); hold(axS2,'on');
for i = 1:size(OUT.sliding.pcDomFreqHz,2)
    plot(axS2, OUT.sliding.centerSec/60, OUT.sliding.pcDomFreqHz(:,i), '-', 'LineWidth', 1);
end
xlabel(axS2,'Window center (min)');
ylabel(axS2,'Dominant frequency (Hz)');
title(axS2,'Sliding-window latent PC dominant frequencies');

axS4 = nexttile(tlo3); hold(axS4,'on');
bar(axS4, OUT.sliding.centerSec/60, OUT.sliding.nEvents);
xlabel(axS4,'Window center (min)');
ylabel(axS4,'# stim events');
title(axS4,'Stim events per window');

axS5 = nexttile(tlo3);
imagesc(axS5, OUT.sliding.centerSec/60, 1:size(OUT.sliding.pcDomFreqHz,2), OUT.sliding.pcDomFreqHz');
xlabel(axS5,'Window center (min)');
ylabel(axS5,'PC');
title(axS5,'Sliding PC dominant frequencies');
colorbar(axS5);
end

%% =======================================================================
function figs = make_mode_plots(OUT, fps)
figs = struct();

Xplot          = OUT.data.Xz;
sortOrd        = OUT.mode.sortOrd;
modeID_sorted  = OUT.mode.modeID_sorted;
modeID         = OUT.mode.modeID;
modeMean       = OUT.mode.modeMean;
domMode_t      = OUT.mode.domMode_t;
latentTS       = OUT.pca.latentTS;
numPCs         = OUT.summary.numPCs;
eventTimesSec  = OUT.data.eventTimesSec;
positions3D    = OUT.data.positions3D;

Te = size(Xplot,2);
tSec = (0:Te-1) / fps;
Xplot_sorted = Xplot(sortOrd, :);
modeColors = get_mode_colors(numPCs);

% force row vector for time strip
domMode_t = domMode_t(:)';
if numel(domMode_t) ~= Te
    error('domMode_t must have one entry per time point.');
end

% ===================== FIGURE 1: MODE MEMBERSHIP =====================
figs.modeMembership = figure( ...
    'Color','w', ...
    'Position', [100 100 1250 760], ...
    'Name', 'Mode membership');

% ---------- layout geometry ----------
leftMain   = 0.08;
bottomMain = 0.10;
widthMain  = 0.68;
heightMain = 0.74;

gapTop     = 0.015;
heightTop  = 0.045;

gapRight   = 0.03;
widthRight = 0.03;

% ---------- main heatmap ----------
ax1 = axes('Position', [leftMain bottomMain widthMain heightMain]);
imagesc(ax1, Xplot_sorted, [0 2]);
xlabel(ax1, 'time (frames)');
ylabel(ax1, 'neurons');
title(ax1, 'Neurons sorted by dominant latent mode');
colormap(ax1, flipud(bone));
axis(ax1, 'tight');

hold(ax1, 'on');
cuts = find(diff(modeID_sorted) ~= 0);
for i = 1:numel(cuts)
    yline(ax1, cuts(i)+0.5, 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
end

% ---------- top strip: dominant mode over time ----------
axTop = axes('Position', [leftMain, bottomMain + heightMain + gapTop, widthMain, heightTop]);
imagesc(axTop, domMode_t);   % 1 x T categorical strip
set(axTop, 'YTick', []);
set(axTop, 'XTick', []);
title(axTop, 'dominant mode over time');
colormap(axTop, modeColors);
caxis(axTop, [1 numPCs]);
axis(axTop, 'tight');

% ---------- right strip: neuron dominant mode ----------
ax2 = axes('Position', [leftMain + widthMain + gapRight, bottomMain, widthRight, heightMain]);
imagesc(ax2, modeID_sorted(:));
set(ax2, 'XTick', []);
ylabel(ax2, 'neurons');
title(ax2, 'mode');
axis(ax2, 'tight');
colormap(ax2, modeColors);
caxis(ax2, [1 numPCs]);

cb = colorbar(ax2, 'eastoutside');
cb.Label.String = 'mode id';
cb.Ticks = 1:numPCs;

% ---------- keep top strip aligned with raster ----------
linkaxes([ax1 axTop], 'x');

% Match x-limits explicitly
xlim(ax1, [1 Te]);
xlim(axTop, [1 Te]);

end

% ===================== FIGURE 2: CLEAN MODE DYNAMICS =====================
figs.modeDynamics = figure('Color','w', 'Position', [120 80 1600 920], ...
    'Name', 'Mode dynamics');

% manual layout for robust control
leftX = 0.06;
leftW = 0.42;
rightX = 0.54;
rightW = 0.40;
bottom0 = 0.08;
top0 = 0.93;
gap = 0.012;
rowH = (top0 - bottom0 - gap*(numPCs-1)) / numPCs;

traceAxes = gobjects(numPCs,1);

for k = 1:numPCs
    y = top0 - k*rowH - (k-1)*gap;
    ax = axes('Position', [leftX y leftW rowH]);
    traceAxes(k) = ax;
    hold(ax, 'on');

    plot(ax, tSec, modeMean(k,:), 'Color', modeColors(k,:), 'LineWidth', 1.15);

    if ~isempty(eventTimesSec)
        yl = ylim(ax);
        for ie = 1:numel(eventTimesSec)
            plot(ax, [eventTimesSec(ie) eventTimesSec(ie)], yl, ':', ...
                'Color', [0.85 0.45 0.45], 'LineWidth', 0.8);
            yl = ylim(ax); % refresh after plot
        end
    end

    ylabel(ax, sprintf('mode %d', k), 'Color', modeColors(k,:));
    ax.YColor = modeColors(k,:);
    xlim(ax, [tSec(1) tSec(end)]);
    grid(ax, 'on');
    box(ax, 'off');

    if k < numPCs
        set(ax, 'XTickLabel', []);
    else
        xlabel(ax, 'time (s)');
    end
end

linkaxes(traceAxes, 'x');

annotation(figs.modeDynamics, 'textbox', [leftX 0.94 leftW 0.04], ...
    'String', 'Mode activity over time', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);

% ===================== RIGHT: one manifold plot with colored segments =====================
axR = axes('Position', [rightX 0.12 rightW 0.78]);
hold(axR, 'on');

% optional faint background trajectory
if size(latentTS,2) >= 3
    plot3(axR, latentTS(:,1), latentTS(:,2), latentTS(:,3), ...
        '-', 'Color', [0.90 0.90 0.90], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
else
    plot(axR, latentTS(:,1), latentTS(:,2), ...
        '-', 'Color', [0.90 0.90 0.90], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

% draw colored segments by dominant mode
if size(latentTS,2) >= 3
    for k = 1:numPCs
        idx = find(domMode_t(:) == k);
        if isempty(idx)
            continue;
        end

        % split into contiguous stretches
        breaks = [1; find(diff(idx) > 1) + 1; numel(idx)+1];

        for b = 1:numel(breaks)-1
            segIdx = idx(breaks(b):breaks(b+1)-1);

            if numel(segIdx) >= 2
                plot3(axR, latentTS(segIdx,1), latentTS(segIdx,2), latentTS(segIdx,3), ...
                    '-', 'Color', modeColors(k,:), 'LineWidth', 1.8, ...
                    'HandleVisibility', 'off');
            else
                scatter3(axR, latentTS(segIdx,1), latentTS(segIdx,2), latentTS(segIdx,3), ...
                    10, modeColors(k,:), 'filled', ...
                    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.9, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    xlabel(axR, 'PC1');
    ylabel(axR, 'PC2');
    zlabel(axR, 'PC3');
    view(axR, 40, 25);

else
    for k = 1:numPCs
        idx = find(domMode_t(:) == k);
        if isempty(idx)
            continue;
        end

        breaks = [1; find(diff(idx) > 1) + 1; numel(idx)+1];

        for b = 1:numel(breaks)-1
            segIdx = idx(breaks(b):breaks(b+1)-1);

            if numel(segIdx) >= 2
                plot(axR, latentTS(segIdx,1), latentTS(segIdx,2), ...
                    '-', 'Color', modeColors(k,:), 'LineWidth', 1.8, ...
                    'HandleVisibility', 'off');
            else
                scatter(axR, latentTS(segIdx,1), latentTS(segIdx,2), ...
                    10, modeColors(k,:), 'filled', ...
                    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.9, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    xlabel(axR, 'PC1');
    ylabel(axR, 'PC2');
end

title(axR, 'Latent trajectory colored by dominant mode');
grid(axR, 'on');
box(axR, 'off');

% --------- create clean dummy handles for legend ----------
legH = gobjects(numPCs,1);
for k = 1:numPCs
    if size(latentTS,2) >= 3
        legH(k) = plot3(axR, nan, nan, nan, ...
            '-', 'Color', modeColors(k,:), 'LineWidth', 2.2);
    else
        legH(k) = plot(axR, nan, nan, ...
            '-', 'Color', modeColors(k,:), 'LineWidth', 2.2);
    end
end

legend(axR, legH, compose('Mode %d', 1:numPCs), 'Location', 'eastoutside');
% ===================== FIGURE 3: BRAIN MODE SCATTER =====================
if ~isempty(positions3D)
    figs.modeBrain = figure('Color','w', 'Position', [140 70 1800 900], ...
        'Name', 'Brain-space mode membership');

    nCols = min(5, numPCs);
    nRows = ceil(numPCs / nCols);
    tloB = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:numPCs
        axB = nexttile(tloB);
        hold(axB, 'on');

        idx = (modeID == k);
        idxOther = ~idx;

        scatter3(axB, positions3D(idxOther,1), positions3D(idxOther,2), positions3D(idxOther,3), ...
            8, [0.85 0.85 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.20, 'MarkerEdgeAlpha', 0.20);

        scatter3(axB, positions3D(idx,1), positions3D(idx,2), positions3D(idx,3), ...
            18, modeColors(k,:), 'filled', ...
            'MarkerFaceAlpha', 0.95, 'MarkerEdgeAlpha', 0.95);

        title(axB, sprintf('Mode %d', k), ...
            'Color', modeColors(k,:), 'FontWeight', 'bold');

        xlabel(axB, 'X');
        ylabel(axB, 'Y');
        zlabel(axB, 'Z');
        axis(axB, 'equal');
        grid(axB, 'on');
        box(axB, 'off');
        view(axB, 3);
    end
end


%% =======================================================================
function modeColors = get_mode_colors(numModes)
% Consistent categorical colors for modes

custom = [
    0.75 0.00 0.00   % dark red
    0.85 0.30 0.00   % orange-red
    0.95 0.55 0.10   % orange
    0.95 0.75 0.20   % yellow-orange
    0.70 0.85 0.20   % yellow-green
    0.25 0.80 0.65   % green-cyan
    0.20 0.65 0.90   % blue
    0.35 0.45 0.90   % indigo
    0.50 0.35 0.80   % purple
    0.35 0.10 0.35   % dark purple
    0.50 0.50 0.50   % gray
    0.75 0.45 0.75   % mauve
    0.45 0.70 0.25   % olive
    0.10 0.60 0.55   % teal
    0.60 0.40 0.20   % brown
    ];

if numModes <= size(custom,1)
    modeColors = custom(1:numModes, :);
else
    base = lines(numModes);
    modeColors = base;
end
end