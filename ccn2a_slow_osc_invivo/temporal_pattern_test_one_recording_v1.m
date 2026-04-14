function OUT = temporal_pattern_test_one_recording_v1(dff_full, fps, varargin)
% TEMPORAL_PATTERN_TEST_ONE_RECORDING_V1
%
% Transpose-PCA version:
%   - neurons are observations
%   - timepoints are variables
%   - PCs are temporal patterns
%
% Interpretation
% --------------
%   tpca.temporalPatterns : [time x K]
%       Temporal motifs extracted by PCA.
%
%   tpca.neuronScores     : [neurons x K]
%       Loading / weight of each neuron onto each temporal motif.
%
%   mode.modeID           : [neurons x 1]
%       Dominant temporal-pattern assignment per neuron.
%
% Inputs
% ------
%   dff_full : [neurons x time]
%   fps      : scalar
%
% Name-value parameters
% ---------------------
%   'EpochFrames'            : [start end] frames, default = full recording
%   'EventFrames'            : event frames in absolute frame index
%   'FiniteFracThresh'       : minimum fraction of finite samples per neuron
%   'UseZScore'              : z-score neurons before PCA
%   'NumPCs'                 : number of temporal PCs to use
%   'BandHz'                 : band for phase analyses
%   'WelchWinSec'            : PSD Welch window length in sec
%   'WelchOverlapFrac'       : PSD overlap fraction
%   'LDSReg'                 : ridge regularization for LDS fit
%   'MinEigMag'              : minimum eigenvalue magnitude for oscillator candidate
%   'MaxEigMag'              : maximum eigenvalue magnitude for oscillator candidate
%   'PhaseAmpQuantileMin'    : mask lowest-amplitude fraction for phase usage
%   'PLVThreshold'           : threshold for phase-lock graph
%   'WinSec'                 : sliding window length in sec
%   'StepSec'                : sliding window step in sec
%
% ACF filter parameters
% ---------------------
%   'EnableAcfPeakFilter'    : true/false
%   'AcfPeakLagRangeSec'     : lag range where peaks are counted
%   'AcfMinProminence'       : MinPeakProminence for ACF peaks
%   'MinNumAcPeaks'          : minimum number of ACF peaks to keep neuron
%   'AcfMaxLagSec'           : maximum lag for ACF calculation
%
% Plotting
% --------
%   'MakeModeMembershipPlot' : true/false
%   'NeuronPositions'        : [neurons x 2] or [neurons x 3]
%   'DoPlots'                : true/false
%   'Verbose'                : true/false
%
% Output
% ------
%   OUT.params
%   OUT.data
%   OUT.filter
%   OUT.tpca
%   OUT.mode
%   OUT.psd
%   OUT.phase
%   OUT.lds
%   OUT.sliding
%   OUT.summary
%   OUT.fig
%
% -------------------------------------------------------------------------

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

addParameter(ip, 'EnableAcfPeakFilter', false, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'AcfPeakLagRangeSec', [10 120], @(x) isnumeric(x) && numel(x)==2 && x(1)>=0 && x(2)>x(1));
addParameter(ip, 'AcfMinProminence', 0.08, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(ip, 'MinNumAcPeaks', 2, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(ip, 'AcfMaxLagSec', 180, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'MakeModeMembershipPlot', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'NeuronPositions', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x) && (size(x,2)==2 || size(x,2)==3)));

addParameter(ip, 'DoPlots', false, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'Verbose', false, @(x) islogical(x) || isnumeric(x));

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
    if epochFrames(2) <= epochFrames(1)
        error('Resolved epoch frame interval is invalid.');
    end
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

Xraw = dff_epoch(keepIdx, :);   % neurons x time

if P.UseZScore
    X = zscore_fill_rows(Xraw);
else
    X = Xraw;
    X(~isfinite(X)) = 0;
end

%% ======================= OPTIONAL NEURON POSITIONS =====================
neuronPositionsRaw = P.NeuronPositions;
neuronPositionsKept = [];

if ~isempty(neuronPositionsRaw)
    if size(neuronPositionsRaw,1) ~= N0
        warning('NeuronPositions row count (%d) does not match input neurons (%d). Ignoring positions.', ...
            size(neuronPositionsRaw,1), N0);
        neuronPositionsRaw = [];
    else
        neuronPositionsKept = neuronPositionsRaw(keepIdx,:);
        if size(neuronPositionsKept,2) == 2
            neuronPositionsKept(:,3) = 0;
        end
    end
end

%% ============================== TIME-PCA ===============================
% rows = neurons (observations)
% cols = timepoints (variables)
%
% coeff  -> temporal patterns [time x K]
% score  -> neuron scores     [neurons x K]

maxPC = min([P.NumPCs, rank(X), size(X,1), size(X,2)]);
[coeff, score, latent, ~, explained, mu] = pca(X, 'NumComponents', maxPC);

numPCs = maxPC;

% Temporal patterns with amplitude scaled by sqrt(latent)
temporalPatterns = coeff(:,1:numPCs) .* sqrt(latent(1:numPCs))';
neuronScores     = score(:,1:numPCs);

%% ======================= NEURON MODE ASSIGNMENT ========================
[modeLoadingAbs, modeID] = max(abs(neuronScores), [], 2);
modeSign = sign(neuronScores(sub2ind(size(neuronScores), (1:size(neuronScores,1))', modeID)));

sortTable = [modeID(:), -modeLoadingAbs(:)];
[~, sortOrd] = sortrows(sortTable, [1 2]);
modeID_sorted = modeID(sortOrd);

modeCount = zeros(numPCs,1);
for k = 1:numPCs
    modeCount(k) = nnz(modeID == k);
end

%% ====================== TEMPORAL-PATTERN PSD ===========================
[pcDomFreq, pcPeakPower, fAxis, PxxPC] = compute_pc_psd( ...
    temporalPatterns, fps, P.WelchWinSec, P.WelchOverlapFrac, P.BandHz);

%% ======================= BANDPASS + HILBERT ============================
[latentBand, latentPhase, latentAmp, PLV, meanPLV] = ...
    compute_phase_metrics(temporalPatterns, fps, P.BandHz, P.PhaseAmpQuantileMin);

adj = false(numPCs, numPCs);
goodPLV = isfinite(PLV);
adj(goodPLV) = PLV(goodPLV) >= P.PLVThreshold;
adj(1:numPCs+1:end) = false;
nPhaseComponents = count_connected_components(adj);

%% =============================== LDS ===================================
[A, lambda, pairs, pairFreqHz] = fit_latent_lds( ...
    temporalPatterns, fps, P.BandHz, P.LDSReg, P.MinEigMag, P.MaxEigMag);

pcFreqCV = std(pcDomFreq, 'omitnan') / mean(pcDomFreq, 'omitnan');

%% ============================= SLIDING =================================
sliding = compute_sliding_latent_metrics( ...
    temporalPatterns, fps, P.WinSec, P.StepSec, ...
    P.BandHz, P.WelchWinSec, P.WelchOverlapFrac, P.PhaseAmpQuantileMin, ...
    P.LDSReg, P.MinEigMag, P.MaxEigMag, P.PLVThreshold);

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
summary.meanSlidingPairs    = mean(sliding.nPairsLDS, 'omitnan');
summary.maxSlidingPairs     = max(sliding.nPairsLDS, [], 'omitnan');
summary.modeCount           = modeCount;

if size(pairs,1) <= 1 && isfinite(meanPLV) && meanPLV > 0.8
    summary.interpretation = "One dominant temporal oscillator + noise.";
elseif size(pairs,1) >= 2 && (pcFreqCV > 0.15 || (~isnan(meanPLV) && meanPLV < 0.7) || nPhaseComponents >= 2)
    summary.interpretation = "Multiple temporal oscillatory motifs.";
elseif summary.maxSlidingPairs >= 2
    summary.interpretation = "Mixed global fit; transient temporal oscillator motifs.";
else
    summary.interpretation = "Weakly coupled / overlapping temporal motifs.";
end

if P.Verbose
    fprintf('\n=== temporal_pattern_test_one_recording_v1 ===\n');
    fprintf('Input neurons               : %d\n', N0);
    fprintf('Epoch frames                : %d to %d (%d frames)\n', epochFrames(1), epochFrames(2), Te);
    fprintf('Kept neurons                : %d\n', numel(keepIdx));
    fprintf('Temporal PCs                : %d\n', numPCs);
    fprintf('PC dominant freqs (Hz)      : %s\n', mat2str(pcDomFreq, 4));
    fprintf('PC frequency CV             : %.4f\n', pcFreqCV);
    fprintf('Mean off-diagonal PLV       : %.4f\n', meanPLV);
    fprintf('Phase components            : %d\n', nPhaseComponents);
    fprintf('LDS oscillatory pairs       : %d\n', size(pairs,1));
    fprintf('Mean sliding pairs          : %.3f\n', summary.meanSlidingPairs);
    fprintf('Max sliding pairs           : %.3f\n', summary.maxSlidingPairs);
    fprintf('Interpretation              : %s\n', summary.interpretation);
end

%% ============================== OUTPUT =================================
OUT = struct();
OUT.params = P;

OUT.data = struct();
OUT.data.keepIdx             = keepIdx;
OUT.data.dropIdx             = dropIdx;
OUT.data.epochFrames         = epochFrames;
OUT.data.eventFramesAbs      = eventFrames;
OUT.data.eventFramesLocal    = eventFramesLocal;
OUT.data.eventTimesSec       = eventTimesSec;
OUT.data.Xraw                = Xraw;
OUT.data.Xz                  = X;
OUT.data.neuronPositionsRaw  = neuronPositionsRaw;
OUT.data.neuronPositionsKept = neuronPositionsKept;

OUT.filter = struct();
OUT.filter.keepBase      = keepBase;
OUT.filter.keepAcf       = keepAcf;
OUT.filter.keepFinal     = keep;
OUT.filter.nAcPeaks_all  = nAcPeaks_all;

OUT.tpca = struct();
OUT.tpca.coeff            = coeff;
OUT.tpca.score            = score;
OUT.tpca.temporalPatterns = temporalPatterns;
OUT.tpca.neuronScores     = neuronScores;
OUT.tpca.latent           = latent;
OUT.tpca.explained        = explained;
OUT.tpca.mu               = mu;

OUT.mode = struct();
OUT.mode.modeID          = modeID;
OUT.mode.modeLoadingAbs  = modeLoadingAbs;
OUT.mode.modeSign        = modeSign;
OUT.mode.sortOrd         = sortOrd;
OUT.mode.modeID_sorted   = modeID_sorted;
OUT.mode.modeCount       = modeCount;

OUT.psd = struct();
OUT.psd.fAxis        = fAxis;
OUT.psd.PxxPC        = PxxPC;
OUT.psd.pcDomFreqHz  = pcDomFreq;
OUT.psd.pcPeakPower  = pcPeakPower;

OUT.phase = struct();
OUT.phase.latentBand  = latentBand;
OUT.phase.latentPhase = latentPhase;
OUT.phase.latentAmp   = latentAmp;
OUT.phase.PLV         = PLV;

OUT.lds = struct();
OUT.lds.A           = A;
OUT.lds.lambda      = lambda;
OUT.lds.pairs       = pairs;
OUT.lds.pairFreqHz  = pairFreqHz;

OUT.sliding = sliding;
OUT.summary = summary;
OUT.fig = struct();

if P.DoPlots
    OUT.fig = make_time_summary_plots(OUT, fps);

    if P.MakeModeMembershipPlot
        modeFigs = make_time_mode_plots(OUT, fps);
        close;
        OUT.fig.modeMembership = modeFigs.modeMembership;
        OUT.fig.modeDynamics   = modeFigs.modeDynamics;
        OUT.fig.neuronMap      = modeFigs.neuronMap;
    end
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
function [pcDomFreq, pcPeakPower, fAxis, PxxPC] = compute_pc_psd(XtimeByPC, fps, winSec, overlapFrac, bandHz)
nPC = size(XtimeByPC,2);

nwin = max(16, round(winSec * fps));
nwin = min(nwin, size(XtimeByPC,1));
nover = round(nwin * overlapFrac);
nfft = max(256, 2^nextpow2(nwin));

pcDomFreq  = nan(nPC,1);
pcPeakPower = nan(nPC,1);
PxxPC = [];
fAxis = [];

for i = 1:nPC
    xi = XtimeByPC(:,i);
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
function [latentBand, latentPhase, latentAmp, PLV, meanPLV] = compute_phase_metrics(XtimeByPC, fps, bandHz, phaseAmpQuantileMin)

Wn = bandHz / (fps/2);
if Wn(2) >= 1
    error('Upper band edge must be below Nyquist.');
end

[b,a] = butter(3, Wn, 'bandpass');

latentBand = nan(size(XtimeByPC));
for i = 1:size(XtimeByPC,2)
    xi = XtimeByPC(:,i);
    xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
    latentBand(:,i) = filtfilt(b, a, xi);
end

analytic = hilbert(latentBand);
latentPhase = angle(analytic);
latentAmp   = abs(analytic);

phaseMask = true(size(latentAmp));
for i = 1:size(latentAmp,2)
    ai = latentAmp(:,i);
    goodAi = ai(isfinite(ai));
    if isempty(goodAi)
        phaseMask(:,i) = false(size(ai));
        continue;
    end
    thr = quantile(goodAi, phaseAmpQuantileMin);
    phaseMask(:,i) = ai >= thr;
end

nPC = size(XtimeByPC,2);
PLV = nan(nPC,nPC);

for i = 1:nPC
    for j = 1:nPC
        good = isfinite(latentPhase(:,i)) & isfinite(latentPhase(:,j)) & ...
               phaseMask(:,i) & phaseMask(:,j);
        if nnz(good) >= 20
            dphi = latentPhase(good,i) - latentPhase(good,j);
            PLV(i,j) = abs(mean(exp(1i*dphi)));
        end
    end
end

maskOff = triu(true(nPC),1);
vals = PLV(maskOff);
vals = vals(isfinite(vals));
if isempty(vals)
    meanPLV = NaN;
else
    meanPLV = mean(vals);
end
end

%% =======================================================================
function [A, lambda, pairs, pairFreqHz] = fit_latent_lds(XtimeByPC, fps, bandHz, ldsReg, minEigMag, maxEigMag)

X1 = XtimeByPC(1:end-1,:)';
X2 = XtimeByPC(2:end,:)';

A = (X2 * X1') / (X1 * X1' + ldsReg * eye(size(X1,1)));

[~,D] = eig(A);
lambda = diag(D);

eigMag  = abs(lambda);
eigAng  = angle(lambda);
eigFreq = abs(eigAng) * fps / (2*pi);

isComplex = abs(imag(lambda)) > 1e-10;
isCandidate = isComplex & eigMag >= minEigMag & eigMag <= maxEigMag & ...
    eigFreq >= bandHz(1) & eigFreq <= bandHz(2);

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
        used(ii) = true;
        used(jj) = true;
        pairs = [pairs; oscIdx(ii) oscIdx(jj)]; %#ok<AGROW>
    else
        used(ii) = true;
    end
end

pairFreqHz = nan(size(pairs,1),1);
for i = 1:size(pairs,1)
    thisFreq = eigFreq(pairs(i,:));
    pairFreqHz(i) = mean(thisFreq,'omitnan');
end
end

%% =======================================================================
function sliding = compute_sliding_latent_metrics(XtimeByPC, fps, slideWinSec, slideStepSec, ...
    bandHz, welchWinSec, welchOverlapFrac, phaseAmpQuantileMin, ...
    ldsReg, minEigMag, maxEigMag, plvThreshold)

T = size(XtimeByPC,1);
winF  = max(20, round(slideWinSec*fps));
stepF = max(1, round(slideStepSec*fps));

starts = 1:stepF:max(1, T-winF+1);
nW = numel(starts);

sliding = struct();
sliding.startFrames = starts(:);
sliding.centerFrames = starts(:) + floor(winF/2);
sliding.centerSec = (sliding.centerFrames - 1) / fps;
sliding.meanPLV   = nan(nW,1);
sliding.pcFreqCV  = nan(nW,1);
sliding.nPairsLDS = nan(nW,1);

for w = 1:nW
    s = starts(w);
    e = min(T, s + winF - 1);
    Xi = XtimeByPC(s:e,:);

    [pcDomFreq, ~, ~, ~] = compute_pc_psd(Xi, fps, min(welchWinSec, (e-s+1)/fps), welchOverlapFrac, bandHz);
    [~,~,~,~,meanPLV] = compute_phase_metrics(Xi, fps, bandHz, phaseAmpQuantileMin);
    [~,~,pairs,~] = fit_latent_lds(Xi, fps, bandHz, ldsReg, minEigMag, maxEigMag);

    sliding.meanPLV(w)  = meanPLV;
    sliding.pcFreqCV(w) = std(pcDomFreq,'omitnan') / mean(pcDomFreq,'omitnan');
    sliding.nPairsLDS(w)= size(pairs,1);
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
        stack = [stack nbrs]; %#ok<AGROW>
    end
end
end

%% =======================================================================
function figs = make_time_summary_plots(OUT, fps)
% MAKE_TIME_SUMMARY_PLOTS
% One-stop plotting function for time-PCA / temporal motif decomposition.
%
% Main figure:
%   - cumulative explained variance
%   - temporal motif traces
%   - PSDs of temporal motifs
%   - neuron loading matrix
%   - dominant temporal motif in anatomy
%   - dominant motif counts
%
% Supplementary figures:
%   - phase diagnostics
%   - sliding diagnostics
%   - LDS eigenvalue diagnostic
%
% -----------------------------------------------------------------------

figs = struct();

% ============================== UNPACK ==================================
temporalPatterns = OUT.tpca.temporalPatterns;   % [time x K]
neuronScores     = OUT.tpca.neuronScores;       % [neurons x K]
explained        = OUT.tpca.explained(:);       % [%]
fAxis            = OUT.psd.fAxis;
PxxPC            = OUT.psd.PxxPC;
pcDomFreq        = OUT.psd.pcDomFreqHz(:);
PLV              = OUT.phase.PLV;
lambda           = OUT.lds.lambda;
pairs            = OUT.lds.pairs;
eventTimesSec    = OUT.data.eventTimesSec(:);
modeID           = OUT.mode.modeID(:);
modeCount        = OUT.summary.modeCount(:);
pos              = OUT.data.neuronPositionsKept;

numPCs = size(temporalPatterns,2);
tSec   = (0:size(temporalPatterns,1)-1) / fps;

modeColors = get_mode_colors(numPCs);

% ====================== SORT NEURONS FOR DISPLAY ========================
[domAbs, domID] = max(abs(neuronScores), [], 2);
sortTable = [domID(:), -domAbs(:)];
[~, sortOrd] = sortrows(sortTable, [1 2]);

scoresSorted = neuronScores(sortOrd, :);
domID_sorted = domID(sortOrd);

% z-score columns only for heatmap display
scoresDisp = scoresSorted;
for k = 1:size(scoresDisp,2)
    x = scoresDisp(:,k);
    mu = mean(x,'omitnan');
    sd = std(x,0,'omitnan');
    if ~isfinite(sd) || sd == 0
        sd = 1;
    end
    scoresDisp(:,k) = (x - mu) ./ sd;
end
scoresDisp(~isfinite(scoresDisp)) = 0;

% ====================== EXPLAINED VARIANCE CURVES =======================
nExpPlot = min(numPCs, 10);
expVals  = explained(1:nExpPlot);
cumExp   = cumsum(expVals);

% ====================== WHICH MOTIFS TO PLOT ============================
nMotifPlot = min(numPCs, 6);

% ====================== SUMMARY FIGURE ==================================
figs.summary = figure( ...
    'Color','w', ...
    'Position',[50 35 1850 1050], ...
    'Name','Temporal motif decomposition');

tlo = tiledlayout(3,4, 'TileSpacing','compact', 'Padding','compact');

% -----------------------------------------------------------------------
% A) cumulative explained variance
% -----------------------------------------------------------------------
ax1 = nexttile(tlo, 1); hold(ax1,'on');

bar(ax1, 1:nExpPlot, expVals, ...
    'FaceColor',[0.72 0.82 0.93], ...
    'EdgeColor','none', ...
    'BarWidth',0.75);

plot(ax1, 1:nExpPlot, cumExp, '-o', ...
    'Color',[0.10 0.25 0.55], ...
    'LineWidth',2.2, ...
    'MarkerSize',5, ...
    'MarkerFaceColor',[0.10 0.25 0.55]);

xlabel(ax1,'Temporal PC');
ylabel(ax1,'% variance');
title(ax1,'Explained variance');
legend(ax1, {'Individual variance','Cumulative variance'}, ...
    'Location','northwest', 'Box','off');
grid(ax1,'on');
box(ax1,'off');

% -----------------------------------------------------------------------
% B) temporal motif traces
% -----------------------------------------------------------------------
ax2 = nexttile(tlo, [1 2]); hold(ax2,'on');

spacing = 4.5;
tickPos = spacing * (0:(nMotifPlot-1));

for k = 1:nMotifPlot
    x = temporalPatterns(:,k);
    s = std(x,0,'omitnan');
    if ~isfinite(s) || s == 0
        s = 1;
    end
    xdisp = x ./ s;

    y0 = tickPos(nMotifPlot - k + 1);   % motif 1 at top

    plot(ax2, tSec, xdisp + y0, ...
        'Color', modeColors(k,:), ...
        'LineWidth', 1.15);
end

if ~isempty(eventTimesSec)
    yl = ylim(ax2);
    for i = 1:numel(eventTimesSec)
        plot(ax2, [eventTimesSec(i) eventTimesSec(i)], yl, ':', ...
            'Color', [0.82 0.82 0.82], 'LineWidth', 0.6);
    end
end

yticks(ax2, tickPos);
yticklabels(ax2, compose('Motif %d', nMotifPlot:-1:1));
ylim(ax2, [-spacing*0.5, tickPos(end) + spacing*1.2]);

xlabel(ax2,'Time (s)');
ylabel(ax2,'Temporal motifs');
title(ax2,'Temporal motif traces');
grid(ax2,'on');
box(ax2,'off');
axis(ax2,'tight');

% -----------------------------------------------------------------------
% C) PSDs of temporal motifs
% -----------------------------------------------------------------------
ax3 = nexttile(tlo, [1 1]); hold(ax3,'on');

for k = 1:nMotifPlot
    % main line (this goes to legend)
    plot(ax3, fAxis, PxxPC(k,:), ...
        'Color', modeColors(k,:), ...
        'LineWidth', 1.5);

    % peak marker (excluded from legend)
    if isfinite(pcDomFreq(k))
        [~, idxPk] = min(abs(fAxis - pcDomFreq(k)));
        if ~isempty(idxPk) && idxPk >= 1 && idxPk <= numel(fAxis)
            plot(ax3, fAxis(idxPk), PxxPC(k,idxPk), 'o', ...
                'Color', modeColors(k,:), ...
                'MarkerFaceColor', modeColors(k,:), ...
                'MarkerSize', 4, ...
                'HandleVisibility','off');  % <-- important
        end
    end
end

xlim(ax3, [0 max(0.12, OUT.params.BandHz(2)*1.25)]);
xlabel(ax3,'Frequency (Hz)');
ylabel(ax3,'Power');
title(ax3,'PSD of temporal motifs');
grid(ax3,'on');
box(ax3,'off');

legtxt = arrayfun(@(k) sprintf('Motif %d (%.1f%% | %.3f Hz)', ...
    k, explained(k), pcDomFreq(k)), 1:nMotifPlot, 'UniformOutput', false);

legend(ax3, legtxt, 'Location','northeast', 'Box','off');

% -----------------------------------------------------------------------
% D) neuron loading matrix
% -----------------------------------------------------------------------
ax4 = nexttile(tlo, [2 2]);
imagesc(ax4, scoresDisp);
xlabel(ax4,'Temporal motif');
ylabel(ax4,'Neurons (sorted)');
title(ax4,'Neuron loading matrix');
colormap(ax4, parula);

c = colorbar(ax4);
c.Label.String = 'z-scored neuron loading';

hold(ax4,'on');
cuts = find(diff(domID_sorted) ~= 0);
for i = 1:numel(cuts)
    yline(ax4, cuts(i)+0.5, '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 1);
end
box(ax4,'off');

for k = 1:numPCs
    patch(ax4, [k-0.5 k+0.5 k+0.5 k-0.5], [0.5 0.5 1.5 1.5], modeColors(k,:), ...
        'EdgeColor','none', 'FaceAlpha', 0.95);
end

% -----------------------------------------------------------------------
% E) dominant motif assignment in anatomy
% -----------------------------------------------------------------------
ax5 = nexttile(tlo, [2 1]); hold(ax5,'on');

if ~isempty(pos)
    if size(pos,2) < 3
        pos(:,3) = 0;
    end

    use3D = size(pos,2) >= 3 && any(abs(pos(:,3)) > eps);

    if use3D
        scatter3(ax5, pos(:,1), pos(:,2), pos(:,3), ...
            10, [0.86 0.86 0.86], 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);

        for k = 1:numPCs
            idx = (modeID == k);
            if any(idx)
                scatter3(ax5, pos(idx,1), pos(idx,2), pos(idx,3), ...
                    22, modeColors(k,:), 'filled');
            end
        end

        xlabel(ax5,'X');
        ylabel(ax5,'Y');
        zlabel(ax5,'Z');
        view(ax5, 35, 22);
    else
        scatter(ax5, pos(:,1), pos(:,2), ...
            10, [0.86 0.86 0.86], 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);

        for k = 1:numPCs
            idx = (modeID == k);
            if any(idx)
                scatter(ax5, pos(idx,1), pos(idx,2), ...
                    22, modeColors(k,:), 'filled');
            end
        end

        xlabel(ax5,'X');
        ylabel(ax5,'Y');
    end

    title(ax5,'Dominant temporal motif in anatomy');
    grid(ax5,'on');
    box(ax5,'off');

else
    axis(ax5,'off');
    text(ax5, 0.02, 0.95, 'No neuron positions provided.', ...
        'VerticalAlignment','top', 'FontName','Arial');
end

% -----------------------------------------------------------------------
% F) dominant motif counts  <-- replaces old text panel
% -----------------------------------------------------------------------
ax6 = nexttile(tlo, [1 1]); hold(ax6,'on');

for k = 1:numPCs
    bar(ax6, k, modeCount(k), ...
        'FaceColor', modeColors(k,:), ...
        'EdgeColor', 'none', ...
        'BarWidth', 0.8);
end

xlabel(ax6,'Temporal motif');
ylabel(ax6,'# neurons');
title(ax6,'Dominant motif counts');
grid(ax6,'on');
box(ax6,'off');
xlim(ax6, [0.3 numPCs+0.7]);

% Add counts on top of bars
yl = ylim(ax6);
for k = 1:numPCs
    text(ax6, k, modeCount(k) + 0.03*(yl(2)-yl(1)), ...
        sprintf('%d', modeCount(k)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',9);
end

% -----------------------------------------------------------------------
% compact summary annotation instead of a dedicated tile
% -----------------------------------------------------------------------
txt = sprintf([ ...
    'Kept neurons: %d\n' ...
    'Temporal motifs: %d\n' ...
    'Peak freqs (Hz): %s\n' ...
    'Timescale spread (CV): %.3f\n' ...
    'Mean PLV: %.3f\n' ...
    'Phase components: %d\n' ...
    'LDS pairs: %d\n' ...
    'Mean sliding pairs: %.3f\n' ...
    'Max sliding pairs: %.3f\n' ...
    'Interpretation: %s'], ...
    OUT.summary.nKeptNeurons, ...
    OUT.summary.numPCs, ...
    mat2str(OUT.summary.pcDomFreqHz,4), ...
    OUT.summary.pcFreqCV, ...
    OUT.summary.meanPLV, ...
    OUT.summary.nPhaseComponents, ...
    OUT.summary.nOscillatorPairsLDS, ...
    OUT.summary.meanSlidingPairs, ...
    OUT.summary.maxSlidingPairs, ...
    char(OUT.summary.interpretation));

annotation(figs.summary, 'textbox', [0.72 0.100 0.90 0.122], ...
    'String', txt, ...
    'FitBoxToText', 'on', ...
    'Interpreter', 'none', ...
    'FontName', 'Courier', ...
    'FontSize', 9, ...
    'EdgeColor', 'none', ...
    'BackgroundColor', 'w');

% ====================== PHASE DIAGNOSTICS ===============================
figs.phase = figure( ...
    'Color','w', ...
    'Position',[80 80 1450 900], ...
    'Name','Temporal motif phase diagnostics');

tlo2 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

axA = nexttile(tlo2); hold(axA,'on');
for k = 1:nMotifPlot
    plot(axA, tSec, OUT.phase.latentBand(:,k), ...
        'Color', modeColors(k,:), 'LineWidth', 1.0);
end
if ~isempty(eventTimesSec)
    yl = ylim(axA);
    for i = 1:numel(eventTimesSec)
        plot(axA, [eventTimesSec(i) eventTimesSec(i)], yl, ':', ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 0.6);
    end
end
xlabel(axA,'Time (s)');
ylabel(axA,'Bandpassed value');
title(axA,'Bandpassed temporal motifs');
grid(axA,'on');
box(axA,'off');

axB = nexttile(tlo2); hold(axB,'on');
for k = 1:min(numPCs,4)
    ph = unwrap(OUT.phase.latentPhase(:,k));
    plot(axB, tSec, ph, 'Color', modeColors(k,:), 'LineWidth', 1.0);
end
xlabel(axB,'Time (s)');
ylabel(axB,'Unwrapped phase (rad)');
title(axB,'Unwrapped phase');
grid(axB,'on');
box(axB,'off');

axC = nexttile(tlo2);
imagesc(axC, tSec, 1:numPCs, OUT.phase.latentAmp');
xlabel(axC,'Time (s)');
ylabel(axC,'Temporal motif');
title(axC,'Instantaneous amplitude');
colorbar(axC);
axis(axC,'tight');

axD = nexttile(tlo2);
imagesc(axD, PLV, [0 1]);
axis(axD,'image');
xlabel(axD,'Temporal motif');
ylabel(axD,'Temporal motif');
title(axD,'PLV matrix');
colorbar(axD);

% ====================== SLIDING DIAGNOSTICS =============================
figs.sliding = figure( ...
    'Color','w', ...
    'Position',[100 100 1500 850], ...
    'Name','Sliding temporal motif diagnostics');

tlo3 = tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');

axS1 = nexttile(tlo3); hold(axS1,'on');
plot(axS1, OUT.sliding.centerSec/60, OUT.sliding.meanPLV, '-o', ...
    'Color',[0.1 0.4 0.8], 'LineWidth', 1.2);
xlabel(axS1,'Window center (min)');
ylabel(axS1,'Mean PLV');
title(axS1,'Sliding PLV');
grid(axS1,'on');
box(axS1,'off');

axS2 = nexttile(tlo3); hold(axS2,'on');
plot(axS2, OUT.sliding.centerSec/60, OUT.sliding.pcFreqCV, '-o', ...
    'Color',[0.1 0.55 0.25], 'LineWidth', 1.2);
xlabel(axS2,'Window center (min)');
ylabel(axS2,'Frequency CV');
title(axS2,'Sliding timescale spread');
grid(axS2,'on');
box(axS2,'off');

axS3 = nexttile(tlo3); hold(axS3,'on');
plot(axS3, OUT.sliding.centerSec/60, OUT.sliding.nPairsLDS, '-o', ...
    'Color',[0.75 0.2 0.15], 'LineWidth', 1.2);
xlabel(axS3,'Window center (min)');
ylabel(axS3,'# LDS pairs');
title(axS3,'Sliding oscillator candidates');
grid(axS3,'on');
box(axS3,'off');

% ====================== LDS DIAGNOSTIC =================================
figs.lds = figure( ...
    'Color','w', ...
    'Position',[120 120 700 600], ...
    'Name','Temporal motif LDS diagnostics');

axL = axes(figs.lds); hold(axL,'on');
theta = linspace(0,2*pi,400);
plot(axL, cos(theta), sin(theta), 'k--', 'LineWidth', 0.8);
scatter(axL, real(lambda), imag(lambda), 35, [0.35 0.35 0.35], 'filled');
axis(axL,'equal');
xlabel(axL,'Real(\lambda)');
ylabel(axL,'Imag(\lambda)');
title(axL, sprintf('LDS eigenvalues | candidate pairs = %d', size(pairs,1)));
grid(axL,'on');
box(axL,'off');

end

%% =======================================================================
function modeColors = get_mode_colors(numModes)
base = [
    0.85 0.00 0.00
    0.90 0.35 0.00
    0.98 0.60 0.10
    0.95 0.78 0.15
    0.65 0.85 0.15
    0.25 0.80 0.65
    0.20 0.65 0.95
    0.35 0.45 0.95
    0.55 0.40 0.85
    0.35 0.10 0.40
    0.50 0.50 0.50
    0.20 0.20 0.20
    ];
if numModes <= size(base,1)
    modeColors = base(1:numModes,:);
else
    extra = lines(numModes - size(base,1));
    modeColors = [base; extra];
end
end
%% =======================================================================
function figs = make_time_mode_plots(OUT, fps)
% Redesigned motif plots for time-PCA.
%
% Focus:
%   - motif traces
%   - neuron loadings on motifs
%   - anatomy of dominant motif assignment
%   - motif-specific anatomy small multiples

figs = struct();
figs.modeMembership = [];
figs.modeDynamics   = [];
figs.neuronMap      = [];
figs.scoreSpace     = [];
figs = struct();

temporalPatterns = OUT.tpca.temporalPatterns;   % [time x K]
neuronScores     = OUT.tpca.neuronScores;       % [neurons x K]
modeID           = OUT.mode.modeID(:);
numPCs           = OUT.summary.numPCs;
tSec             = (0:size(temporalPatterns,1)-1) / fps;
eventTimesSec    = OUT.data.eventTimesSec(:);
pos              = OUT.data.neuronPositionsKept;

modeColors = get_mode_colors(numPCs);

% ---------------- dominant motif and sorting
[domAbs, domID] = max(abs(neuronScores), [], 2);
sortTable = [domID(:), -domAbs(:)];
[~, sortOrd] = sortrows(sortTable, [1 2]);

scoresSorted = neuronScores(sortOrd,:);
domID_sorted = domID(sortOrd);

scoresDisp = scoresSorted;
for k = 1:size(scoresDisp,2)
    x = scoresDisp(:,k);
    mu = mean(x,'omitnan');
    sd = std(x,0,'omitnan');
    if ~isfinite(sd) || sd == 0
        sd = 1;
    end
    scoresDisp(:,k) = (x - mu) ./ sd;
end
scoresDisp(~isfinite(scoresDisp)) = 0;

%% ===================== FIGURE 1: LOADING + TRACES ======================
figs.modeMembership = figure( ...
    'Color','w', ...
    'Position',[90 90 1600 950], ...
    'Name','Temporal motif membership');

tlo = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

% -----------------------------------------------------------------------
% A) neuron loading heatmap
% -----------------------------------------------------------------------
ax1 = nexttile(tlo, [2 1]);
imagesc(ax1, scoresDisp);
xlabel(ax1,'Temporal motif');
ylabel(ax1,'Neurons (sorted)');
title(ax1,'Neuron loading matrix');
colormap(ax1, parula);
c = colorbar(ax1);
c.Label.String = 'z-scored neuron loading';

hold(ax1,'on');
cuts = find(diff(domID_sorted) ~= 0);
for i = 1:numel(cuts)
    yline(ax1, cuts(i)+0.5, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
end
box(ax1,'off');

% add motif color strip across top
for k = 1:numPCs
    patch(ax1, [k-0.5 k+0.5 k+0.5 k-0.5], [0.5 0.5 1.5 1.5], modeColors(k,:), ...
        'EdgeColor','none', 'FaceAlpha', 0.95);
end

% -----------------------------------------------------------------------
% B) motif traces
% -----------------------------------------------------------------------
ax2 = nexttile(tlo, [1 1]); hold(ax2,'on');

nPlot = min(numPCs, 6);
offset = 0;
traceOffsets = zeros(nPlot,1);

for k = 1:nPlot
    x = temporalPatterns(:,k);
    s = std(x,0,'omitnan');
    if ~isfinite(s) || s == 0
        s = 1;
    end
    xdisp = x ./ s;

    traceOffsets(k) = offset;
    plot(ax2, tSec, xdisp + offset, 'Color', modeColors(k,:), 'LineWidth', 1.2);
    offset = offset + 3.2;
end

if ~isempty(eventTimesSec)
    yl = ylim(ax2);
    for i = 1:numel(eventTimesSec)
        plot(ax2, [eventTimesSec(i) eventTimesSec(i)], yl, ':', ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 0.6);
    end
end

yticks(ax2, traceOffsets);
yticklabels(ax2, compose('Motif %d', 1:nPlot));
xlabel(ax2,'Time (s)');
ylabel(ax2,'Temporal motifs');
title(ax2,'Temporal motif traces');
grid(ax2,'on');
box(ax2,'off');
axis(ax2,'tight');

% -----------------------------------------------------------------------
% C) dominant motif counts
% -----------------------------------------------------------------------
ax3 = nexttile(tlo, [1 1]); hold(ax3,'on');

counts = zeros(numPCs,1);
for k = 1:numPCs
    counts(k) = nnz(modeID == k);
end

for k = 1:numPCs
    bar(ax3, k, counts(k), 'FaceColor', modeColors(k,:), 'EdgeColor', 'none');
end
xlabel(ax3,'Temporal motif');
ylabel(ax3,'# neurons');
title(ax3,'Dominant motif assignment counts');
grid(ax3,'on');
box(ax3,'off');

%% ===================== FIGURE 2: ANATOMY ===============================
if ~isempty(pos)
    if size(pos,2) < 3
        pos(:,3) = 0;
    end

    figs.modeDynamics = figure( ...
        'Color','w', ...
        'Position',[110 110 1700 950], ...
        'Name','Temporal motif anatomy');

    tlo2 = tiledlayout(2, ceil((numPCs+1)/2), 'TileSpacing','compact', 'Padding','compact');

    % -------------------------------------------------------------------
    % A) combined anatomy map
    % -------------------------------------------------------------------
    axA = nexttile(tlo2); hold(axA,'on');

    scatter(axA, pos(:,1), pos(:,2), ...
        10, [0.85 0.85 0.85], 'filled', ...
        'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);

    for k = 1:numPCs
        idx = (modeID == k);
        if any(idx)
            scatter(axA, pos(idx,1), pos(idx,2), ...
                20, modeColors(k,:), 'filled');
        end
    end

    xlabel(axA,'X');
    ylabel(axA,'Y');
    title(axA,'Dominant temporal motif in anatomy');
    %view(axA, 35, 22);
    grid(axA,'on');
    box(axA,'off');

    % -------------------------------------------------------------------
    % B) motif-specific small multiples
    % -------------------------------------------------------------------
    for k = 1:numPCs
        ax = nexttile(tlo2); hold(ax,'on');
        idx = (modeID == k);

        scatter(ax, pos(:,1), pos(:,2), ...
            8, [0.85 0.85 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35);

        if any(idx)
            scatter(ax, pos(idx,1), pos(idx,2), ...
                18, modeColors(k,:), 'filled');
        end

        xlabel(ax,'X');
        ylabel(ax,'Y');
        title(ax, sprintf('Motif %d', k), 'Color', modeColors(k,:));
        %view(ax, 35, 22);
        grid(ax,'on');
        box(ax,'off');
    end

    % keep compatibility with your old field naming
    figs.neuronMap = figs.modeDynamics;

else
    figs.modeDynamics = [];
    figs.neuronMap = [];
end

%% ===================== FIGURE 3: SCORE SPACE ===========================
figs.scoreSpace = figure( ...
    'Color','w', ...
    'Position',[130 130 1700 900], ...
    'Name','Neuron score space by temporal motif');

Y = neuronScores(:, 1:min(3,size(neuronScores,2)));

if size(Y,2) >= 3
    tlo3 = tiledlayout(2, ceil(numPCs/2), 'TileSpacing','compact', 'Padding','compact');

    for k = 1:numPCs
        ax = nexttile(tlo3); hold(ax,'on');
        idx = (modeID == k);

        scatter(ax, Y(~idx,1), Y(~idx,2), ...
            10, [0.85 0.85 0.85], 'filled');
        scatter(ax, Y(idx,1), Y(idx,2), ...
            18, modeColors(k,:), 'filled');

        xlabel(ax,'Score 1');
        ylabel(ax,'Score 2');
        title(ax, sprintf('Neurons assigned to motif %d', k));
        grid(ax,'on');
        box(ax,'off');
    end
else
    figs.scoreSpace = [];
end

end
