function OUT = oscillator_test_pipeline_v1(dff_full, fps, varargin)
%OSCILLATOR_TEST_PIPELINE_V1
% Test whether population dynamics are better explained by one oscillator
% or multiple oscillatory modes in latent space.
%
% INPUT
%   dff_full : [neurons x time] matrix
%   fps      : sampling rate in Hz
%
% NAME-VALUE PARAMETERS
%   'EpochFrames'        : [start end] frame range (default full recording)
%   'KeepMask'           : logical [neurons x 1] mask of neurons to keep (default all finite-enough)
%   'FiniteFracThresh'   : minimum fraction finite per neuron (default 0.95)
%   'UseZScore'          : z-score each neuron before PCA (default true)
%   'NumPCs'             : number of latent PCs to keep (default 6)
%   'BandHz'             : band for phase analyses (default [0.01 0.1])
%   'WelchWinSec'        : window length for PSD (default 120 s)
%   'WelchOverlapFrac'   : overlap fraction for PSD (default 0.5)
%   'LDSReg'             : ridge regularization for latent dynamics fit (default 1e-3)
%   'MinEigMag'          : minimum |lambda| for oscillatory mode candidate (default 0.90)
%   'MaxEigMag'          : maximum |lambda| for oscillatory mode candidate (default 1.05)
%   'DoPlots'            : true/false (default true)
%   'Verbose'            : true/false (default true)
%
% OUTPUT
%   OUT.data
%   OUT.pca
%   OUT.psd
%   OUT.phase
%   OUT.lds
%   OUT.summary
%   OUT.fig
%
% NOTES
%   - Complex conjugate eigenvalue pairs of the latent dynamics matrix A
%     are interpreted as candidate oscillatory modes.
%   - Phase-locking is computed after bandpass filtering latent PCs in BandHz.
%   - This is an exploratory dynamical test, not a formal model comparison.

%% -------------------- Parse inputs --------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'dff_full', @(x) isnumeric(x) && ndims(x)==2);
addRequired(ip, 'fps', @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'EpochFrames', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'KeepMask', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(ip, 'FiniteFracThresh', 0.95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(ip, 'UseZScore', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'NumPCs', 6, @(x) isnumeric(x) && isscalar(x) && x >= 3);

addParameter(ip, 'BandHz', [0.01 0.1], @(x) isnumeric(x) && numel(x)==2 && x(1)>0 && x(2)>x(1));
addParameter(ip, 'WelchWinSec', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'WelchOverlapFrac', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);

addParameter(ip, 'LDSReg', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(ip, 'MinEigMag', 0.90, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'MaxEigMag', 1.05, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'DoPlots', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(ip, dff_full, fps, varargin{:});
P = ip.Results;

[N0, T0] = size(dff_full);

%% -------------------- Epoch selection --------------------
if isempty(P.EpochFrames)
    epochFrames = [1 T0];
else
    epochFrames = round(P.EpochFrames(:)');
    epochFrames(1) = max(1, epochFrames(1));
    epochFrames(2) = min(T0, epochFrames(2));
end

tIdx = epochFrames(1):epochFrames(2);
dff_epoch = dff_full(:, tIdx);
T = size(dff_epoch, 2);

%% -------------------- Neuron selection --------------------
fracFinite = mean(isfinite(dff_epoch), 2);

if isempty(P.KeepMask)
    keep = fracFinite >= P.FiniteFracThresh;
else
    keep = logical(P.KeepMask(:));
    if numel(keep) ~= N0
        error('KeepMask must have one entry per neuron in dff_full.');
    end
    keep = keep & (fracFinite >= P.FiniteFracThresh);
end

keepIdx = find(keep);
dropIdx = find(~keep);

if isempty(keepIdx)
    error('No neurons survive selection.');
end

Xraw = dff_epoch(keepIdx, :);

if P.UseZScore
    X = zscore_fill_rows(Xraw);
else
    X = Xraw;
    X(~isfinite(X)) = 0;
end

if P.Verbose
    fprintf('\n=== oscillator_test_pipeline_v1 ===\n');
    fprintf('Input neurons               : %d\n', N0);
    fprintf('Epoch frames                : %d to %d (%d frames)\n', epochFrames(1), epochFrames(2), T);
    fprintf('Kept neurons                : %d\n', numel(keepIdx));
    fprintf('Band for phase analyses     : [%.3f %.3f] Hz\n', P.BandHz(1), P.BandHz(2));
end

%% -------------------- PCA --------------------
numPCs = min([P.NumPCs, size(X,1), size(X,2)]);
% observations = time, variables = neurons
maxPC = min([numPCs, rank(X'), size(X,1), size(X,2)]);
[coeff, score, latent, ~, explained, mu] = pca(X', 'NumComponents', maxPC);
numPCs = maxPC;
latentTS = score(:, 1:numPCs);

% score: [time x PCs]
latentTS = score(:, 1:numPCs);  % T x numPCs

%% -------------------- PSD of latent PCs --------------------
[pcDomFreq, pcPeakPower, fAxis, PxxPC] = compute_pc_psd(latentTS, fps, P.WelchWinSec, P.WelchOverlapFrac, P.BandHz);

%% -------------------- Bandpass + Hilbert phase --------------------
latentBand = bandpass_latent(latentTS, fps, P.BandHz);
latentHilb = hilbert(latentBand);
latentPhase = angle(latentHilb);
latentAmp = abs(latentHilb);

% PLV matrix
% PLV matrix
PLV = nan(numPCs, numPCs);
for i = 1:numPCs
    for j = 1:numPCs
        dphi = latentPhase(:,i) - latentPhase(:,j);
        good = isfinite(dphi);
        if nnz(good) < 10
            PLV(i,j) = NaN;
        else
            PLV(i,j) = abs(mean(exp(1i * dphi(good))));
        end
    end
end

%% -------------------- Fit linear latent dynamics --------------------
% X1 -> X2
X1 = latentTS(1:end-1, :)';   % [pcs x time-1]
X2 = latentTS(2:end,   :)';   % [pcs x time-1]

% ridge-stabilized least-squares: A = X2*X1'/(X1*X1' + reg*I)
reg = P.LDSReg;
A = (X2 * X1') / (X1 * X1' + reg * eye(size(X1,1)));

[V,D] = eig(A);
lambda = diag(D);

eigMag = abs(lambda);
eigAng = angle(lambda);                 % radians / sample
eigFreqHz = abs(eigAng) * fps / (2*pi); % cycles / second

% Candidate oscillatory modes = complex eigenvalues with magnitude near unit circle
isComplex = abs(imag(lambda)) > 1e-10;
isCandidate = isComplex & eigMag >= P.MinEigMag & eigMag <= P.MaxEigMag & ...
    eigFreqHz >= P.BandHz(1) & eigFreqHz <= P.BandHz(2);

% Count unique conjugate pairs
oscIdx = find(isCandidate);
pairMask = false(size(oscIdx));
pairs = [];

for ii = 1:numel(oscIdx)
    if pairMask(ii), continue; end
    k = oscIdx(ii);
    lam = lambda(k);

    % find best conjugate partner among remaining candidates
    remIdx = oscIdx(~pairMask);
    remVals = lambda(remIdx);
    [~, jjRel] = min(abs(remVals - conj(lam)));
    partner = remIdx(jjRel);

    if partner ~= k
        pairMask(oscIdx == k) = true;
        pairMask(oscIdx == partner) = true;
        pairs = [pairs; k partner]; %#ok<AGROW>
    else
        pairMask(oscIdx == k) = true;
    end
end

% summarize pair frequencies/magnitudes
pairFreqHz = nan(size(pairs,1),1);
pairMag = nan(size(pairs,1),1);
for i = 1:size(pairs,1)
    pairFreqHz(i) = mean(eigFreqHz(pairs(i,:)), 'omitnan');
    pairMag(i) = mean(eigMag(pairs(i,:)), 'omitnan');
end

%% -------------------- Summary metrics --------------------
% frequency diversity across latent PCs
pcFreqCV = std(pcDomFreq, 'omitnan') / mean(pcDomFreq, 'omitnan');

% mean off-diagonal PLV
maskOff = triu(true(numPCs),1);
plvVals = PLV(maskOff);
plvVals = plvVals(isfinite(plvVals));
if isempty(plvVals)
    meanPLV = NaN;
else
    meanPLV = mean(plvVals);
end

% phase-lock graph components heuristic
plvThresh = 0.7;
adj = false(numPCs, numPCs);
goodPLV = isfinite(PLV);
adj(goodPLV) = PLV(goodPLV) >= plvThresh;
adj(1:numPCs+1:end) = false;
nPhaseComponents = count_connected_components(adj);

summary = struct();
summary.nKeptNeurons = numel(keepIdx);
summary.numPCs = numPCs;
summary.pcDomFreqHz = pcDomFreq;
summary.pcFreqCV = pcFreqCV;
summary.meanPLV = meanPLV;
summary.plvThreshold = plvThresh;
summary.nPhaseComponents = nPhaseComponents;
summary.nOscillatorPairsLDS = size(pairs,1);
summary.ldsPairFreqHz = pairFreqHz;
summary.ldsPairMagnitude = pairMag;

% coarse interpretation
if size(pairs,1) <= 1 && meanPLV > 0.8
    summary.interpretation = "More consistent with one dominant oscillator + noise.";
elseif size(pairs,1) >= 2 && (pcFreqCV > 0.15 || meanPLV < 0.7 || nPhaseComponents >= 2)
    summary.interpretation = "More consistent with multiple oscillatory modes.";
else
    summary.interpretation = "Mixed evidence; possible weakly coupled or partially overlapping oscillators.";
end

if P.Verbose
    fprintf('Latent PC dominant freqs (Hz): %s\n', mat2str(pcDomFreq, 4));
    fprintf('PC dominant frequency CV     : %.4f\n', pcFreqCV);
    fprintf('Mean off-diagonal PLV        : %.4f\n', meanPLV);
    fprintf('Phase-lock components        : %d (threshold %.2f)\n', nPhaseComponents, plvThresh);
    fprintf('LDS oscillatory eigen-pairs  : %d\n', size(pairs,1));
    if ~isempty(pairFreqHz)
        fprintf('LDS pair freqs (Hz)          : %s\n', mat2str(pairFreqHz, 4));
    end
    fprintf('Interpretation               : %s\n', summary.interpretation);
end

%% -------------------- Output --------------------
OUT = struct();

OUT.params = P;

OUT.data = struct();
OUT.data.keepIdx = keepIdx;
OUT.data.dropIdx = dropIdx;
OUT.data.epochFrames = epochFrames;
OUT.data.Xraw = Xraw;
OUT.data.Xz = X;

OUT.pca = struct();
OUT.pca.coeff = coeff;
OUT.pca.score = score;
OUT.pca.latentTS = latentTS;
OUT.pca.latent = latent;
OUT.pca.explained = explained;
OUT.pca.mu = mu;

OUT.psd = struct();
OUT.psd.fAxis = fAxis;
OUT.psd.PxxPC = PxxPC;
OUT.psd.pcDomFreqHz = pcDomFreq;
OUT.psd.pcPeakPower = pcPeakPower;

OUT.phase = struct();
OUT.phase.latentBand = latentBand;
OUT.phase.latentHilb = latentHilb;
OUT.phase.latentPhase = latentPhase;
OUT.phase.latentAmp = latentAmp;
OUT.phase.PLV = PLV;

OUT.lds = struct();
OUT.lds.A = A;
OUT.lds.V = V;
OUT.lds.D = D;
OUT.lds.lambda = lambda;
OUT.lds.eigMag = eigMag;
OUT.lds.eigAng = eigAng;
OUT.lds.eigFreqHz = eigFreqHz;
OUT.lds.isCandidate = isCandidate;
OUT.lds.pairs = pairs;
OUT.lds.pairFreqHz = pairFreqHz;
OUT.lds.pairMag = pairMag;

OUT.summary = summary;
OUT.fig = struct();

%% -------------------- Plots --------------------
if P.DoPlots
    OUT.fig = make_oscillator_plots(OUT, fps);
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
% X: [time x pcs]
% Butterworth fallback, zero-phase
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
% simple undirected graph connected components count
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
function figs = make_oscillator_plots(OUT, fps)
figs = struct();

latentTS = OUT.pca.latentTS;
explained = OUT.pca.explained;
fAxis = OUT.psd.fAxis;
PxxPC = OUT.psd.PxxPC;
pcDomFreq = OUT.psd.pcDomFreqHz;
PLV = OUT.phase.PLV;
lambda = OUT.lds.lambda;
eigMag = OUT.lds.eigMag;
eigFreqHz = OUT.lds.eigFreqHz;
isCand = OUT.lds.isCandidate;
pairs = OUT.lds.pairs;

numPCs = size(latentTS,2);
tSec = (0:size(latentTS,1)-1) / fps;

% -------- Summary figure --------
figs.summary = figure('Color','w','Position',[90 70 1450 900], 'Name','Oscillator test summary');
tlo = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

% 3D trajectory
ax1 = nexttile(tlo, [2 1]); hold(ax1,'on');
if size(latentTS,2) >= 3
    plot3(ax1, latentTS(:,1), latentTS(:,2), latentTS(:,3), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    xlabel(ax1,'PC1'); ylabel(ax1,'PC2'); zlabel(ax1,'PC3');
    title(ax1, sprintf('PC trajectory (PC1-3 = %.1f%%)', sum(explained(1:min(3,end)))));
    view(ax1,45,30); grid(ax1,'on'); box(ax1,'off');
else
    plot(ax1, latentTS(:,1), latentTS(:,2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    xlabel(ax1,'PC1'); ylabel(ax1,'PC2');
    title(ax1,'PC trajectory');
end

% scree
ax2 = nexttile(tlo);
bar(ax2, explained(1:min(10,numel(explained))));
xlabel(ax2,'PC'); ylabel(ax2,'% variance explained');
title(ax2,'Explained variance');

% PSDs
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

% eigenvalues
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

% PLV matrix
ax5 = nexttile(tlo);
imagesc(ax5, PLV, [0 1]);
axis(ax5,'image');
xlabel(ax5,'PC'); ylabel(ax5,'PC');
title(ax5,'Phase-locking value (PLV)');
colorbar(ax5);

% summary text
ax6 = nexttile(tlo);
axis(ax6,'off');
txt = {
    sprintf('Kept neurons: %d', OUT.summary.nKeptNeurons)
    sprintf('Latent PCs: %d', OUT.summary.numPCs)
    sprintf('PC dom freqs (Hz): %s', mat2str(OUT.summary.pcDomFreqHz, 3))
    sprintf('PC freq CV: %.3f', OUT.summary.pcFreqCV)
    sprintf('Mean PLV: %.3f', OUT.summary.meanPLV)
    sprintf('Phase components: %d', OUT.summary.nPhaseComponents)
    sprintf('LDS oscillator pairs: %d', OUT.summary.nOscillatorPairsLDS)
    sprintf('LDS pair freqs (Hz): %s', mat2str(OUT.summary.ldsPairFreqHz, 3))
    ['Interpretation: ' char(OUT.summary.interpretation)]
    };
text(ax6, 0.01, 0.98, txt, 'VerticalAlignment','top', 'FontName','Courier');

% -------- phase-difference figure --------
figs.phase = figure('Color','w','Position',[120 100 1200 800], 'Name','Phase coupling details');
tlo2 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

axA = nexttile(tlo2); hold(axA,'on');
for i = 1:numPCs
    plot(axA, tSec, OUT.phase.latentBand(:,i), 'LineWidth', 1);
end
xlabel(axA,'Time (s)');
ylabel(axA,'Bandpassed latent signal');
title(axA,'Bandpassed PCs');
legend(axA, compose('PC%d', 1:numPCs), 'Location','best');

axB = nexttile(tlo2); hold(axB,'on');
for i = 1:min(numPCs,4)
    plot(axB, tSec, unwrap(OUT.phase.latentPhase(:,i)), 'LineWidth', 1);
end
xlabel(axB,'Time (s)');
ylabel(axB,'Unwrapped phase (rad)');
title(axB,'Unwrapped phases (first PCs)');
legend(axB, compose('PC%d', 1:min(numPCs,4)), 'Location','best');

axC = nexttile(tlo2);
imagesc(axC, OUT.phase.latentAmp');
xlabel(axC,'Time');
ylabel(axC,'PC');
title(axC,'Instantaneous amplitude');
colorbar(axC);

axD = nexttile(tlo2);
imagesc(axD, OUT.phase.PLV, [0 1]);
axis(axD,'image');
xlabel(axD,'PC'); ylabel(axD,'PC');
title(axD,'PLV matrix');
colorbar(axD);
end