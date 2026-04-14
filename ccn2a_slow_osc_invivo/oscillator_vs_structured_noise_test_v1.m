function OUT = oscillator_vs_structured_noise_test_v1(baseFolder, varargin)
% OSCILLATOR_VS_STRUCTURED_NOISE_TEST_V1
%
% Purpose
% -------
% Test whether ultraslow activity is better described as:
%   (1) a genuine oscillator / quasi-oscillator
%   (2) structured slow noise with apparent rhythmicity
%
% Strategy
% --------
% Run the analysis PER RECORDING (recommended first), then summarize across fish.
% For each recording, compare the real data against surrogate controls using:
%   - PSD peak height in the target band
%   - PSD peak sharpness / Q-like factor
%   - ACF peak prominence and peak count
%   - Sliding-window frequency stability
%   - Instantaneous-frequency regularity
%
% Surrogates
% ----------
%   A. phase-randomized surrogates (preserve PSD globally, destroy temporal phase structure)
%   B. AR(1) surrogates (preserve slow autocorrelation, produce structured noise null)
%
% Inputs
% ------
% baseFolder : root containing fish/session folders
%
% Name-value parameters
% ---------------------
% 'ResultsFileName'      : default 'dffs_repact_respcells.mat'
% 'MetaFileName'         : default 'metadata_multimodal.mat'
% 'DffVarCandidates'     : variable names for dff matrix
% 'FpsVarCandidates'     : variable names for fps
% 'EpochFrames'          : [] or [start end]
% 'FiniteFracThresh'     : default 0.95
% 'UseZScore'            : default true
% 'EnableAcfPeakFilter'  : default true
% 'AcfPeakLagRangeSec'   : default [10 120]
% 'AcfMinProminence'     : default 0.08
% 'MinNumAcPeaks'        : default 2
% 'AcfMaxLagSec'         : default 180
% 'BandHz'               : default [0.01 0.10]
% 'WelchWinSec'          : default 120
% 'WelchOverlapFrac'     : default 0.5
% 'NumPCs'               : default 6
% 'WinSec'               : default 300
% 'StepSec'              : default 60
% 'NumSurrogates'        : default 100
% 'RepresentativeSignal' : 'pc1' | 'mean' | 'median_neuron'
% 'DoPlots'              : default true
% 'SaveFigures'          : default true
% 'Verbose'              : default true
%
% Output
% ------
% OUT.summaryTable       : one row per recording
% OUT.groupTable         : aggregate summary of real-vs-surrogate evidence
% OUT.AllRes             : detailed per-recording results
%
% Notes
% -----
% The code is intentionally per-fish first. That is the right level to test
% whether oscillator-like structure is reproducible before pooling.
%
% Publication-useful figures produced:
%   Fig 1. Representative per-recording real vs surrogate comparison
%   Fig 2. Per-recording effect-size summary across fish
%   Fig 3. Fraction of recordings exceeding surrogate thresholds
%   Fig 4. Real vs surrogate distributions for key metrics
%
% ---------------------------------------------------------------------

%% ============================ PARSE ===================================
ip = inputParser;
addRequired(ip, 'baseFolder', @(x) ischar(x) || isstring(x));

addParameter(ip, 'ResultsFileName', 'dffs_repact_respcells.mat', @(x)ischar(x)||isstring(x));
addParameter(ip, 'MetaFileName', 'metadata_multimodal.mat', @(x)ischar(x)||isstring(x));
addParameter(ip, 'DffVarCandidates', {'dff_new','dff','DFF','dffMat'}, @iscell);
addParameter(ip, 'FpsVarCandidates', {'fps','fs','Fs','origFps'}, @iscell);

addParameter(ip, 'EpochMinutes', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'FiniteFracThresh', 0.95, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
addParameter(ip, 'UseZScore', true, @(x) islogical(x) || isnumeric(x));

addParameter(ip, 'EnableAcfPeakFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'AcfPeakLagRangeSec', [10 120], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'AcfMinProminence', 0.08, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(ip, 'MinNumAcPeaks', 2, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(ip, 'AcfMaxLagSec', 180, @(x) isnumeric(x) && isscalar(x) && x>0);

addParameter(ip, 'BandHz', [0.01 0.10], @(x) isnumeric(x) && numel(x)==2 && x(1)>0 && x(2)>x(1));
addParameter(ip, 'WelchWinSec', 120, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(ip, 'WelchOverlapFrac', 0.5, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<1);
addParameter(ip, 'NumPCs', 6, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(ip, 'WinSec', 300, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(ip, 'StepSec', 60, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(ip, 'NumSurrogates', 100, @(x) isnumeric(x) && isscalar(x) && x>=20);
addParameter(ip, 'RepresentativeSignal', 'pc1', @(x) ischar(x) || isstring(x));

addParameter(ip, 'DoPlots', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'SaveFigures', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(ip, baseFolder, varargin{:});
P = ip.Results;
P.baseFolder = char(P.baseFolder);

%% =========================== OUTPUT DIR ===============================
outDir = fullfile(P.baseFolder, 'oscillator_vs_noise_results_v1');
figDir = fullfile(outDir, 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% ====================== FIND VALID FOLDERS ============================
allFolders = strsplit(genpath(P.baseFolder), pathsep);
validFolders = {};
resultsFiles = {};
metaFiles = {};
fishLabel = {};

for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    rf = fullfile(allFolders{i}, P.ResultsFileName);
    mf = fullfile(allFolders{i}, P.MetaFileName);
    if exist(rf,'file')==2 && exist(mf,'file')==2
        validFolders{end+1} = allFolders{i}; %#ok<AGROW>
        resultsFiles{end+1} = rf; %#ok<AGROW>
        metaFiles{end+1} = mf; %#ok<AGROW>
        [~, fishLabel{end+1}] = fileparts(allFolders{i}); %#ok<AGROW>
    end
end

if isempty(validFolders)
    error('No valid folders found under %s', P.baseFolder);
end

if P.Verbose
    fprintf('Found %d recordings. Running per-recording oscillator vs surrogate test.\n', numel(validFolders));
end

%% ============================ RUN ====================================
AllRes = struct([]);
rows = [];

for iRec = 1:numel(validFolders)
    if P.Verbose
        fprintf('\n============================================================\n');
        fprintf('Recording %d / %d\n', iRec, numel(validFolders));
        fprintf('Folder: %s\n', validFolders{iRec});
    end

    try
        R = load(resultsFiles{iRec});
        [dff, fps] = extract_dff_and_fps(R, P.DffVarCandidates, P.FpsVarCandidates, resultsFiles{iRec});

        if isempty(P.EpochMinutes)
            epochFrames = [1 size(dff,2)];
        else
            startMin = P.EpochMinutes(1);
            endMin   = P.EpochMinutes(2);
        
            startFrame = max(1, round(startMin * 60 * fps));
            endFrame   = min(size(dff,2), round(endMin * 60 * fps));
        
            if endFrame <= startFrame
                warning('Invalid epoch (too short). Using full recording.');
                epochFrames = [1 size(dff,2)];
            else
                epochFrames = [startFrame endFrame];
            end
        end

        Res = analyze_one_recording(dff, fps, epochFrames, P);
        Res.recordingLabel = fishLabel{iRec};
        Res.folder = validFolders{iRec};
        AllRes = [AllRes; Res]; %#ok<AGROW>

        row = struct();
        row.recordingLabel = string(fishLabel{iRec});
        row.folder = string(validFolders{iRec});
        row.nNeurons = size(dff,1);
        row.nFrames = size(dff,2);
        row.fps = fps;
        row.durationSec = size(dff,2)/fps;
        row.nKept = Res.nKept;

        row.real_peakFreq = Res.real.peakFreq;
        row.real_peakPower = Res.real.peakPower;
        row.real_qFactor = Res.real.qFactor;
        row.real_acfProm = Res.real.acfProm;
        row.real_acfPeakCount = Res.real.acfPeakCount;
        row.real_freqCVwindows = Res.real.freqCVwindows;
        row.real_instFreqCV = Res.real.instFreqCV;

        row.z_peakPower_phaseRand = Res.zscores.phaseRand.peakPower;
        row.z_qFactor_phaseRand = Res.zscores.phaseRand.qFactor;
        row.z_acfProm_phaseRand = Res.zscores.phaseRand.acfProm;
        row.z_freqCVwindows_phaseRand = Res.zscores.phaseRand.freqCVwindows_neg;
        row.z_instFreqCV_phaseRand = Res.zscores.phaseRand.instFreqCV_neg;

        row.z_peakPower_ar1 = Res.zscores.ar1.peakPower;
        row.z_qFactor_ar1 = Res.zscores.ar1.qFactor;
        row.z_acfProm_ar1 = Res.zscores.ar1.acfProm;
        row.z_freqCVwindows_ar1 = Res.zscores.ar1.freqCVwindows_neg;
        row.z_instFreqCV_ar1 = Res.zscores.ar1.instFreqCV_neg;

        row.pass_peakPower_phaseRand = Res.pass.phaseRand.peakPower;
        row.pass_qFactor_phaseRand = Res.pass.phaseRand.qFactor;
        row.pass_acfProm_phaseRand = Res.pass.phaseRand.acfProm;
        row.pass_freqStability_phaseRand = Res.pass.phaseRand.freqStability;
        row.pass_instFreqReg_phaseRand = Res.pass.phaseRand.instFreqReg;

        row.pass_peakPower_ar1 = Res.pass.ar1.peakPower;
        row.pass_qFactor_ar1 = Res.pass.ar1.qFactor;
        row.pass_acfProm_ar1 = Res.pass.ar1.acfProm;
        row.pass_freqStability_ar1 = Res.pass.ar1.freqStability;
        row.pass_instFreqReg_ar1 = Res.pass.ar1.instFreqReg;

        row.oscillatorEvidenceScore = Res.oscillatorEvidenceScore;
        row.classification = string(Res.classification);

        rows = [rows; row]; %#ok<AGROW>

        if P.DoPlots && P.SaveFigures && isfield(Res, 'fig')
            fn = fullfile(figDir, [fishLabel{iRec} '_representative.png']);
            exportgraphics(Res.fig.representative, fn, 'Resolution', 200);
            close(Res.fig.representative);
        end

    catch ME
        warning('Failed on folder:\n%s\n\n%s', validFolders{iRec}, getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

if isempty(AllRes)
    error('No recordings completed successfully.');
end

summaryTable = struct2table(rows);
groupTable = build_group_summary(summaryTable);

writetable(summaryTable, fullfile(outDir, 'per_recording_summary.csv'));
writetable(groupTable, fullfile(outDir, 'group_summary.csv'));
save(fullfile(outDir, 'summary_struct.mat'), 'AllRes', 'summaryTable', 'groupTable', 'P', '-v7.3');

OUT = struct();
OUT.AllRes = AllRes;
OUT.summaryTable = summaryTable;
OUT.groupTable = groupTable;
OUT.outDir = outDir;
OUT.fig = struct();

if P.DoPlots
    OUT.fig = make_group_figures(summaryTable, groupTable, figDir, P.SaveFigures);
end

end

%% =====================================================================
function Res = analyze_one_recording(dff_full, fps, epochFrames, P)
% -------- epoch and QC --------
tIdx = epochFrames(1):epochFrames(2);
dff_epoch = dff_full(:, tIdx);
fracFinite = mean(isfinite(dff_epoch), 2);
keepBase = fracFinite >= P.FiniteFracThresh;

keepAcf = true(size(keepBase));
acfPeakCountAll = nan(size(keepBase));
if P.EnableAcfPeakFilter
    maxLag = round(P.AcfMaxLagSec * fps);
    lagAxis = (-maxLag:maxLag)/fps;
    lagMask = lagAxis >= P.AcfPeakLagRangeSec(1) & lagAxis <= P.AcfPeakLagRangeSec(2);

    for i = 1:size(dff_epoch,1)
        if ~keepBase(i)
            keepAcf(i) = false;
            continue;
        end
        xi = dff_epoch(i,:);
        xi = xi(isfinite(xi));
        if numel(xi) < maxLag + 5
            keepAcf(i) = false;
            acfPeakCountAll(i) = 0;
            continue;
        end
        xi = xi - mean(xi);
        aci = xcorr(xi, maxLag, 'coeff');
        aciWin = aci(lagMask);
        pks = findpeaks(aciWin, 'MinPeakProminence', P.AcfMinProminence);
        acfPeakCountAll(i) = numel(pks);
        keepAcf(i) = acfPeakCountAll(i) >= P.MinNumAcPeaks;
    end
end

keep = keepBase & keepAcf;
Xraw = dff_epoch(keep,:);
if isempty(Xraw)
    error('No neurons survived filtering in selected epoch.');
end

if P.UseZScore
    X = zscore_fill_rows(Xraw);
else
    X = Xraw;
    X(~isfinite(X)) = 0;
end

% -------- representative signal --------
rep = get_representative_signal(X, fps, P);

% -------- real metrics --------
realM = compute_signal_metrics(rep, fps, P);

% -------- surrogates --------
S_phase = nan(P.NumSurrogates, 6);
S_ar1 = nan(P.NumSurrogates, 6);
for s = 1:P.NumSurrogates
    x_phase = make_phase_randomized(rep);
    x_ar1 = make_ar1_surrogate(rep);

    m1 = compute_signal_metrics(x_phase, fps, P);
    m2 = compute_signal_metrics(x_ar1, fps, P);

    S_phase(s,:) = [m1.peakPower, m1.qFactor, m1.acfProm, m1.acfPeakCount, m1.freqCVwindows, m1.instFreqCV];
    S_ar1(s,:)   = [m2.peakPower, m2.qFactor, m2.acfProm, m2.acfPeakCount, m2.freqCVwindows, m2.instFreqCV];
end

% -------- z-like evidence --------
Zphase = compute_zscores(realM, S_phase);
Zar1   = compute_zscores(realM, S_ar1);

% -------- thresholded pass/fail --------
passPhase = struct();
passPhase.peakPower    = realM.peakPower    > prctile(S_phase(:,1),95);
passPhase.qFactor      = realM.qFactor      > prctile(S_phase(:,2),95);
passPhase.acfProm      = realM.acfProm      > prctile(S_phase(:,3),95);
passPhase.freqStability= realM.freqCVwindows < prctile(S_phase(:,5),5);
passPhase.instFreqReg  = realM.instFreqCV   < prctile(S_phase(:,6),5);

passAR1 = struct();
passAR1.peakPower    = realM.peakPower    > prctile(S_ar1(:,1),95);
passAR1.qFactor      = realM.qFactor      > prctile(S_ar1(:,2),95);
passAR1.acfProm      = realM.acfProm      > prctile(S_ar1(:,3),95);
passAR1.freqStability= realM.freqCVwindows < prctile(S_ar1(:,5),5);
passAR1.instFreqReg  = realM.instFreqCV   < prctile(S_ar1(:,6),5);

% -------- combined score --------
oscScore = mean([
    Zphase.peakPower, Zphase.qFactor, Zphase.acfProm, Zphase.freqCVwindows_neg, Zphase.instFreqCV_neg, ...
    Zar1.peakPower,   Zar1.qFactor,   Zar1.acfProm,   Zar1.freqCVwindows_neg,   Zar1.instFreqCV_neg], 'omitnan');

nStrongPass = sum(struct2array(passPhase)) + sum(struct2array(passAR1));
if nStrongPass >= 8 && oscScore > 1.5
    classification = "Oscillator-like";
elseif nStrongPass >= 4 && oscScore > 0.5
    classification = "Quasi-oscillatory / structured rhythmic dynamics";
else
    classification = "Structured slow noise / weak rhythmic structure";
end

Res = struct();
Res.nKept = size(X,1);
Res.repSignal = rep;
Res.real = realM;
Res.surrogates.phaseRand = S_phase;
Res.surrogates.ar1 = S_ar1;
Res.zscores.phaseRand = Zphase;
Res.zscores.ar1 = Zar1;
Res.pass.phaseRand = passPhase;
Res.pass.ar1 = passAR1;
Res.oscillatorEvidenceScore = oscScore;
Res.classification = classification;
Res.fig = struct();
Res.fig.representative = make_representative_figure(rep, realM, S_phase, S_ar1, fps, P, classification);
end

%% =====================================================================
function rep = get_representative_signal(X, fps, P)
switch lower(string(P.RepresentativeSignal))
    case "pc1"
        numPC = min(P.NumPCs, min(size(X,1), size(X,2)));
        [~, score] = pca(X', 'NumComponents', numPC);
        rep = score(:,1);
    case "mean"
        rep = mean(X,1,'omitnan')';
    case "median_neuron"
        [~, score] = pca(X', 'NumComponents', min(P.NumPCs, min(size(X,1), size(X,2))));
        [domFreqs,~,~,~] = compute_pc_psd(score(:,1:min(3,size(score,2))), fps, P.WelchWinSec, P.WelchOverlapFrac, P.BandHz);
        targetF = median(domFreqs,'omitnan');
        neuronFreqs = nan(size(X,1),1);
        for i = 1:size(X,1)
            m = compute_signal_metrics(X(i,:)', fps, P);
            neuronFreqs(i) = m.peakFreq;
        end
        [~,ix] = min(abs(neuronFreqs - targetF));
        rep = X(ix,:)';
    otherwise
        error('Unknown RepresentativeSignal: %s', P.RepresentativeSignal);
end
rep = rep(:);
rep = fillmissing(rep, 'linear', 'EndValues', 'nearest');
rep = detrend(rep);
end

%% =====================================================================
function M = compute_signal_metrics(x, fps, P)
x = x(:);
x = fillmissing(x, 'linear', 'EndValues', 'nearest');

% PSD and peak power
nwin = max(16, round(P.WelchWinSec * fps));
nwin = min(nwin, numel(x));
nover = round(nwin * P.WelchOverlapFrac);
nfft = max(256, 2^nextpow2(nwin));
[Pxx, f] = pwelch(x, hann(nwin), nover, nfft, fps);
mask = f >= P.BandHz(1) & f <= P.BandHz(2);
fsub = f(mask);
psub = Pxx(mask);
if isempty(psub)
    peakPower = NaN; peakFreq = NaN; qFactor = NaN;
else
    [peakPower, ix] = max(psub);
    peakFreq = fsub(ix);

    halfH = peakPower/2;
    above = psub >= halfH;
    if any(above)
        bw = max(fsub(above)) - min(fsub(above));
        if bw <= 0
            qFactor = NaN;
        else
            qFactor = peakFreq / bw;
        end
    else
        qFactor = NaN;
    end
end

% ACF peak prominence / peak count
maxLag = round(P.AcfMaxLagSec * fps);
ac = xcorr(x-mean(x), maxLag, 'coeff');
lagAxis = (-maxLag:maxLag)/fps;
lagMask = lagAxis >= P.AcfPeakLagRangeSec(1) & lagAxis <= P.AcfPeakLagRangeSec(2);
acw = ac(lagMask);
if isempty(acw) || all(~isfinite(acw))
    acfProm = NaN;
    acfPeakCount = 0;
else
    [pks,~,~,proms] = findpeaks(acw, 'MinPeakProminence', 0);
    if isempty(proms)
        acfProm = 0;
        acfPeakCount = 0;
    else
        acfProm = max(proms);
        acfPeakCount = numel(pks);
    end
end

% Sliding-window freq stability
winF = max(20, round(P.WinSec*fps));
stepF = max(1, round(P.StepSec*fps));
starts = 1:stepF:max(1, numel(x)-winF+1);
winFreqs = nan(numel(starts),1);
for w = 1:numel(starts)
    s = starts(w);
    e = min(numel(x), s+winF-1);
    xi = x(s:e);
    [Pw, fw] = pwelch(xi, hann(min(nwin,numel(xi))), round(min(nwin,numel(xi))*P.WelchOverlapFrac), max(256,2^nextpow2(min(nwin,numel(xi)))), fps);
    m = fw >= P.BandHz(1) & fw <= P.BandHz(2);
    if any(m)
        [~,ix2] = max(Pw(m));
        f2 = fw(m);
        winFreqs(w) = f2(ix2);
    end
end
if nnz(isfinite(winFreqs)) >= 2 && mean(winFreqs,'omitnan')>0
    freqCVwindows = std(winFreqs,'omitnan') / mean(winFreqs,'omitnan');
else
    freqCVwindows = NaN;
end

% Instantaneous frequency regularity
xb = bandpass_1d(x, fps, P.BandHz);
hb = hilbert(xb);
ph = unwrap(angle(hb));
instF = [NaN; diff(ph)] * fps / (2*pi);
instF = instF(isfinite(instF));
instF = instF(instF >= P.BandHz(1)/4 & instF <= P.BandHz(2)*4);
if numel(instF)>=5 && mean(instF,'omitnan')>0
    instFreqCV = std(instF,'omitnan') / mean(instF,'omitnan');
else
    instFreqCV = NaN;
end

M = struct();
M.peakFreq = peakFreq;
M.peakPower = peakPower;
M.qFactor = qFactor;
M.acfProm = acfProm;
M.acfPeakCount = acfPeakCount;
M.freqCVwindows = freqCVwindows;
M.instFreqCV = instFreqCV;
M.f = f;
M.Pxx = Pxx;
M.winFreqs = winFreqs;
M.bandpassed = xb;
M.phase = ph;
end

%% =====================================================================
function xS = make_phase_randomized(x)
x = x(:);
N = numel(x);
X = fft(x);
amp = abs(X);
phi = angle(X);

% randomize positive frequencies, keep conjugate symmetry
Xnew = X;
if mod(N,2)==0
    pos = 2:(N/2);
    neg = N:-1:(N/2+2);
else
    pos = 2:((N+1)/2);
    neg = N:-1:((N+3)/2);
end
randPh = 2*pi*rand(numel(pos),1) - pi;
Xnew(pos) = amp(pos) .* exp(1i*randPh);
Xnew(neg) = amp(neg) .* exp(-1i*flipud(randPh));
Xnew(1) = X(1);
if mod(N,2)==0
    Xnew(N/2+1) = X(N/2+1);
end
xS = real(ifft(Xnew));
end

%% =====================================================================
function xS = make_ar1_surrogate(x)
x = x(:);
x = x - mean(x);
if numel(x) < 3
    xS = x;
    return;
end
phi = sum(x(1:end-1).*x(2:end)) / sum(x(1:end-1).^2);
phi = max(min(phi,0.995),-0.995);
res = x(2:end) - phi*x(1:end-1);
sig = std(res,'omitnan');
if ~isfinite(sig) || sig==0, sig = std(x,'omitnan'); end
xS = zeros(size(x));
xS(1) = x(1);
for t = 2:numel(x)
    xS(t) = phi*xS(t-1) + sig*randn;
end
end

%% =====================================================================
function Z = compute_zscores(realM, S)
mu = mean(S,1,'omitnan');
sd = std(S,0,1,'omitnan');
sd(sd==0 | ~isfinite(sd)) = NaN;
Z = struct();
Z.peakPower = (realM.peakPower - mu(1)) / sd(1);
Z.qFactor = (realM.qFactor - mu(2)) / sd(2);
Z.acfProm = (realM.acfProm - mu(3)) / sd(3);
Z.acfPeakCount = (realM.acfPeakCount - mu(4)) / sd(4);
Z.freqCVwindows = (realM.freqCVwindows - mu(5)) / sd(5);
Z.instFreqCV = (realM.instFreqCV - mu(6)) / sd(6);
% flipped signs so higher is better for stability / regularity
Z.freqCVwindows_neg = -(realM.freqCVwindows - mu(5)) / sd(5);
Z.instFreqCV_neg = -(realM.instFreqCV - mu(6)) / sd(6);
end

%% =====================================================================
function fig = make_representative_figure(x, realM, Sphase, Sar1, fps, P, classification)
t = (0:numel(x)-1)/fps;
fig = figure('Color','w','Position',[80 80 1600 950], 'Name','Representative oscillator-vs-noise test');
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tlo); hold(ax1,'on');
plot(ax1, t, x, 'k');
plot(ax1, t, realM.bandpassed, 'r');
xlabel(ax1,'Time (s)'); ylabel(ax1,'Signal');
title(ax1, sprintf('Representative signal | %s', classification));
legend(ax1, {'Raw rep signal','Bandpassed'}, 'Location','best');

ax2 = nexttile(tlo); hold(ax2,'on');
plot(ax2, realM.f, realM.Pxx, 'k', 'LineWidth', 1.5);
xline(ax2, P.BandHz(1), '--'); xline(ax2, P.BandHz(2), '--');
if isfinite(realM.peakFreq)
    plot(ax2, realM.peakFreq, realM.peakPower, 'ro', 'MarkerFaceColor', 'r');
end
xlim(ax2, [0 max(0.15, P.BandHz(2)*1.5)]);
xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'Power');
title(ax2, sprintf('PSD | peak=%.3f Hz | Q=%.2f', realM.peakFreq, realM.qFactor));

ax3 = nexttile(tlo); hold(ax3,'on');
plot(ax3, t, realM.phase, 'b');
xlabel(ax3,'Time (s)'); ylabel(ax3,'Unwrapped phase (rad)');
title(ax3, sprintf('Phase evolution | inst-f CV=%.2f', realM.instFreqCV));

ax4 = nexttile(tlo); hold(ax4,'on');
histogram(ax4, Sphase(:,1), 25, 'FaceAlpha', 0.5, 'DisplayName','Phase-rand');
histogram(ax4, Sar1(:,1),   25, 'FaceAlpha', 0.5, 'DisplayName','AR(1)');
xline(ax4, realM.peakPower, 'k', 'LineWidth', 2, 'DisplayName','Real');
xlabel(ax4,'Peak power'); ylabel(ax4,'Count'); title(ax4,'Real vs surrogate peak power'); legend(ax4,'Location','best');

ax5 = nexttile(tlo); hold(ax5,'on');
histogram(ax5, Sphase(:,2), 25, 'FaceAlpha', 0.5, 'DisplayName','Phase-rand');
histogram(ax5, Sar1(:,2),   25, 'FaceAlpha', 0.5, 'DisplayName','AR(1)');
xline(ax5, realM.qFactor, 'k', 'LineWidth', 2, 'DisplayName','Real');
xlabel(ax5,'Q-like factor'); ylabel(ax5,'Count'); title(ax5,'Real vs surrogate peak sharpness'); legend(ax5,'Location','best');

ax6 = nexttile(tlo); hold(ax6,'on');
histogram(ax6, Sphase(:,5), 25, 'FaceAlpha', 0.5, 'DisplayName','Phase-rand');
histogram(ax6, Sar1(:,5),   25, 'FaceAlpha', 0.5, 'DisplayName','AR(1)');
xline(ax6, realM.freqCVwindows, 'k', 'LineWidth', 2, 'DisplayName','Real');
xlabel(ax6,'Sliding-window freq CV'); ylabel(ax6,'Count'); title(ax6,'Lower is better: frequency stability'); legend(ax6,'Location','best');
end

%% =====================================================================
function groupTable = build_group_summary(T)
g = struct();
g.nRecordings = height(T);

g.meanOscillatorEvidenceScore = mean(T.oscillatorEvidenceScore, 'omitnan');
g.fracClass_oscillatorLike = mean(T.classification == "Oscillator-like");
g.fracClass_quasi = mean(T.classification == "Quasi-oscillatory / structured rhythmic dynamics");
g.fracClass_structuredNoise = mean(T.classification == "Structured slow noise / weak rhythmic structure");

g.fracPass_peakPower_phaseRand = mean(T.pass_peakPower_phaseRand);
g.fracPass_qFactor_phaseRand = mean(T.pass_qFactor_phaseRand);
g.fracPass_acfProm_phaseRand = mean(T.pass_acfProm_phaseRand);
g.fracPass_freqStability_phaseRand = mean(T.pass_freqStability_phaseRand);
g.fracPass_instFreqReg_phaseRand = mean(T.pass_instFreqReg_phaseRand);

g.fracPass_peakPower_ar1 = mean(T.pass_peakPower_ar1);
g.fracPass_qFactor_ar1 = mean(T.pass_qFactor_ar1);
g.fracPass_acfProm_ar1 = mean(T.pass_acfProm_ar1);
g.fracPass_freqStability_ar1 = mean(T.pass_freqStability_ar1);
g.fracPass_instFreqReg_ar1 = mean(T.pass_instFreqReg_ar1);

g.meanReal_peakFreq = mean(T.real_peakFreq, 'omitnan');
g.meanReal_qFactor = mean(T.real_qFactor, 'omitnan');
g.meanReal_freqCVwindows = mean(T.real_freqCVwindows, 'omitnan');
g.meanReal_instFreqCV = mean(T.real_instFreqCV, 'omitnan');

g.note = "Per-fish first. Use fractions of surrogate exceedance as publication summary.";
groupTable = struct2table(g);
end

%% =====================================================================
function figs = make_group_figures(T, G, figDir, saveFigures)
figs = struct();

% ---------- Fig 2: per-recording effect sizes ----------
figs.effectSizes = figure('Color','w','Position',[80 80 1500 900], 'Name','Per-recording surrogate effect sizes');
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile(tlo); boxplot([T.z_peakPower_phaseRand, T.z_peakPower_ar1], 'Labels', {'PeakPow vs PR','PeakPow vs AR1'}); ylabel('z-like effect'); title('Peak power');
nexttile(tlo); boxplot([T.z_qFactor_phaseRand, T.z_qFactor_ar1], 'Labels', {'Q vs PR','Q vs AR1'}); ylabel('z-like effect'); title('Peak sharpness');
nexttile(tlo); boxplot([T.z_acfProm_phaseRand, T.z_acfProm_ar1], 'Labels', {'ACF prom vs PR','ACF prom vs AR1'}); ylabel('z-like effect'); title('ACF prominence');
nexttile(tlo); boxplot([T.z_freqCVwindows_phaseRand, T.z_freqCVwindows_ar1], 'Labels', {'Freq stab vs PR','Freq stab vs AR1'}); ylabel('higher = more stable than surrogate'); title('Windowed stability');
nexttile(tlo); boxplot([T.z_instFreqCV_phaseRand, T.z_instFreqCV_ar1], 'Labels', {'Inst.reg vs PR','Inst.reg vs AR1'}); ylabel('higher = more regular than surrogate'); title('Instantaneous frequency');
nexttile(tlo); histogram(T.oscillatorEvidenceScore, 12); xlabel('Oscillator evidence score'); ylabel('# recordings'); title('Overall evidence score');

% ---------- Fig 3: fraction passing ----------
figs.passFractions = figure('Color','w','Position',[100 100 1400 700], 'Name','Fraction of recordings exceeding surrogate thresholds');
labels = {'Peak power','Q factor','ACF prom','Freq stability','Inst-f regularity'};
PR = [mean(T.pass_peakPower_phaseRand), mean(T.pass_qFactor_phaseRand), mean(T.pass_acfProm_phaseRand), mean(T.pass_freqStability_phaseRand), mean(T.pass_instFreqReg_phaseRand)];
AR = [mean(T.pass_peakPower_ar1), mean(T.pass_qFactor_ar1), mean(T.pass_acfProm_ar1), mean(T.pass_freqStability_ar1), mean(T.pass_instFreqReg_ar1)];
bar([PR(:), AR(:)]);
set(gca,'XTickLabel',labels,'XTickLabelRotation',25);
ylabel('Fraction of recordings passing');
legend({'vs phase-rand','vs AR(1)'}, 'Location','best');
title('How often does the real data beat surrogate nulls?');

% ---------- Fig 4: real metric distributions ----------
figs.realMetrics = figure('Color','w','Position',[120 120 1500 850], 'Name','Real metric distributions');
tlo2 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile(tlo2); histogram(T.real_peakFreq, 12); xlabel('Peak frequency (Hz)'); ylabel('# recordings'); title('Preferred frequencies');
nexttile(tlo2); histogram(T.real_qFactor, 12); xlabel('Q-like factor'); ylabel('# recordings'); title('Peak sharpness');
nexttile(tlo2); histogram(T.real_acfProm, 12); xlabel('ACF peak prominence'); ylabel('# recordings'); title('Temporal periodicity');
nexttile(tlo2); histogram(T.real_freqCVwindows, 12); xlabel('Windowed freq CV'); ylabel('# recordings'); title('Lower = more stable');
nexttile(tlo2); histogram(T.real_instFreqCV, 12); xlabel('Instantaneous freq CV'); ylabel('# recordings'); title('Lower = more regular');
nexttile(tlo2);
counts = [sum(T.classification=="Oscillator-like"), sum(T.classification=="Quasi-oscillatory / structured rhythmic dynamics"), sum(T.classification=="Structured slow noise / weak rhythmic structure")];
bar(counts); set(gca,'XTickLabel',{'Oscillator-like','Quasi-osc','Structured noise'}, 'XTickLabelRotation',20); ylabel('# recordings'); title('Classification');

if saveFigures
    exportgraphics(figs.effectSizes, fullfile(figDir, 'group_effect_sizes.png'), 'Resolution', 200);
    exportgraphics(figs.passFractions, fullfile(figDir, 'group_pass_fractions.png'), 'Resolution', 200);
    exportgraphics(figs.realMetrics, fullfile(figDir, 'group_real_metric_distributions.png'), 'Resolution', 200);
end
end

%% =====================================================================
function [dff, fps] = extract_dff_and_fps(R, dffNames, fpsNames, filePath)
dff = [];
fps = [];
vars = fieldnames(R);

for i = 1:numel(dffNames)
    vn = dffNames{i};
    if isfield(R, vn)
        cand = R.(vn);
        if isnumeric(cand) && ismatrix(cand) && size(cand,1) > 1 && size(cand,2) > 10
            dff = cand;
            break;
        end
    end
end

if isempty(dff)
    nested = {'dff_new','dff','DFF','DV_DFFmovwindow'};
    for i = 1:numel(vars)
        x = R.(vars{i});
        if isstruct(x)
            for j = 1:numel(nested)
                if isfield(x, nested{j})
                    cand = x.(nested{j});
                    if isnumeric(cand) && ismatrix(cand) && size(cand,1) > 1 && size(cand,2) > 10
                        dff = cand;
                        break;
                    end
                end
            end
        end
        if ~isempty(dff), break; end
    end
end

if isempty(dff)
    error('Could not find dF/F matrix in %s', filePath);
end

for i = 1:numel(fpsNames)
    vn = fpsNames{i};
    if isfield(R, vn) && isnumeric(R.(vn)) && isscalar(R.(vn))
        fps = double(R.(vn));
        break;
    end
end

if isempty(fps)
    for i = 1:numel(vars)
        x = R.(vars{i});
        if isstruct(x)
            for j = 1:numel(fpsNames)
                vn = fpsNames{j};
                if isfield(x, vn) && isnumeric(x.(vn)) && isscalar(x.(vn))
                    fps = double(x.(vn));
                    break;
                end
            end
        end
        if ~isempty(fps), break; end
    end
end

if isempty(fps)
    error('Could not find fps in %s', filePath);
end
end

%% =====================================================================
function Xz = zscore_fill_rows(X)
mu = mean(X, 2, 'omitnan');
sd = std(X, 0, 2, 'omitnan');
sd(sd == 0 | ~isfinite(sd)) = 1;
Xz = (X - mu) ./ sd;
Xz(~isfinite(Xz)) = 0;
end

%% =====================================================================
function xbp = bandpass_1d(x, fps, bandHz)
Wn = bandHz / (fps/2);
if Wn(2) >= 1
    error('Upper band edge must be below Nyquist.');
end
[b,a] = butter(3, Wn, 'bandpass');
xbp = filtfilt(b, a, x(:));
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
