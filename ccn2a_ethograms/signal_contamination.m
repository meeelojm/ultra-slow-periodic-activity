%% Inputs
fs    = 120;                  % sampling rate (Hz)
fband = [0 30];               % band of interest (Hz)
sig_row = horzcat(hea_sig,ope_sig,mou_sig,eye_sig);
X     = sig_row;              % T x 4 (heart, operculum, mouth, eye)

chanNames = {'Heart','Operculum','Mouth','Eye'};

%% 1) Clean & detrend
Xc = X;
for k = 1:size(X,2)
    x = X(:,k);
    x(~isfinite(x)) = 0;
    % Detrend (robust enough for long traces)
    Xc(:,k) = detrend(double(x));
end

fs = 120;                               % Hz
t  = (0:size(Xc,1)-1)'/fs;              % seconds
labels = ["Heart","Operculum","Mouth","Eye"];

tt = array2timetable(Xc, ...
    'RowTimes', seconds(t), ...
    'VariableNames', cellstr(labels));
figure
stackedplot(tt);
title('Ethogram channels'); xlabel('Time (s)');

%% 2) Autocorrelation per channel
fs = 120;                      % Hz (already defined)
labels = ["Heart","Operculum","Mouth","Eye"];

maxLagSec = 10;                % show up to 10 s lag (tweak as you like)
maxLags   = round(maxLagSec*fs);

N = size(Xc,1);
acfS   = [];
tauSec = (0:maxLags)'/fs;      % nonnegative lags in seconds
dom    = struct('lag_s',nan,'freq_hz',nan);

% 95% white-noise CI (helps judge significance)
ci = 1.96/sqrt(N);

for k = 1:size(Xc,2)
    % normalized autocorrelation up to maxLags
    [acf,lags] = xcorr(Xc(:,k), maxLags, 'coeff');
    acf = acf(lags>=0);              % keep nonnegative lags
    acfS = [acfS,acf];

    % simple dominant period estimate = first peak after 0
    [pks,locs] = findpeaks(acf(2:end), tauSec(2:end), ...
                           'MinPeakProminence', 0.02);
    if ~isempty(locs)
        dom(k).lag_s   = locs(1);
        dom(k).freq_hz = 1/locs(1);
    end
end

tt = array2timetable(Xc, ...
    'RowTimes', seconds(t), ...
    'VariableNames', cellstr(labels));
figure
stackedplot(acfS);
title('Ethogram ACS'); xlabel('Time (s)');

%% 2) PSD per channel & dominant frequency in [8 15] Hz
nfftPSD = 1024;               % gives ~0.0293 Hz resolution at 120 Hz
domf = nan(4,1);
bw   = nan(4,1);              % occupied bandwidth (99% default)

for k = 1:4
    x = acfS(:,k);
    [Pxx,F] = pwelch(x, nfftPSD, [], nfftPSD, fs);  % column vector
    in = (F >= fband(1)) & (F <= fband(2));
    if any(in)
        [~,ix] = max(Pxx(in));
        fseg   = F(in);
        pseg   = Pxx(in);
        domf(k) = fseg(ix);              % dominant freq (Hz)
        bw(k)   = obw(pseg, fseg);       % occupied bandwidth (Hz)
    end
end

fprintf('\nDominant frequencies (8–15 Hz band):\n');
for k = 1:4
    fprintf('  %s: domf = %.4f Hz, bandwidth ≈ %.2f Hz\n', chanNames{k}, domf(k), bw(k));
end

figure
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
for k = 1:4
    nexttile; hold on; grid on;
    [Pxx,F] = pwelch(acfS(:,k), nfftPSD, [], nfftPSD, fs);
    plot(F, 10*log10(Pxx), 'LineWidth', 1);
    % band markers
    xline(fband(1),'r:'); xline(fband(2),'r:');
    % domf marker (if available)
    if isfinite(domf(k))
        xline(domf(k),'k-','LineWidth',1.2);
        txt = sprintf('domf=%.2f Hz, BW≈%.2f', domf(k), bw(k));
    else
        txt = 'no peak in band';
    end
    title(sprintf('Ch %d (%s): %s', k, chanNames{k}, txt));
    xlim([0 30]); ylabel('PSD (dB)'); 
end
xlabel('Frequency (Hz)');
sgtitle('PSD per channel (8–15 Hz band shown)');

figure
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
for k = 1:4
    nexttile; hold on; grid on;
    [Pxx,F] = pwelch(Xc(:,k), nfftPSD, [], nfftPSD, fs);
    plot(F, 10*log10(Pxx), 'LineWidth', 1);
    % band markers
    xline(fband(1),'r:'); xline(fband(2),'r:');
    % domf marker (if available)
    if isfinite(domf(k))
        xline(domf(k),'k-','LineWidth',1.2);
        txt = sprintf('domf=%.2f Hz, BW≈%.2f', domf(k), bw(k));
    else
        txt = 'no peak in band';
    end
    title(sprintf('Ch %d (%s): %s', k, chanNames{k}, txt));
    xlim([0 30]); ylabel('PSD (dB)'); 
end
xlabel('Frequency (Hz)');
sgtitle('PSD per channel (8–15 Hz band shown)');

%% 3) Make the 3-panel figure
figure('Color','w','Position',[100 100 1500 420]);

% ---------- Panel A: long raw view (downsampled for speed) ----------
subplot(1,3,1); hold on;
T = size(Xc,1);
dec = max(1, floor(T/5e5));             % decimate just for plotting density
t  = (0:T-1)'/fs;
for k = 1:4
    plot(t(1:dec:end), Xc(1:dec:end,k) + 20*(k-1), 'b');  % vertical offsets
end
yticks(20*(0:3));
yticklabels(chanNames);
xlabel('Time (s)');
title('Raw signals (offset)');

% ---------- Panel B: short filtered view (8–15 Hz) ----------
subplot(1,3,2); hold on;
d = designfilt('bandpassiir','FilterOrder',4, ...
               'HalfPowerFrequency1',fband(1), ...
               'HalfPowerFrequency2',fband(2), ...
               'SampleRate',fs);
winSec = 15;                         % show first ~15 s
Nshow  = min(T, round(fs*winSec));
ts     = t(1:Nshow);
for k = 1:4
    xf = filtfilt(d, Xc(:,k));
    % normalize for display so amplitudes are comparable panel-wise
    xx = xf(1:Nshow);
    xx = (xx - min(xx));  xx = xx / max(1e-12, max(xx));
    plot(ts, xx + (k-1), 'LineWidth', 1);  % stack at 0..3
end
yticks(0:3); yticklabels(chanNames);
xlabel('Time (s)'); ylim([-0.2 3.2]);
title(sprintf('Bandpass %.0f–%.0f Hz (normalized view)', fband(1), fband(2)));

% ---------- Panel C: PSD with band and domf ----------
subplot(1,3,3); 
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
for k = 1:4
    nexttile; hold on; grid on;
    [Pxx,F] = pwelch(Xc(:,k), nfftPSD, [], nfftPSD, fs);
    plot(F, 10*log10(Pxx), 'LineWidth', 1);
    % band markers
    xline(fband(1),'r:'); xline(fband(2),'r:');
    % domf marker (if available)
    if isfinite(domf(k))
        xline(domf(k),'k-','LineWidth',1.2);
        txt = sprintf('domf=%.2f Hz, BW≈%.2f', domf(k), bw(k));
    else
        txt = 'no peak in band';
    end
    title(sprintf('Ch %d (%s): %s', k, chanNames{k}, txt));
    xlim([0 30]); ylabel('PSD (dB)'); 
end
xlabel('Frequency (Hz)');
sgtitle('PSD per channel (8–15 Hz band shown)');

