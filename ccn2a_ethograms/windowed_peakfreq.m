function [t_win, f_win] = windowed_peakfreq(t, y, fs, win_sec, hop_sec, band)
% Sliding window Welch peak frequency; returns centers and peak Hz
    N = numel(t);
    win = max(1, round(win_sec*fs));
    hop = max(1, round(hop_sec*fs));
    idxs = 1:hop:(N-win+1);
    t_win = t(idxs + floor(win/2));
    f_win = nan(numel(idxs),1);
    nfft  = 2^nextpow2(max(win,1024));
    for k=1:numel(idxs)
        seg = y(idxs(k):idxs(k)+win-1);
        [Pxx,F] = pwelch(seg, hamming(min(numel(seg),512)), [], nfft, fs);
        m = (F>=band(1)) & (F<=band(2));
        if any(m)
            [~,ix] = max(Pxx(m));
            fw = F(m); f_win(k) = fw(ix);
        end
    end
end