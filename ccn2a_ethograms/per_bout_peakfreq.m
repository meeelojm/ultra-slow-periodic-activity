function f_list = per_bout_peakfreq(t, y, fs, bouts, band)
% Welch PSD per bout; return list of peak freqs (Hz)
    f_list = nan(0,1);
    if isempty(bouts), return; end
    for i = 1:size(bouts,1)
        seg = y(bouts(i,1):bouts(i,2));
        N   = numel(seg);
        if N < round(0.3*fs), continue; end   % skip very short bouts
        nfft = 2^nextpow2(max(N,512));
        [Pxx,F] = pwelch(seg, hamming(min(N,512)), [], nfft, fs);
        m = (F>=band(1)) & (F<=band(2));
        if any(m)
            [~,ix] = max(Pxx(m));
            f = F(m); f_list(end+1,1) = f(ix); %#ok<AGROW>
        end
    end
end