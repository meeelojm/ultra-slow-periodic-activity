function out = domfreq_of_trace(y, fs, cfg)
% y: column vector (time); fs: Hz
% cfg: struct with fields num_lags, windowSize, maxpx, minPeakThreshold, freqRange

    sig_row = double(y(:)).';        % 1 × T (rows=signals, cols=time)
    start_frame = 1; 
    end_frame   = numel(sig_row);

    [acs, domf, Pxxs, Fss, isoinds] = calcsort_autocorr_freq_analysis_v3( ...
        sig_row(:, start_frame:end_frame), ...
        cfg.num_lags, fs, cfg.windowSize, cfg.maxpx, cfg.minPeakThreshold, cfg.freqRange);

    out = struct('domf',domf, 'acs',acs, 'Pxxs',Pxxs, 'Fss',Fss, 'isoinds',isoinds);
end