function r_ts = bout_rate_ts(t, bouts, win_sec)
% Bouts/min as a continuous time series (boxcar kernel of length win_sec)
    r_ts = zeros(size(t));
    if isempty(bouts), return; end
    starts = t(bouts(:,1));
    win = win_sec;
    % naive O(T * nB) OK for moderate sizes
    for i=1:numel(t)
        r_ts(i) = 60 * sum(starts >= (t(i)-win) & starts <= t(i)) / max(win,eps);
    end
end