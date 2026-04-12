function [t, ang] = load_eth(fpath)
    % Load tail angle and time; repair mismatches robustly.
    S = load(fpath, 'tai_ang_uni_frames', 'fra_tim_uni', 'fra_rat');
    ang = S.tai_ang_uni_frames(:);                 % column vector
    t0  = []; fs = [];
    if isfield(S,'fra_tim_uni'), t0 = S.fra_tim_uni(:); end
    if isfield(S,'fra_rat'),     fs = S.fra_rat;        end

    if ~isempty(t0) && numel(t0)==numel(ang)
        t = t0;                                       % perfect: use as-is
    elseif ~isempty(t0) && abs(numel(t0)-numel(ang))<=5
        % small off-by-few: crop both to min length
        L = min(numel(t0), numel(ang));
        t = t0(1:L);
        ang = ang(1:L);
        fprintf('[%s] cropped to L=%d (t=%d, ang=%d)\n', fpath, L, numel(t0), numel(ang));
    elseif ~isempty(fs)
        % rebuild a uniform timebase to MATCH ang length
        t_start = ~isempty(t0) * t0(1);              % if t0 exists, keep its start
        t = t_start + (0:numel(ang)-1)'/fs;
        fprintf('[%s] rebuilt t from fra_rat to match ang (L=%d)\n', fpath, numel(ang));
    elseif ~isempty(t0)
        % no fs, but have t0: just crop to min (last resort)
        L = min(numel(t0), numel(ang));
        t = t0(1:L);
        ang = ang(1:L);
        warning('[%s] no fra_rat; cropped to L=%d', fpath, L);
    else
        error('[%s] No time info (need fra_tim_uni or fra_rat).', fpath);
    end

    % ensure strictly increasing time (unique X for interp1)
    [t, ia] = unique(t, 'stable'); 
    ang = ang(ia);
end
