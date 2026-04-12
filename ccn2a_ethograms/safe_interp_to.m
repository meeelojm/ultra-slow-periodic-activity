function yq = safe_interp_to(t, y, t_q)
% Deduplicate t and interpolate y -> t_q
    [tu, ia] = unique(t(:), 'stable');
    yu = y(ia);
    yq = interp1(