function yq = resamp_to(t, y, t_common)
    % Deduplicate time just in case, then interp.
    [tu, ia] = unique(t, 'stable');
    yu = y(ia);
    yq = interp1(tu, yu, t_common, 'linear');  % no extrapolation since we use intersection
end
