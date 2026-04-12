function [t, y] = pick_series(S, uni_field, sig_field)
% Return time+series preferring *_uni_frames if present; else *_sig with its time.
    if isfield(S,uni_field) && ~isempty(S.(uni_field))
        y = S.(uni_field)(:);
        t = pick_time_for_field(S, numel(y));
    elseif isfield(S,sig_field) && ~isempty(S.(sig_field))
        y = S.(sig_field)(:);
        t = pick_time_for_field(S, numel(y));
    else
        error('Neither %s nor %s present.', uni_field, sig_field);
    end
end