
function [t, y, fs] = load_etho_var_robust(fpath, field)
% Load a variable and return a matching time vector t and sampling rate fs.
    S = load(fpath);
    assert(isfield(S, field), 'Field %s not found in %s', field, fpath);

    y = S.(field)(:);
    n = numel(y);

    % try matching-length time vectors first
    cand = {'fra_tim_uni_frames','tim_uni_frames','fra_tim_uni','fra_tim','time','t'};
    t = [];
    for k = 1:numel(cand)
        nm = cand{k};
        if isfield(S, nm) && numel(S.(nm)) == n
            t = S.(nm)(:);
            break;
        end
    end

    fs = [];
    if isfield(S,'fra_rat_get') && ~isempty(S.fra_rat_get), fs = double(S.fra_rat_get); end
    if isempty(fs) && isfield(S,'fra_rat') && ~isempty(S.fra_rat), fs = double(S.fra_rat); end

    if isempty(t)
        % build from fps
        assert(~isempty(fs) && fs>0, 'Need time or fps in %s', fpath);
        t = (0:n-1)'/fs;
    else
        dt = median(diff(t));
        if isempty(fs) || ~isfinite(fs) || fs<=0, fs = 1/dt; end
    end
end