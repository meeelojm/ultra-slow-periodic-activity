function t = pick_time_for_field(S, ylen)
% Return a time vector that matches this series length or build from fps
    cand = {'fra_tim_uni_frames','tim_uni_frames','fra_tim_uni','fra_tim','time','t'};
    t = [];
    for k = 1:numel(cand)
        nm = cand{k};
        if isfield(S, nm) && numel(S.(nm)) == ylen
            t = S.(nm)(:);
            break;
        end
    end
    if isempty(t)
        fs = [];
        if isfield(S,'fra_rat_get') && ~isempty(S.fra_rat_get), fs = double(S.fra_rat_get); end
        if isempty(fs) && isfield(S,'fra_rat') && ~isempty(S.fra_rat), fs = double(S.fra_rat); end
        assert(~isempty(fs) && fs>0, 'No matching time vector or fps found.');
        t = (0:ylen-1)'/fs;
    end
end