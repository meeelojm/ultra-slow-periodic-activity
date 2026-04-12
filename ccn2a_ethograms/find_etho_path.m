function fpath = find_etho_path(base_dirs, fish, cond)
% Return full path to the newest ethogram_*.mat for a fish/condition.
    subs_try = { [cond 'act'], cond };   % try 'highact' then 'high'
    for b = 1:numel(base_dirs)
        root = base_dirs{b};
        if ~isfolder(root), continue; end
        for s = 1:numel(subs_try)
            pat = fullfile(root, fish, subs_try{s}, 'suite2p', 'ethogram_*.mat');
            D = dir(pat);
            if ~isempty(D)
                [~,ix] = max([D.datenum]);
                fpath = fullfile(D(ix).folder, D(ix).name);
                return;
            end
        end
    end
    error('No ethogram file found for %s / %s in provided base_dirs.', fish, cond);
end