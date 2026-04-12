function r = get_fish_root(baseRoots, fish)
% return the first root that contains fish folder
    for b = 1:numel(baseRoots)
        p = fullfile(baseRoots{b}, fish);
        if isfolder(p), r = fullfile(baseRoots{b}, fish); return; end
    end
    r = fullfile(baseRoots{1}, fish); % fallback
end
