function S = load_ethogram_file(baseRoots, fish, cond, debug)
%LOAD_ETHOGRAM_FILE Find and load the latest ethogram MAT file for one fish/condition.
%
% Example:
%   S = load_ethogram_file(baseRoots, 'f010', 'high', true);

    if nargin < 4 || isempty(debug)
        debug = false;
    end

    condAliases = {[cond 'act'], cond};
    tailPattern = fullfile('suite2p', 'ethogram_*.mat');

    [matPath, tried] = find_latest(baseRoots, fish, condAliases, tailPattern, debug);
    tmp = load(matPath);

    S = struct();
    S.path   = matPath;
    S.tried  = tried;
    S.fish   = fish;
    S.cond   = cond;
    S.data   = tmp;
end