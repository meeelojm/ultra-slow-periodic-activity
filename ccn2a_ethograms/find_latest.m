function [fullpath, tried] = find_latest(baseRoots, fish, condAliases, tail, debug)
%FIND_LATEST Find the newest file matching a pattern across possible roots/aliases.
%
% Inputs
% ------
% baseRoots    : cell array of root folders
% fish         : e.g. 'f010'
% condAliases  : e.g. {'highact','high'}
% tail         : e.g. fullfile('suite2p','ethogram_*.mat')
% debug        : true/false
%
% Outputs
% -------
% fullpath     : full path to newest matching file
% tried        : string array of tried patterns

    if nargin < 5 || isempty(debug)
        debug = false;
    end

    if isstring(fish)
        fish = char(fish);
    end

    if ischar(condAliases) || isstring(condAliases)
        condAliases = {char(condAliases)};
    end

    tried = strings(0,1);
    fullpath = '';

    for r = 1:numel(baseRoots)
        root = baseRoots{r};

        if ~exist(root, 'dir')
            if debug
                fprintf('Skipping missing root: %s\n', root);
            end
            continue
        end

        for a = 1:numel(condAliases)
            sub = condAliases{a};
            patt = fullfile(root, fish, sub, tail);
            tried(end+1,1) = string(patt); %#ok<AGROW>

            d = dir(patt);

            if debug
                fprintf('Scanning: %s -> %d match(es)\n', patt, numel(d));
            end

            if ~isempty(d)
                [~, ix] = max([d.datenum]);
                fullpath = fullfile(d(ix).folder, d(ix).name);

                if debug
                    fprintf('Found: %s\n', fullpath);
                end
                return
            end
        end
    end

    if debug
        fprintf('No file matched any of these patterns:\n');
        fprintf('  %s\n', tried);
    end

    error('File not found for fish "%s" and pattern "%s".', fish, tail);
end