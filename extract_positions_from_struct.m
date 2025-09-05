%% -------- helper: extract positions from a results struct (works for structs or objects) --------
function pos = extract_positions_from_struct(S)
    pos = [];
    try
        % 1) Try exact paths first (handles nested S.S.results... too)
        paths = {
            {'results','positions_nodoubles'}
            {'S','results','positions_nodoubles'}
            {'results','positions'}
            {'S','results','positions'}
            {'results','DV_Centroids'}
            {'results','centroids'}
            {'positions_nodoubles'}
            {'positions'}
        };
        for p = 1:numel(paths)
            [ok,val] = try_path(S, paths{p});
            if ok && ~isempty(val)
                pos = val; break;
            end
        end

        % 2) Fuzzy search under results if still empty
        if isempty(pos)
            [ok,R] = try_path(S, {'results'});
            if ok && ~isempty(R)
                fns = fieldnames(struct(R));   % struct() handles many objects; safe if already struct
                low = lower(fns);
                % prefer fields containing both 'position' and some 'nodouble' spelling
                idx = find(contains(low,'position') & ...
                           (contains(low,'nodouble') | contains(low,'no_double') | ...
                            contains(low,'nodoubles') | contains(low,'nodbl')));
                if isempty(idx), idx = find(contains(low,'position') | contains(low,'centroid')); end
                if ~isempty(idx), pos = R.(fns{idx(1)}); end
            end
        end

        % 3) Top-level fuzzy fallback
        if isempty(pos)
            fns = fieldnames(S);
            low = lower(fns);
            idx = find(contains(low,'positions_nodoubles') | contains(low,'positions') | contains(low,'centroids'),1);
            if ~isempty(idx), pos = S.(fns{idx}); end
        end

        % 4) Shape & type
        if isempty(pos)
            warning('extract_positions_from_struct: could not find positions. Top-level fields: %s', strjoin(fieldnames(S),', '));
            return;
        end
        pos = double(pos);
        if size(pos,2) < 2, error('Positions must have at least 2 columns (X,Y).'); end
        if size(pos,2) == 2, pos(:,3) = 0; end
        if size(pos,2) > 3, pos = pos(:,1:3); end

    catch ME
        warning('extract_positions_from_struct: %s', ME.message);
        pos = [];
    end
end

% --- helper: safely follow a nested field/property path using subsref (works for objects too)
function [ok,val] = try_path(obj, pathCells)
    val = obj; ok = true;
    try
        for i = 1:numel(pathCells)
            s = struct('type','.', 'subs', pathCells{i});
            val = builtin('subsref', val, s);
        end
    catch
        ok = false; val = [];
    end
end
