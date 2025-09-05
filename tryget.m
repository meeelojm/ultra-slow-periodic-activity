%% --------- helper: nested field getter ---------
function x = tryget(S, pathcells)
% pathcells is a cellstr list of nested field names, e.g., {'results','dff'}
    x = [];
    try
        switch numel(pathcells)
            case 1
                if isfield(S, pathcells{1}), x = S.(pathcells{1}); end
            case 2
                if isfield(S, pathcells{1}) && isfield(S.(pathcells{1}), pathcells{2})
                    x = S.(pathcells{1}).(pathcells{2});
                end
            otherwise
                % extend if you need deeper nesting
        end
    catch
        x = [];
    end
end