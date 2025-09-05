function resFile = findResultsForFish(fid, validFolders)
    resFile = '';
    for ii = 1:numel(validFolders)
        [~, nm] = fileparts(validFolders{ii});
        if strcmpi(nm, fid)
            D = dir(fullfile(validFolders{ii}, '*_Results_dff.mat'));
            if ~isempty(D)
                resFile = fullfile(validFolders{ii}, D(1).name);
                return
            end
        end
    end
end
