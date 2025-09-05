% at end of the script:
function resFile = findRes(fid, validFolders)
    resFile = '';
    for ii = 1:numel(validFolders)
        [~, nm] = fileparts(validFolders{ii});
        if strcmpi(nm, fid)
            D = dir(fullfile(validFolders{ii}, '*_Results_dff.mat'));
            if ~isempty(D)
                resFile = fullfile(D(1).folder, D(1).name);
                return
            end
        end
    end
end
