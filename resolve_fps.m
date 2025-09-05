function fps = resolve_fps(resultsMat, metaMat)
% Try in this order:
% 1) fps inside Suite2p results (dffs_repact_respcells.mat)
% 2) clear fields inside metadata_multimodal.mat
% 3) any .ini in the folder (scan text)
% 4) fallback default

    fps = NaN;

    % 1) Suite2p results
    try
        R = load(resultsMat, 'fps');
        if isfield(R,'fps') && isfinite(R.fps) && R.fps > 0
            fps = double(R.fps);
            fprintf('  fps from results: %.4f\n', fps);
            return;
        end
    catch
    end

    % 2) metadata_multimodal.mat
    try
        M = load(metaMat);
        cand = {'fps','frameRate','frame_rate','imaging_fps','scanimage_fps','acq_fps','framerate','FrameRate'};
        for k = 1:numel(cand)
            if isfield(M, cand{k})
                v = double(M.(cand{k}));
                if isscalar(v) && isfinite(v) && v > 0
                    fps = v;
                    fprintf('  fps from metadata (%s): %.4f\n', cand{k}, fps);
                    return;
                end
            end
        end
    catch
    end

    % 3) parse any .ini in the same folder
    try
        folder = fileparts(resultsMat);
        iniList = dir(fullfile(folder, '*.ini'));
        if ~isempty(iniList)
            fps = parse_ini_for_fps(fullfile(iniList(1).folder, iniList(1).name));
            if isfinite(fps) && fps > 0
                fprintf('  fps from ini: %.4f\n', fps);
                return;
            end
        end
    catch
    end

    % 4) fallback default
    fps = 14.64;
    fprintf('  fps fallback default: %.4f\n', fps);
end

function fps = parse_ini_for_fps(iniPath)
% Very permissive: looks for numbers after tokens like fps, frameRate, etc.
    fps = NaN;
    try
        txt = fileread(iniPath);
        % examples it matches: fps=12.5, frameRate : 12.5, imaging_fps 12.5
        pat = '(?i)\b(fps|frame[_ ]?rate|imaging[_ ]?fps|scanimage[_ ]?fps)\b[^0-9\-\.]*([0-9]+(?:\.[0-9]+)?)';
        tok = regexp(txt, pat, 'tokens', 'once');
        if ~isempty(tok)
            fps = str2double(tok{2});
        end
    catch
    end
end
