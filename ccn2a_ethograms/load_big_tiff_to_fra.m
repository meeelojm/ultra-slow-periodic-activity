function [fra, info] = load_big_tiff_to_fra(tifFile, frameRange)
%LOAD_BIG_TIFF_TO_FRA  Read (Big)TIFF stack into a 3-D array 'fra'.
%   [FRA, INFO] = LOAD_BIG_TIFF_TO_FRA(TIFFILE)
%   [FRA, INFO] = LOAD_BIG_TIFF_TO_FRA(TIFFILE, FRAME RANGE)  % e.g. 1:1000

    % --- arguments with safer validation ---
    arguments
        tifFile (1,1) string
        frameRange = []   % allow empty by default
    end

    info = imfinfo(tifFile);
    T = numel(info);

    % If not given, take all frames
    if isempty(frameRange)
        frameRange = 1:T;
    end

    frameRange = frameRange(:)';  % ensure row vector
    Tsel = numel(frameRange);

    h  = info(1).Height;
    w  = info(1).Width;
    bs = info(1).BitDepth;  % bits per pixel

    switch bs
        case 8,  cls = 'uint8';
        case 16, cls = 'uint16';
        otherwise
            error('Unsupported BitDepth=%d. Only 8/16 bpp supported.', bs);
    end

    % Memory estimate
    bytesPerPixel = bs/8;
    needGB = h*w*Tsel*bytesPerPixel/1024^3;
    if needGB > 16
        warning('This read will allocate ~%.1f GB. Ensure you have enough RAM.', needGB);
    end

    fra = zeros(h, w, Tsel, cls);

    % Waitbar
    msg = sprintf('Loading %d frames...', Tsel);
    hWait = waitbar(0, msg, 'Name','TIFF Import', ...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
    setappdata(hWait,'canceling',false);
    updateEvery = max(1, round(Tsel/100));

    t = Tiff(tifFile, 'r');
    clean = onCleanup(@() (isvalid(hWait) && delete(hWait)));

    try
        for i = 1:Tsel
            if getappdata(hWait,'canceling'), error('User canceled.'); end
            k = frameRange(i);
            t.setDirectory(k);
            fra(:,:,i) = t.read();
            if mod(i, updateEvery) == 0 || i == Tsel
                waitbar(i/Tsel, hWait, sprintf('Loading frame %d / %d', i, Tsel));
                drawnow limitrate
            end
        end
    catch ME
        t.close();
        rethrow(ME);
    end
    t.close();
end
