function save_stack_big_tiff(fra, tifFile)
% SAVE_STACK_BIG_TIFF  Write H×W×T stack to BigTIFF using Tiff class.
%   - Converts to 8-bit if input isn't uint8/uint16
%   - Shows a waitbar with Cancel
%   - Robust cleanup of file handle and waitbars

    if nargin<2 || isempty(tifFile)
        tifFile = fullfile(pwd, sprintf('registered_fra_%s.tif', datestr(now,'yyyymmdd_HHMMSS')));
    end

    [h,w,T] = size(fra);

    % Cast/scale if needed
    if ~isa(fra,'uint8') && ~isa(fra,'uint16')
        fra = uint8(255 * mat2gray(fra));  % simple 8-bit scaling
    end
    if isa(fra,'uint16')
        bits = 16;
    else
        bits = 8;
    end

    % Kill any old waitbars
    delete(findall(0,'Type','figure','Tag','TMWWaitbar'));

    % Create waitbar with cancel button
    hWait = waitbar(0,'Saving frames to BigTIFF...','Name','Saving Progress', ...
        'CreateCancelBtn', @(src,evt) setappdata(hWait,'canceling',true));
    setappdata(hWait,'canceling',false);

    % Open BigTIFF
    t = Tiff(tifFile,'w8');  % 'w8' => BigTIFF
    % Robust cleanup (local helper below)
    c_closeTiff = onCleanup(@() close_tiff_quiet(t));
    c_killBars  = onCleanup(@() delete(findall(0,'Type','figure','Tag','TMWWaitbar')));

    % Common tags
    tagstruct.ImageLength      = h;
    tagstruct.ImageWidth       = w;
    tagstruct.Photometric      = 1;
    tagstruct.BitsPerSample    = bits;
    tagstruct.SamplesPerPixel  = 1;
    tagstruct.PlanarConfiguration = 1;
    tagstruct.Compression      = 1;    % or LZW for compression
    tagstruct.RowsPerStrip     = min(h, 256);              % safe strip height
    tagstruct.Software         = 'MATLAB';

    % Write first frame + directories
    updateEvery = max(1, round(T/100));  % ~100 updates
    for k = 1:T
        if k==1
            t.setTag(tagstruct);
        else
            t.writeDirectory();
            t.setTag(tagstruct);
        end
        % Ensure 2-D slice
        frame = fra(:,:,k);
        if ~isa(frame, class(fra))
            frame = cast(frame, class(fra));
        end
        t.write(frame);

        % progress / cancel
        if mod(k, updateEvery)==0 || k==T
            if getappdata(hWait,'canceling')
                warning('save_stack_big_tiff:Canceled','User canceled save at frame %d/%d.',k,T);
                break;  % onCleanup will still close TIFF; partial file will contain written frames
            end
            waitbar(k/T, hWait, sprintf('Saving frames to BigTIFF... %d/%d',k,T));
        end
    end

    % Close explicitly (also handled by onCleanup)
    try t.close(); catch, end
    try delete(hWait); catch, end

    fprintf('Saved %s (%dx%dx%d, %d-bit)\n', tifFile, w, h, T, bits);

    % ---- local helper
    function close_tiff_quiet(tt)
        try
            if ~isempty(tt) && isa(tt,'Tiff')
                tt.close();
            end
        catch
        end
    end
end
