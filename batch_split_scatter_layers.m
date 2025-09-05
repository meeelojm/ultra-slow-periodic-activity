function batch_split_scatter_layers(srcDir, outDir, pattern, dpi, bgVec)
% Split heavy scatter figures into:
%  (A) vector PDF (axes + labels + ORANGE points only)
%  (B) raster PNG (GRAY background points only, transparent)
%
% Usage:
%   srcDir  = 'C:\path\to\figs';
%   outDir  = 'C:\path\to\exports';
%   pattern = 'scatter_*.fig';
%   dpi     = 300;        % 150–300 works well
%   bgVec   = 'white';    % background for the vector PDF ('white' or 'none')
%
% Notes:
% - If you can, tag objects in your plotting code:
%       set(hGray,   'Tag','bg_gray');
%       set(hOrange, 'Tag','fg_orange');
%   This function will use tags if present; otherwise it falls back to color heuristics.

if nargin < 3 || isempty(pattern), pattern = 'scatter_*.fig'; end
if nargin < 4 || isempty(dpi),     dpi     = 300;             end
if nargin < 5 || isempty(bgVec),   bgVec   = 'white';         end
if ~exist(outDir,'dir'), mkdir(outDir); end

files = dir(fullfile(srcDir, pattern));
fprintf('Found %d figure(s). Exporting to %s\n', numel(files), outDir);

for k = 1:numel(files)
    inPath = fullfile(files(k).folder, files(k).name);
    [~, base] = fileparts(inPath);
    pdfPath = fullfile(outDir, [base '_VECTOR_orange_axes.pdf']);
    pngPath = fullfile(outDir, [base '_RASTER_gray.png']);

    try
        fig = openfig(inPath, 'invisible');
        set(fig, 'InvertHardcopy','off');  % preserve colors
        ax  = findobj(fig,'Type','axes');  % all axes

        % --- Detect graphics objects ---
        % Preferred: tags
        hGray   = findobj(fig,'-regexp','Tag','^bg_gray$');
        hOrange = findobj(fig,'-regexp','Tag','^fg_orange$');

        % Fallback: color/alpha heuristics
        if isempty(hGray) || isempty(hOrange)
            gObjs = findall(fig, '-property','Visible'); % all drawable
            % Likely gray points (~0.6 gray, low alpha)
            hGray = [hGray; find_gray_like(gObjs)];
            % Likely orange points (~[1 .5 0], higher alpha)
            hOrange = [hOrange; find_orange_like(gObjs)];
            % Deduplicate
            hGray   = unique(hGray);
            hOrange = setdiff(unique(hOrange), hGray); % ensure disjoint
        end

        % ---- PASS A: VECTOR PDF (axes + labels + orange only) ----
        toggle_all(gObjs, 'off');     % hide everything
        set([ax; hOrange], 'Visible','on');

        % keep tick labels & titles visible
        for a = reshape(ax,1,[])
            set(get(a,'Title'), 'Visible','on');
            set(get(a,'XLabel'), 'Visible','on');
            set(get(a,'YLabel'), 'Visible','on');
            set(get(a,'ZLabel'), 'Visible','on');
        end

        % vector export (painters)
        set(fig,'Renderer','painters','Color', bg_map(bgVec));
        pos = get(fig,'Position'); w = pos(3); h = pos(4);
        set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 w h]/dpi, 'PaperSize',[w h]/dpi);
        print(fig, pdfPath, '-dpdf', ['-r' num2str(dpi)]);  % vector

        % ---- PASS B: RASTER PNG (gray mass only, transparent) ----
        toggle_all(gObjs, 'off');
        % hide axes entirely for a clean layer
        set(ax,'Visible','off');
        set([get(ax,'Title'); get(ax,'XLabel'); get(ax,'YLabel'); get(ax,'ZLabel')],'Visible','off');

        set(hGray,'Visible','on');
        set(fig,'Renderer','opengl');  % consistent rasterization

        % transparent PNG
        try
            exportgraphics(fig, pngPath, 'BackgroundColor','none', 'Resolution', dpi);
        catch
            % fallback: print (no true transparency -> premultiply over white)
            set(fig,'Color','white');
            print(fig, pngPath, ['-r' num2str(dpi)], '-dpng');
        end

        close(fig);
        fprintf('✓ %s\n', base);

    catch ME
        fprintf(2,'✗ %s — %s\n', base, ME.message);
        if exist('fig','var') && ishghandle(fig), close(fig); end
    end
end

fprintf('Done.\n');

% ---------- helpers ----------
    function toggle_all(hh, vis)
        for hi = reshape(hh,1,[])
            try, set(hi,'Visible',vis); end %#ok<TRYNC>
        end
    end

    function c = bg_map(b)
        if ischar(b) && strcmpi(b,'none')
            c = 'white';  % PDFs prefer a solid bg
        else
            c = b;
        end
    end

end

function h = find_gray_like(objs)
% Heuristic: gray-ish face/edge color & low alpha
h = [];
for o = reshape(objs,1,[])
    try
        col = get_color(o);
        aF = get_alpha(o,'MarkerFaceAlpha');
        aE = get_alpha(o,'MarkerEdgeAlpha');
        if ~isempty(col)
            isGray = max(abs(col - mean(col))) < 0.08 && mean(col) > 0.45 && mean(col) < 0.75;
            isFaint = (~isempty(aF) && aF <= 0.2) || (~isempty(aE) && aE <= 0.2);
            if isGray && (isFaint || isScatterLike(o))
                h = [h; o]; %#ok<AGROW>
            end
        end
    catch, end
end
end

function h = find_orange_like(objs)
% Heuristic: orange-ish & moderate/high alpha
h = [];
for o = reshape(objs,1,[])
    try
        col = get_color(o);
        aF  = get_alpha(o,'MarkerFaceAlpha'); aE = get_alpha(o,'MarkerEdgeAlpha');
        if ~isempty(col)
            isOrange = (col(1) > 0.85) && (col(2) > 0.35) && (col(3) < 0.20);
            isOpaque = (isempty(aF) || aF >= 0.5) && (isempty(aE) || aE >= 0.5);
            if isOrange && isOpaque
                h = [h; o]; %#ok<AGROW>
            end
        end
    catch, end
end
end

function tf = isScatterLike(o)
t = get(o,'Type');
tf = any(strcmpi(t, {'scatter','line','patch','surface'}));
end

function a = get_alpha(o, prop)
a = [];
if isprop(o, prop)
    try, a = get(o, prop); end
end
end

function col = get_color(o)
col = [];
% Try typical properties
props = {'MarkerFaceColor','MarkerEdgeColor','CData','Color'};
for p = 1:numel(props)
    if isprop(o, props{p})
        v = get(o, props{p});
        if isnumeric(v) && numel(v) >= 3
            v = v(1:3);
            if all(isfinite(v)) && all(v>=0 & v<=1)
                col = v(:).';
                return
            end
        end
    end
end
end
