function batch_flatten_scatter_to_svg(srcDir, outDir, pattern, dpi, bg)
% Batch "flatten" .fig files into lightweight SVGs by embedding a PNG.
% Works on MATLAB versions that do NOT support SVG in exportgraphics.
%
% Usage:
%   srcDir  = 'C:\path\to\figs';
%   outDir  = 'C:\path\to\exports';
%   pattern = 'scatter_*.fig';
%   dpi     = 300;         % 150–300 good
%   bg      = 'none';      % or 'white'
%   batch_flatten_scatter_to_svg(srcDir, outDir, pattern, dpi, bg)

if nargin < 3 || isempty(pattern), pattern = 'scatter_*.fig'; end
if nargin < 4 || isempty(dpi),     dpi     = 300;             end
if nargin < 5 || isempty(bg),      bg      = 'none';          end
if ~exist(outDir,'dir'), mkdir(outDir); end

files = dir(fullfile(srcDir, pattern));
fprintf('Found %d figure(s). Exporting to %s\n', numel(files), outDir);

for k = 1:numel(files)
    inPath = fullfile(files(k).folder, files(k).name);
    [~, base] = fileparts(inPath);
    pngPath = fullfile(outDir, [base '.png']);
    svgPath = fullfile(outDir, [base '.svg']);

    try
        fig = openfig(inPath, 'invisible');
        set(fig,'Renderer','opengl');                   % consistent rasterization

        % Export a high-DPI PNG (this is the flattened image)
        try
            exportgraphics(fig, pngPath, 'BackgroundColor', bg, 'Resolution', dpi);
        catch
            % fallback for very old MATLAB: use print
            set(fig, 'PaperPositionMode','auto');
            print(fig, pngPath, ['-r' num2str(dpi)], '-dpng');
        end

        % Get pixel dimensions for the wrapper
        pos = get(fig,'Position');   % [left bottom width height] in px
        wpx = max(1, round(pos(3)));
        hpx = max(1, round(pos(4)));

        % Create an SVG that embeds the PNG as base64
        png2svg_embedded(pngPath, svgPath, wpx, hpx);

        close(fig);
        fprintf('✓ %s\n', base);

    catch ME
        fprintf(2,'✗ %s — %s\n', base, ME.message);
        if exist('fig','var') && ishghandle(fig), close(fig); end
    end
end

fprintf('Done.\n');
end

function png2svg_embedded(pngPath, svgPath, width_px, height_px)
% Wrap a PNG file into an SVG container using base64 data URI

% read file bytes
fid = fopen(pngPath,'r');
if fid < 0, error('Cannot open %s', pngPath); end
bytes = fread(fid, Inf, '*uint8'); fclose(fid);

% base64 encode (available in MATLAB R2016b+)
b64 = matlab.net.base64encode(bytes);

% write SVG
fid = fopen(svgPath,'w');
if fid < 0, error('Cannot write %s', svgPath); end
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n');
fprintf(fid,'<svg xmlns="http://www.w3.org/2000/svg" ');
fprintf(fid,'xmlns:xlink="http://www.w3.org/1999/xlink" ');
fprintf(fid,'width="%dpx" height="%dpx" viewBox="0 0 %d %d">\n', ...
        width_px, height_px, width_px, height_px);
fprintf(fid,'  <image x="0" y="0" width="%d" height="%d" xlink:href="data:image/png;base64,%s"/>\n', ...
        width_px, height_px, b64);
fprintf(fid,'</svg>\n');
fclose(fid);
end
