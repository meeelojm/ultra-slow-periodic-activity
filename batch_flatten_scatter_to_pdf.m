function batch_flatten_scatter_to_pdf(srcDir, outDir, pattern, dpi, bg)
% Flatten .fig files to lightweight (rasterized) PDFs.
% Usage:
%   srcDir  = 'C:\path\to\figs';
%   outDir  = 'C:\path\to\exports';
%   pattern = 'scatter_*.fig';
%   dpi     = 300;          % 150–300 good
%   bg      = 'white';      % 'white' or 'none' (PDF shows 'none' as transparent)

if nargin < 3 || isempty(pattern), pattern = 'scatter_*.fig'; end
if nargin < 4 || isempty(dpi),     dpi     = 300;             end
if nargin < 5 || isempty(bg),      bg      = 'white';         end
if ~exist(outDir,'dir'), mkdir(outDir); end

files = dir(fullfile(srcDir, pattern));
fprintf('Found %d figure(s). Exporting to %s\n', numel(files), outDir);

for k = 1:numel(files)
    inPath = fullfile(files(k).folder, files(k).name);
    [~, base] = fileparts(inPath);
    pdfPath = fullfile(outDir, [base '.pdf']);

    try
        fig = openfig(inPath, 'invisible');
        set(fig, 'Renderer','opengl', 'InvertHardcopy','off');  % keep bg

        % Preferred: exportgraphics (raster inside PDF)
        ok = false;
        try
            exportgraphics(fig, pdfPath, ...
                'ContentType','image', ...
                'BackgroundColor', bg, ...
                'Resolution', dpi);
            ok = true;
        catch
            % Older MATLABs may not support PDF in exportgraphics
        end

        if ~ok
            % Fallback: print -dpdf -opengl with page sized to figure @ DPI
            pos = get(fig, 'Position');              % [x y w h] in pixels
            wpx = max(1, round(pos(3)));
            hpx = max(1, round(pos(4)));
            win = wpx / dpi;                         % inches
            hin = hpx / dpi;

            set(fig, 'PaperUnits','inches', ...
                     'PaperPosition',[0 0 win hin], ...
                     'PaperSize',[win hin], ...
                     'Color', bg_if_needed(bg));     % set figure bg

            print(fig, pdfPath, '-dpdf', ['-r' num2str(dpi)], '-opengl');
        end

        close(fig);
        fprintf('✓ %s\n', base);

    catch ME
        fprintf(2, '✗ %s — %s\n', base, ME.message);
        if exist('fig','var') && ishghandle(fig), close(fig); end
    end
end
fprintf('Done.\n');

    function c = bg_if_needed(b)
        % map 'none' to actual figure color for print fallback
        if ischar(b) && strcmpi(b,'none')
            c = 'white';
        else
            c = b;
        end
    end
end
