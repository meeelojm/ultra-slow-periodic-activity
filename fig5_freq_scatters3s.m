%% ===== 11 XY scatters by dominant-frequency bin (CCN2A, MATLAB blue) =====
% needs: domFreqs_by_fish, useFish, resultsFiles (if positions_by_fish not given)

% --- Fig 2 binning (reuse if available) ---
if exist('binEdges11','var') && numel(binEdges11)>=2
    binEdges11 = binEdges11(:).';
elseif exist('edges','var') && numel(edges)>=2
    binEdges11 = edges(:).';
else
    binEdges11 = linspace(0, 0.10, 12);          % fallback: 11 bins in [0,0.10]
end
binCenters11 = binEdges11(1:end-1) + diff(binEdges11)/2;

% colors
col_blue = [0 0.4470 0.7410];                    % MATLAB blue
gray80   = 0.6*[1 1 1];

% --- get positions per fish if not provided ---
if ~exist('positions_by_fish','var') || isempty(positions_by_fish)
    positions_by_fish = cell(numel(resultsFiles),1);

    for f = 1:size(useFish,1)
        L = load(resultsFiles{f});
        if isfield(L,'S') && isstruct(L.S), L = L.S; end
        if isfield(L,'results') && isstruct(L.results), R = L.results; else, R = L; end

        if isfield(R,'position') && ~isempty(R.position)
            pos = double(R.position);
        else
            error('positions_nodoubles not found in %s', resultsFiles{f});
        end

        % enforce Nx3
        if size(pos,2) < 2
            error('Positions must have at least 2 columns (X,Y).');
        elseif size(pos,2) == 2
            pos(:,3) = 0;
        elseif size(pos,2) > 3
            pos = pos(:,1:3);
        end

        positions_by_fish{f} = pos;
    end
end

% --- stack across selected fish ---
ALL_pos = []; ALL_dom = [];
for f = find(useFish(:))'
    dom = domFreqs_by_fish{f}(:);
    pos = positions_by_fish{f};
    if isempty(dom) || isempty(pos), continue; end
    if size(pos,2) > 3, pos = pos(:,1:3); end
    n = min(numel(dom), size(pos,1));
    ALL_dom = [ALL_dom; dom(1:n)];
    ALL_pos = [ALL_pos; pos(1:n,:)];
end
if isempty(ALL_dom) || isempty(ALL_pos), error('No data to plot.'); end

% --- scalebar params ---
scale_um     = 50;     % draw a 50 µm bar
um_per_unit  = 1;      % change if your coords are not µm
if exist('outRoot','var') && isfolder(outRoot)
    outDir = fullfile(outRoot, 'figs_xy_bins');
else
    outDir = 'C:\Users\emiliajm\Documents\MATLAB\ccn2a_slow_osc_invivo';
end
if ~exist(outDir,'dir'), mkdir(outDir); end

% ===== put this BEFORE the for-loop =====
dens_scale = [];              % global [lo hi] used by all panels

% keep CLim consistent across panels
climFixed = [];

for b = 1:11
    lfs  = binEdges11(b);
    lffs = binEdges11(b+1);

    if b < 11
        inBin  = ALL_dom >= lfs & ALL_dom <  lffs;
    else
        inBin  = ALL_dom >= lfs & ALL_dom <= lffs;
    end
    outBin = ~inBin & isfinite(ALL_dom);

    fig = figure('Color','w','Renderer','painters','Position',[100 80 520 600]);
    ax  = axes('Parent',fig); hold(ax,'on');

    % ---------- OUTSIDE (gray background cloud) ----------
    okOut = outBin & isfinite(ALL_pos(:,1)) & isfinite(ALL_pos(:,2));
    if any(okOut)
        scatter(ax, ALL_pos(okOut,1), ALL_pos(okOut,2), 18, gray80, 'filled', ...
                'MarkerFaceAlpha',0.08, 'MarkerEdgeAlpha',0.00);
    end

    % ---------- INSIDE (density-colored MATLAB-blue) ----------
    if any(inBin)
        xi = ALL_pos(inBin,1); yi = ALL_pos(inBin,2);
        ok = isfinite(xi) & isfinite(yi);
        xi = xi(ok); yi = yi(ok);

        if numel(xi) < 3
            scatter(ax, xi, yi, 24, [0 0.4470 0.7410], 'filled', ...
                    'MarkerFaceAlpha',0.9, 'MarkerEdgeAlpha',0.9);
        else
            sx = iqr(xi)/1.349; if sx==0, sx = std(xi); end
            sy = iqr(yi)/1.349; if sy==0, sy = std(yi); end
            bw = 1.06 * [sx sy] * numel(xi)^(-1/5);
            dens = ksdensity([xi yi], [xi yi], 'Bandwidth', bw);

            if isempty(dens_scale)
                lo = prctile(dens,5); hi = prctile(dens,95);
                if ~(isfinite(lo)&&isfinite(hi)) || hi<=lo, lo=min(dens); hi=max(dens); end
                dens_scale = [lo hi];
            end
            lo = dens_scale(1); hi = dens_scale(2);

            dclr = (dens - lo) ./ max(eps, hi - lo);
            dclr = min(max(dclr,0),1);
            cmap = [linspace(0.85,0,256)', linspace(0.92,0.4470,256)', linspace(0.98,0.7410,256)'];
            idx  = 1 + round(dclr*(size(cmap,1)-1));
            colors = cmap(idx,:);

            scatter(ax, xi, yi, 24, colors, 'filled', ...
                    'MarkerFaceAlpha',0.9, 'MarkerEdgeAlpha',0.9);
        end
    end

    % styling, scalebar, save ... (unchanged)
    set(ax,'ZDir','reverse','YDir','reverse','XDir','reverse');
    axis(ax,'equal','tight'); axis(ax,'off'); view(ax,90,90);
    add_scalebar(ax, scale_um, um_per_unit, [0 0.4470 0.7410]);

    outFile = fullfile(outDir, sprintf('scatter_bin_%s_%sHz.png', ...
             strrep(sprintf('%.3f',lfs),'.','p'), strrep(sprintf('%.3f',lffs),'.','p')));
    exportgraphics(ax, outFile, 'Resolution', 300, 'BackgroundColor','none');
    % close(fig);
end


%% ===== helpers =====
function cmap = matlab_blue_map(n)
% light-to-MATLAB-blue gradient
    c1 = [0.85 0.92 0.98];                  % very light blue
    c2 = [0 0.4470 0.7410];                 % MATLAB blue
    t  = linspace(0,1,max(2,n))';
    cmap = [c1(1)+(c2(1)-c1(1))*t, ...
            c1(2)+(c2(2)-c1(2))*t, ...
            c1(3)+(c2(3)-c1(3))*t];
end

function add_scalebar(ax, scale_um, um_per_unit, col)
% draws a horizontal scalebar of length "scale_um" at bottom-left
    L = scale_um / um_per_unit;             % length in data units
    axlim = axis(ax);
    xspan = axlim(2)-axlim(1);
    yspan = axlim(4)-axlim(3);

    x0 = axlim(1) + 0.08*xspan;
    y0 = axlim(3) + 0.08*yspan;
    plot(ax, [x0 x0+L], [y0 y0], '-', 'Color', col, 'LineWidth', 2.5);
    text(ax, x0+L/2, y0 - 0.02*yspan, sprintf('%d \\mum', scale_um), ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'Color', col, 'FontWeight','bold');
end
