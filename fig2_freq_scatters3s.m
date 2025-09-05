%% ===== 11 XY scatters by dominant-frequency bin, with 50 µm scale bars =====
% needs: domFreqs_by_fish, positions_by_fish (or resultsFiles), useFish, edges
% If positions_by_fish is missing, we’ll try to load positions from each MAT.

col_orange = [1 0.5 0];
gray80     = 0.6*[1 1 1];

% ---- get positions per fish if not provided ----
if ~exist('positions_by_fish','var') || isempty(positions_by_fish)
    positions_by_fish = cell(numel(resultsFiles),1);

    for f = find(useFish(:))'
        L = load(resultsFiles{f});        % struct of variables in the MAT

        % unwrap if saved as variable 'S'
        if isfield(L,'S') && isstruct(L.S)
            L = L.S;
        end

        % pull the 'results' container if present
        if isfield(L,'results') && isstruct(L.results)
            R = L.results;
        else
            R = L;  % fallback: maybe results fields are at top level
        end

        % grab positions_nodoubles explicitly
        if isfield(R,'positions_nodoubles') && ~isempty(R.positions_nodoubles)
            pos = double(R.positions_nodoubles);
        else
            error('positions_nodoubles not found in %s. Available fields in results: %s', ...
                  resultsFiles{f}, strjoin(fieldnames(R), ', '));
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

%% ---- stack across selected fish ----
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

% ---- define exactly 11 bins over [0, 0.10] Hz ----
binEdges11   = linspace(0, 0.10, 12);
binCenters11 = binEdges11(1:end-1) + diff(binEdges11)/2;

% ---- plotting params for scalebar ----
scale_um     = 50;     % draw a 50 µm bar
um_per_unit  = 1;      % change if your coordinates are not in µm (e.g., um_per_unit = 0.78)
outDir = '\\home.ansatt.ntnu.no\emiliajm\Documents\MATLAB\huc_slow_osc';
if ~exist(outDir,'dir'), mkdir(outDir); end

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

    % outside: gray 10%
    if any(outBin)
        scatter(ax, ALL_pos(outBin,1), ALL_pos(outBin,2), 18, gray80, 'filled', ...
                'MarkerFaceAlpha',0.10, 'MarkerEdgeAlpha',0.10);
    end

    % inside: density-colored (scatter_kde), ~70% opacity
    if any(inBin)
        h = scatter_kde(ALL_pos(inBin,1), ALL_pos(inBin,2), ...
            'filled','MarkerSize',20, ...
            'MarkerFaceAlpha',0.70,'MarkerEdgeAlpha',0.70);
        colormap(ax, orange_map(256));

        if b == 1
            climFixed = caxis(ax);
        else
            set(ax,'CLim',climFixed);
        end
    end

    % orientation & styling
    set(ax,'ZDir','reverse','YDir','reverse','XDir','reverse');
    axis(ax,'equal','tight');
    axis(ax,'off'); box(ax,'off'); grid(ax,'off');
    view(ax,90,90);

    % title & safe filename
    ttl      = sprintf('scatter_bin_%.3f_%.3fHz', lfs, lffs);
    ttl_file = regexprep(ttl, '\.', 'p');   % e.g., 0p009
    %title(ax, ttl, 'Interpreter','none');

    outFile = fullfile(outDir, [ttl_file '.png']);

    % ---- Save (axes handle!), with fallback ----
    try
        exportgraphics(ax, outFile, 'Resolution', 300, 'BackgroundColor','none');
    catch
        set(fig,'PaperPositionMode','auto','Color','white'); % fallback uses white bg
        print(fig, outFile, '-dpng', '-r300');
    end

    close(fig);
end


%% -------- helper: extract positions from a results struct --------
