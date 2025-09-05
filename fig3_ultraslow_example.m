%% ---- stack across selected fish (unchanged) ----
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

%% ---- band filter: 0.02–0.08 Hz ----
fmin = 0.05; fmax = 0.07;
inBand  = isfinite(ALL_dom) & ALL_dom >= fmin & ALL_dom <= fmax;
outBand = isfinite(ALL_dom) & ~inBand;

%% ---- plotting params ----
scale_um     = 50;   % 50 µm scale bar
um_per_unit  = 1;    % change if your coords are not in µm
if ~exist('gray80','var') || isempty(gray80), gray80 = 0.8*[1 1 1]; end

%% ---- plot (replicates your per-bin panel style, but single band) ----
figure('Color','w','Renderer','painters','Position',[100 80 520 600]); hold on;

% outside band: light gray, ~10% opacity
if any(outBand)
    scatter(ALL_pos(outBand,1), ALL_pos(outBand,2), 18, gray80, 'filled', ...
            'MarkerFaceAlpha',0.10, 'MarkerEdgeAlpha',0.10);
end

% inside band: density-colored (scatter_kde), ~70% opacity
if any(inBand)
    h = scatter_kde(ALL_pos(inBand,1), ALL_pos(inBand,2), ...
        'filled','MarkerSize',20, 'MarkerFaceAlpha',0.70,'MarkerEdgeAlpha',0.70);
    % keep your orange density shades
    if exist('orange_map','file')
        colormap(gca, orange_map(256));
    else
        colormap(gca, parula(256)); % fallback if orange_map not on path
    end
else
    warning('No neurons in the %.2f–%.2f Hz band.', fmin, fmax);
end

% same orientation & cosmetics as your original
set(gca,'ZDir','reverse','YDir','reverse','XDir','reverse');
axis equal tight; axis off; box off; grid off;

% 50 µm L-shaped scale bar (your helper)
if exist('draw_scalebar_xy','file')
    draw_scalebar_xy(gca, scale_um, um_per_unit);
end

view(90,90);   % top view
title(sprintf('%.2f–%.2f Hz (N=%d)', fmin, fmax, nnz(inBand)));

%% Params you already have
% Assumes you already have:
% resultsFiles (cellstr of MAT paths), useFish (logical), 
% fps, windowSize, maxPeakCount, minPeakHeight, frequencyRange,
% numLagsFreq, startFrameFreq, endFrameFreq

fmin = 0.05; 
fmax = 0.07;

AUTOCORRS = cell(numel(resultsFiles),1);
DOMFREQS  = cell(numel(resultsFiles),1);
POWERS    = cell(numel(resultsFiles),1);
FBINS     = cell(numel(resultsFiles),1);
ISOIDX    = cell(numel(resultsFiles),1);
N_INBAND  = zeros(numel(resultsFiles),1);

fishList = find(useFish(:))';

for ii = 1:numel(fishList)
    f = fishList(ii);
    matPath = resultsFiles{f};
    if ~isfile(matPath), warning('Fish %d: missing file: %s', f, matPath); continue; end

    % -------- Load this fish --------
    S = load(matPath);
    if ~isfield(S,'results') || ~isfield(S.results,'DV_DFFmovwindow')
        warning('Fish %d: expected S.results.DV_DFFmovwindow in %s', f, matPath);
        continue;
    end
    dff_red = double(S.results.DV_DFFmovwindow);
    % Make [cells x frames]
    if size(dff_red,1) > size(dff_red,2), dff_red = dff_red.'; end

    % Sanity: restrict to requested window
    dfw = dff_red(:, startFrameFreq:endFrameFreq);
    if isempty(dfw) || size(dfw,2) < 2
        warning('Fish %d: empty/short window after indexing.', f);
        continue;
    end

    % -------- First pass: dom freqs for mask --------
    [autoc_all, dom_all, pow_all, fBins_all, iso_all] = ...
        calcsort_autocorr_freq_analysis_v3( ...
            dfw, numLagsFreq, fps, ...
            windowSize, maxPeakCount, minPeakHeight, frequencyRange);

    dom_all = dom_all(:);
    n = min(size(dfw,1), numel(dom_all));
    dom_all = dom_all(1:n);
    dfw     = dfw(1:n,:);

    % Optional: intersect with "good" (isoIdx) if you normally do that
    good = true(n,1);
    if ~isempty(iso_all)
        good = false(n,1);
        keep = iso_all(iso_all>=1 & iso_all<=n);
        good(keep) = true;
    end

    inBand = isfinite(dom_all) & dom_all >= fmin & dom_all <= fmax & good;
    N_INBAND(f) = nnz(inBand);
    if N_INBAND(f) == 0
        fprintf('Fish %d: no neurons in %.2f–%.2f Hz. Skipping.\n', f, fmin, fmax);
        continue;
    end

    % -------- Second pass: run ONLY on in-band neurons --------
    [autocorrs, domFreqs, powerSpecs, freqBins, isoIdx] = ...
        calcsort_autocorr_freq_analysis_v3( ...
            dfw(inBand, :), numLagsFreq, fps, ...
            windowSize, maxPeakCount, minPeakHeight, frequencyRange);

    % Label the figure that the function just produced (if any)
    try
        fig = gcf;
        if ishghandle(fig)
            set(fig, 'Name', sprintf('Fish %d — domFreq %.2f–%.2f Hz (N=%d)', f, fmin, fmax, N_INBAND(f)));
            sgtitle(sprintf('Fish %d — domFreq %.2f–%.2f Hz (N=%d)', f, fmin, fmax, N_INBAND(f)), 'FontWeight','bold');
        end
    catch
        % If the function doesn’t draw a figure, nothing to do.
    end

    % Store outputs (optional)
    AUTOCORRS{f} = autocorrs;
    DOMFREQS{f}  = domFreqs;
    POWERS{f}    = powerSpecs;
    FBINS{f}     = freqBins;
    ISOIDX{f}    = isoIdx;
end

% Quick summary in console
fprintf('\nIn-band counts per fish (%.2f–%.2f Hz):\n', fmin, fmax);
disp(table(fishList(:), N_INBAND(fishList(:)), 'VariableNames', {'Fish','N_inBand'}));
