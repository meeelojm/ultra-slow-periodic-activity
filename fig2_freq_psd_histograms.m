%% ---------------------- Config ----------------------
baseFolder = 'Z:\mh-kin\yaksi5\anna\Data\Processed 2P data\mGluR multimodal\cell_detect_data\';
folderPattern = '^[A-Za-z0-9]{4}_Wt$';   % or 'Het'
origStartFrame = 200; origFps = 2.5;     % where your analysis window starts in the *raw* movies

% Frequency analysis params (use your preferred ones)
windowSize     = 256;
maxPeakCount   = 15;
minPeakHeight  = 2;
frequencyRange = [0.01 0.10];            % wider for storage; you can restrict when plotting

%% ---------------------- 1) Find fish folders & files ----------------------
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    [~, nm] = fileparts(allFolders{i});
    if ~isempty(regexp(nm, folderPattern, 'once'))
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
    end
end
fprintf('Found %d matching folders.\n', numel(validFolders));

resultsFiles = {};
fishLabel    = {};
for i = 1:numel(validFolders)
    D = dir(fullfile(validFolders{i}, '*_Results_dff.mat'));       % ⚠️ pattern
    if ~isempty(D)
        resultsFiles{end+1} = fullfile(validFolders{i}, D(1).name); %#ok<SAGROW>
        [~, fishLabel{end+1}] = fileparts(validFolders{i});         %#ok<SAGROW>
    else
        warning('No *_Results_dff.mat in %s', validFolders{i});
    end
end

%% ---------------------- 2) Per-fish analysis & storage ----------------------
nF = numel(resultsFiles);

% Cells (as you asked)
domFreqs_by_fish  = cell(nF,1);     % each: [nNeurons x 1]
powerSpecs_by_fish= cell(nF,1);     % each: [nNeurons x nFreqs]
freqBins_by_fish  = cell(nF,1);     % each: [1 x nFreqs]
fps_by_fish       = nan(nF,1);
type_by_fish      = strings(nF,1);  % "in_vivo" / "explant" (guessed from path)

% Optional convenient struct mirror
fish(nF) = struct('id','','path','','type','','fps',NaN, ...
                  'domFreqs',[],'powerSpecs',[],'freqBins',[]);

for f = 1:nF
    S = load(resultsFiles{f});

    % --------- Pull data from the MAT (⚠️ rename if your fields differ) ---------
    % Try common names safely:
    % -------- DFF extraction (always here) --------
    if isfield(S,'results') && isfield(S.results,'DV_DFFmovwindow')
        dff_red = S.results.DV_DFFmovwindow;       % numeric matrix
    else
        error('Expected S.results.DV_DFFmovwindow in %s', resultsFiles{f});
    end
    
    % make sure it’s double and oriented as [cells x frames]
    dff_red = double(dff_red);
    if size(dff_red,1) > size(dff_red,2)
        dff_red = dff_red.';                         % transpose if [frames x cells]
    end

    % FPS
    fps = 2.5;
    % if isfield(S,'fps'), fps = S.fps; end
    % if isnan(fps) && isfield(S,'ops') && isfield(S.ops,'fs'), fps = S.ops.fs; end
    % if isnan(fps) && isfield(S,'meta') && isfield(S.meta,'fps'), fps = S.meta.fps; end
    % if isnan(fps)
    %     fps = 2.5;
    %     %error('Could not determine fps for %s', resultsFiles{f});
    % end
    fps_by_fish(f) = fps;

    % Analysis window
    startSec       = origStartFrame / origFps;
    startFrameFreq = max(1, round(startSec * fps));
    endFrameFreq   = min(size(dff_red,2), startFrameFreq + round(240*fps)); % 4 min
    numLagsFreq    = round(120 * fps);                                      % 2 min

    % --------- Your function (returns per-neuron measures) ---------
    [autocorrs, domFreqs, powerSpecs, freqBins, isoIdx] = ...
        calcsort_autocorr_freq_analysis_v3( ...
            dff_red(:, startFrameFreq:endFrameFreq), ...
            numLagsFreq, fps, ...
            windowSize, maxPeakCount, minPeakHeight, frequencyRange);

    % --------- Store (cells + struct) ---------
    domFreqs_by_fish{f}   = domFreqs(:);
    powerSpecs_by_fish{f} = powerSpecs;
    freqBins_by_fish{f}   = freqBins;

    fish(f).id        = fishLabel{f};
    fish(f).path      = resultsFiles{f};
    fish(f).fps       = fps;
    fish(f).domFreqs  = domFreqs(:);
    fish(f).powerSpecs= powerSpecs;
    fish(f).freqBins  = freqBins;

    % crude condition from path text; refine if you have a better tag
    pth = lower(resultsFiles{f});
    if contains(pth,'explant')
        fish(f).type = 'explant'; type_by_fish(f) = "explant";
    else
        fish(f).type = 'in_vivo'; type_by_fish(f) = "in_vivo";
    end
end

% Save once so you can reload later without recomputing
save(fullfile(baseFolder,'per_fish_domfreqs_powerSpecs.mat'), ...
     'domFreqs_by_fish','powerSpecs_by_fish','freqBins_by_fish', ...
     'fps_by_fish','type_by_fish','fish','resultsFiles','fishLabel','-v7.3');


%% ================== PLOTTING PARAMS ==================
edges   = 0:0.005:0.10;                   % same x-axis as before
centers = edges(1:end-1) + diff(edges)/2;
col     = [1 0.5 0];                      % orange (change per group if you want)

% choose which fish to include
groupSel = "in_vivo";                      % "in_vivo", "explant", or "all"
if exist('type_by_fish','var') && groupSel~="all"
    useFish = (type_by_fish==groupSel);
else
    useFish = true(numel(domFreqs_by_fish),1); % fall back to all
end

%% ================== 1) % HISTOGRAM + PER-FISH SCATTER ==================
% pooled dom freqs (within edges and finite)
allDom = vertcat(domFreqs_by_fish{useFish});
allDom = allDom(~isnan(allDom) & allDom>=edges(1) & allDom<=edges(end));

figure('Color','w'); hold on;
h = histogram(allDom,'BinEdges',edges,'DisplayStyle','stairs','LineWidth',3, ...
              'EdgeColor',col,'Normalization','percentage'); % fraction                                % show as %
xlabel('Dominant frequency (Hz)'); ylabel('Neuron percentage');
title('Distribution of dominant frequencies (pooled)');

% overlay per-fish percentages as scatter
for f = find(useFish(:))'
    v = domFreqs_by_fish{f};
    v = v(~isnan(v) & v>=edges(1) & v<=edges(end));
    if isempty(v), continue; end
    frac = histcounts(v, edges, 'Normalization','percentage'); % 1x(nBins)
    scatter(centers, frac, 20, col, 'filled', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
end
box on; grid on; xlim([edges(1) edges(end)]);

%% ================== POWER @ DOM FREQ VS FREQUENCY (OWN FIGURE) ==================
% Assumes `edges`, `centers`, `useFish`, `col` exist

normMode = 'area';   % 'area' | 'max' | 'none'

nBins  = numel(centers);
perFishMedP = nan(sum(useFish), nBins);     % rows = fish, cols = freq bins
allBinIdx = []; 
allPeakP  = [];

fishList = find(useFish(:))';
row = 0;

for f = fishList
    dom = domFreqs_by_fish{f};            if isempty(dom), continue; end
    dom = dom(:);                          % n x 1
    P   = powerSpecs_by_fish{f};           if isempty(P),   continue; end
    F   = freqBins_by_fish{f};             if isempty(F),   continue; end
    F   = F(:).';                          % 1 x m

    % ---- orient/align P with F ----
    if size(P,2) ~= numel(F) && size(P,1) == numel(F), P = P.'; end
    m = min(size(P,2), numel(F));          % common # of freq bins
    P = P(:,1:m);  F = F(1:m);
    n = min(size(P,1), numel(dom));        % common # of neurons
    if n==0, continue; end
    P = P(1:n,:);  dom = dom(1:n);

    % ---- convert F to Hz if needed (rad/sample or cycles/sample) ----
    fps = NaN;
    if exist('fps_by_fish','var') && numel(fps_by_fish) >= f && ~isnan(fps_by_fish(f))
        fps = fps_by_fish(f);
    elseif exist('fish','var') && numel(fish) >= f && isfield(fish,'fps') && ~isempty(fish(f).fps)
        fps = fish(f).fps;
    end
    if ~isnan(fps)
        mxF = max(F);
        if mxF <= pi + 1e-6           % rad/sample -> Hz
            F = F * (fps/(2*pi));
        elseif mxF <= 0.5 + 1e-6      % cycles/sample -> Hz
            F = F * fps;
        else
            % assume already in Hz
        end
    end

    % ---- NORMALIZE POWER (choose ONE via normMode) ----
    switch lower(normMode)
        case 'area'   % unit-area within plotting band; fallback to whole band if empty
            idxF = find(F >= edges(1) & F <= edges(end));
            if ~isempty(idxF)
                Pband = P(:, idxF);
                areas = trapz(F(idxF), Pband, 2);              % integrate per neuron
                P(:, idxF) = Pband ./ max(eps, areas);
            else
                areasAll = trapz(F, P, 2);
                P = P ./ max(eps, areasAll);
            end
        case 'max'    % per-neuron max = 1
            P = P ./ max(eps, max(P, [], 2));
        case 'none'
            % do nothing
    end

    % ---- nearest freq-bin to each neuron's dominant frequency ----
    D = abs(dom - F);                       % n x m (implicit expansion)
    D(~isfinite(D)) = inf;
    [~, idxNearest] = min(D, [], 2);
    idxNearest = idxNearest(:);             % column

    % robust row-wise indexing (avoid sub2ind shape issues)
    peakP = arrayfun(@(i) P(i, idxNearest(i)), 1:n).';   % n x 1

    % ---- bin + collect ----
    keep = isfinite(dom) & isfinite(peakP) & dom>=edges(1) & dom<=edges(end);
    if ~any(keep), continue; end
    b = discretize(dom(keep), edges);

    allBinIdx = [allBinIdx; b(:)];
    allPeakP  = [allPeakP;  peakP(keep)];

    row = row + 1;
    perFishMedP(row,:) = accumarray(b(:), peakP(keep), [nBins,1], @median, NaN).';
end
perFishMedP = perFishMedP(1:row,:);

% ---- pooled median across all neurons ----
if isempty(allBinIdx)
    medPooled = nan(1,nBins);
else
    medPooled = accumarray(allBinIdx, allPeakP, [nBins,1], @median, NaN).';
end

% ------------------- plot (own figure) -------------------
figure('Color','w'); hold on;
% plot(centers, medPooled, '-', 'LineWidth', 3, 'Color', col);   % pooled line (optional)
for r = 1:size(perFishMedP,1)                                   % per-fish dots
    scatter(centers, perFishMedP(r,:), 20, col, 'filled', ...
            'MarkerFaceAlpha',0.55, 'MarkerEdgeAlpha',0.55);
end

xlabel('Dominant frequency (Hz)');

% y-label based on normalization mode
switch lower(normMode)
    case 'area'
        ylabel('Power density at dominant frequency (1/Hz)');
    case 'max'
        ylabel('Relative power at dominant frequency (a.u., max = 1)');
    otherwise
        ylabel('Power at dominant frequency (a.u.)');
end

title(sprintf('Power vs dominant frequency (%s-normalized)', normMode));
grid on; box on; xlim([edges(1) edges(end)]);

% robust y-limit to avoid axis blow-ups from stragglers
yAll = [perFishMedP(:); medPooled(:)];
if all(isnan(yAll))
    ylim([0 1]);
else
    ylim([0, prctile(yAll(~isnan(yAll)), 99)]);
end
%% ================== (c) Mixed-effects: does sharpness (Q) increase with frequency? ==================
% Q-factor per neuron = f0 / FWHM around the dominant-frequency peak

fishList = find(useFish(:))';

Q_all   = [];
f0_all  = [];
fish_all= [];
type_all= strings(0,1);

for f = fishList
    dom = domFreqs_by_fish{f};     if isempty(dom), continue; end
    dom = dom(:);
    P   = powerSpecs_by_fish{f};    if isempty(P),   continue; end
    F   = freqBins_by_fish{f};      if isempty(F),   continue; end
    F   = F(:).';

    % --- orient/align with F ---
    if size(P,2) ~= numel(F) && size(P,1)==numel(F), P = P.'; end
    m = min(size(P,2), numel(F));   P = P(:,1:m);  F = F(1:m);
    n = min(size(P,1), numel(dom)); P = P(1:n,:);  dom = dom(1:n);

    % --- convert F to Hz if needed ---
    fps = NaN;
    if exist('fps_by_fish','var') && numel(fps_by_fish)>=f && ~isnan(fps_by_fish(f))
        fps = fps_by_fish(f);
    elseif exist('fish','var') && numel(fish)>=f && isfield(fish,'fps') && ~isempty(fish(f).fps)
        fps = fish(f).fps;
    end
    if ~isnan(fps)
        mxF = max(F);
        if mxF <= pi + 1e-6,   F = F * (fps/(2*pi));   % rad/sample -> Hz
        elseif mxF <= 0.5 + 1e-6, F = F * fps;        % cycles/sample -> Hz
        end
    end

    % --- restrict to analysis band ---
    inBand = dom>=edges(1) & dom<=edges(end);
    if ~any(inBand), continue; end
    dom = dom(inBand);
    P   = P(inBand,:);

    % --- compute Q per neuron (FWHM around nearest-bin peak) ---
    n = numel(dom);
    Q  = nan(n,1);
    f0 = nan(n,1);
    for i = 1:n
        f0i = dom(i);
        if ~isfinite(f0i), continue; end

        % spectrum for neuron i, normalized to its max (FWHM is scale-free)
        Si = P(i,:);
        mx = max(Si);
        if ~isfinite(mx) || mx<=0, continue; end
        Si = Si ./ mx;

        % nearest index to f0
        [~, k0] = min(abs(F - f0i));

        % find 0.5 crossings left/right with linear interpolation
        % left
        kL = k0;
        while kL>1 && Si(kL) > 0.5, kL = kL-1; end
        fL = NaN;
        if kL < k0
            fL = interp1([Si(kL) Si(kL+1)], [F(kL) F(kL+1)], 0.5, 'linear', 'extrap');
        end
        % right
        kR = k0;
        while kR<numel(F) && Si(kR) > 0.5, kR = kR+1; end
        fR = NaN;
        if kR > k0 && kR<=numel(F)
            fR = interp1([Si(kR-1) Si(kR)], [F(kR-1) F(kR)], 0.5, 'linear', 'extrap');
        end

        if isfinite(fL) && isfinite(fR) && fR>fL
            bw  = fR - fL;         % FWHM (Hz)
            Q(i)= f0i / bw;
            f0(i)= f0i;
        end
    end

    % collect
    keep = isfinite(Q) & isfinite(f0) & Q>0;
    Q_all    = [Q_all;  Q(keep)];
    f0_all   = [f0_all; f0(keep)];
    fish_all = [fish_all; repmat(string(fishLabel{f}), sum(keep), 1)];
    if exist('type_by_fish','var') && numel(type_by_fish)>=f
        type_all = [type_all; repmat(type_by_fish(f), sum(keep), 1)];
    else
        type_all = [type_all; repmat("unknown", sum(keep),1)];
    end
end

% ---- assemble table for LME ----
T = table( categorical(fish_all), f0_all, Q_all, categorical(type_all), ...
           'VariableNames', {'fish','f0','Q','type'});

% log-transform for stability (Q, f0 are positive)
T.logQ  = log(T.Q);
T.logf0 = log(T.f0);

% remove any remaining non-finite rows
T = T(all(isfinite(T{:,{'logQ','logf0'}}),2), :);

% ---- mixed-effects: logQ ~ logf0 + (1|fish) ----
lme0 = fitlme(T, 'logQ ~ 1 + (1|fish)');
lme1 = fitlme(T, 'logQ ~ logf0 + (1|fish)');

% report slope + SE + p, and LRT vs null
coef  = lme1.Coefficients;
ix    = strcmp(coef.Name,'logf0');
beta  = coef.Estimate(ix);  se = coef.SE(ix);  p  = coef.pValue(ix);

fprintf('Mixed-effects result (logQ ~ logf0 + (1|fish)):\n');
fprintf('  slope(logf0) = %.3f ± %.3f, p = %.3g\n', beta, se, p);

cmp = compare(lme0, lme1);   % likelihood-ratio test
fprintf('  LRT vs null: chi2(1)=%.3f, p=%.3g\n', cmp{2,'LRStat'}, cmp{2,'pValue'});

% Optional: add type as a covariate or interaction
% lme2 = fitlme(T, 'logQ ~ logf0 + type + (1|fish)');
% lme3 = fitlme(T, 'logQ ~ logf0*type + (1|fish)');
% compare(lme1,lme2), compare(lme2,lme3);

