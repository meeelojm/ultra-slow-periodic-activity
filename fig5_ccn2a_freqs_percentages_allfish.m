%% ================== CONFIG (editables) ==================
baseFolder = 'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';
saveName   = 'per_fish_domfreqs_powerSpecs.mat';

% analysis window in *seconds* relative to the start of the loaded traces
analysisStartSec = 0;        % if you want to skip an initial chunk, put seconds here
analysisDurSec   = 240;      % analyze 4 minutes per fish (adjust as needed)

% frequency-analysis params (same API as your function)
windowSize     = 1024;
maxPeakCount   = 15;
minPeakHeight  = 2;
frequencyRange = [0.01 0.10];

% plotting params (can reuse your previous plotting blocks)
edges   = 0:0.005:0.10;
centers = edges(1:end-1) + diff(edges)/2;
col = [0 0.4470 0.7410];


%% ================== 1) FIND FOLDERS WITH DATA ==================
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
resultsFiles = {};
metaFiles    = {};
for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    rf = fullfile(allFolders{i}, 'dffs_repact_respcells.mat');
    mf = fullfile(allFolders{i}, 'metadata_multimodal.mat');
    if exist(rf,'file')==2 && exist(mf,'file')==2
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
        resultsFiles{end+1} = rf;            %#ok<SAGROW>
        metaFiles{end+1}    = mf;            %#ok<SAGROW>
    end
end
fprintf('N fish folders: %d\n', numel(validFolders));

%% ================== 2) FPS PER FOLDER (ini/metadata fallback) ==================
fps_by_folder = nan(1, numel(validFolders));
for i = 1:numel(validFolders)
    % try metadata first
    fps_by_folder(i) = try_get_fps_from_metadata(metaFiles{i});
    % then try any .ini next
    if isnan(fps_by_folder(i))
        iniFiles = dir(fullfile(validFolders{i}, '*.ini'));
        if ~isempty(iniFiles)
            fps_by_folder(i) = try_get_fps_from_ini(fullfile(validFolders{i}, iniFiles(1).name));
        end
    end
    % default
    if isnan(fps_by_folder(i))
        fps_by_folder(i) = 14.64;  % your default for this dataset
    end
end

%% ================== 3) PER-FISH ANALYSIS & STORAGE ==================
nF = numel(resultsFiles);

domFreqs_by_fish   = cell(nF,1);     % [nNeurons x 1]
powerSpecs_by_fish = cell(nF,1);     % [nNeurons x nFreqs]
freqBins_by_fish   = cell(nF,1);     % [1 x nFreqs]
fps_by_fish        = nan(nF,1);
type_by_fish       = strings(nF,1);  % optional label if you want

fish(nF) = struct('id','','path','','type','','fps',NaN, ...
                  'domFreqs',[],'powerSpecs',[],'freqBins',[]);

for f = 1:nF
    Sres  = load(resultsFiles{f});
    Smeta = load(metaFiles{f});

    % -------- DFF extraction (robust) --------
    [dff, chosenPath] = find_dff_matrix(Sres);
    if isempty(dff)
        error('Could not locate a DFF-like matrix in %s', resultsFiles{f});
    end
    fprintf('[%02d/%02d] Using "%s" as dff in %s\n', f, nF, chosenPath, resultsFiles{f});

    % orient to [cells x frames]
    dff = double(dff);
    if size(dff,1) > size(dff,2)     % [frames x cells] -> transpose
        dff = dff.';
    end

    % FPS
    fps = fps_by_folder(f);
    fps_by_fish(f) = fps;

    % analysis window in frames
    startFrame = max(1, round(analysisStartSec * fps) + 1);
    endFrame   = min(size(dff,2), startFrame + round(analysisDurSec*fps) - 1);
    numLags    = round(0.5 * analysisDurSec * fps);  % ~ half-window for autocorr
    if endFrame <= startFrame
        warning('Fish %d: analysis window collapsed; skipping.', f);
        continue;
    end

    % -------- your function (returns per-neuron measures) --------
    [autocorrs, domFreqs, powerSpecs, freqBins, isoIdx] = ...
        calcsort_autocorr_freq_analysis_v3( ...
            dff(:, startFrame:endFrame), ...
            numLags, fps, ...
            windowSize, maxPeakCount, minPeakHeight, frequencyRange);

    % store
    domFreqs_by_fish{f}   = domFreqs(:);
    powerSpecs_by_fish{f} = powerSpecs;
    freqBins_by_fish{f}   = freqBins;

    % labels
    [~, fishID] = fileparts(validFolders{f});
    fish(f).id         = fishID;
    fish(f).path       = resultsFiles{f};
    fish(f).fps        = fps;
    fish(f).domFreqs   = domFreqs(:);
    fish(f).powerSpecs = powerSpecs;
    fish(f).freqBins   = freqBins;

    % crude type tag from folder name (edit if you have better labels)
    pth = lower(validFolders{f});
    if contains(pth,'explant')
        fish(f).type = 'explant'; type_by_fish(f) = "explant";
    else
        fish(f).type = 'in_vivo'; type_by_fish(f) = "in_vivo";
    end
end

% save once for quick reloading later
save(fullfile(baseFolder, saveName), ...
     'domFreqs_by_fish','powerSpecs_by_fish','freqBins_by_fish', ...
     'fps_by_fish','type_by_fish','fish','resultsFiles','metaFiles','validFolders','-v7.3');

%% ================== 4) QUICK PLOTS (same style as before) ==================
% choose which fish to include
if exist('binEdges11','var') && numel(binEdges11) >= 2
    edges = binEdges11(:).';                    % <- use Fig 2's edges
elseif ~exist('edges','var') || numel(edges) < 2
    edges = linspace(0, 0.10, 12);              % fallback: 11 bins
end
centers = edges(1:end-1) + diff(edges)/2;

% choose which fish to include
groupSel = "in_vivo";  % "in_vivo", "explant", or "all"
if groupSel ~= "all"
    useFish = (type_by_fish==groupSel);
else
    useFish = true(numel(domFreqs_by_fish),1);
end
% (a) histogram + per-fish scatter
allDom = vertcat(domFreqs_by_fish{useFish});
allDom = allDom(~isnan(allDom) & allDom>=edges(1) & allDom<=edges(end));

figure('Color','w'); hold on;
histogram(allDom,'BinEdges',edges,'DisplayStyle','stairs','LineWidth',3, ...
          'EdgeColor',col,'Normalization','percentage');
% vertical dotted lines at bin edges (to mirror Fig 2)
arrayfun(@(x) xline(x,':k','LineWidth',0.75), edges);

xlabel('Dominant frequency (Hz)'); ylabel('Neuron percentage');
title('Distribution of dominant frequencies (pooled)');

for f = find(useFish(:))'
    v = domFreqs_by_fish{f};
    v = v(~isnan(v) & v>=edges(1) & v<=edges(end));
    if isempty(v), continue; end
    frac = histcounts(v, edges, 'Normalization','percentage');
    scatter(centers, frac, 20, col, 'filled', ...
            'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.25);
end
box off; grid off; xlim([edges(1) edges(end)]);

% (b) power at dom freq vs dom freq (per-fish dots)
normMode = 'area';   % 'area' | 'max' | 'none'
nBins  = numel(centers);
perFishMedP = nan(sum(useFish), nBins);
allBinIdx = []; allPeakP = [];
fishList = find(useFish(:))'; row = 0;

for f = fishList
    dom = domFreqs_by_fish{f};     if isempty(dom), continue; end
    dom = dom(:);
    P   = powerSpecs_by_fish{f};   if isempty(P),   continue; end
    F   = freqBins_by_fish{f};     if isempty(F),   continue; end
    F   = F(:).';

    if size(P,2) ~= numel(F) && size(P,1) == numel(F), P = P.'; end
    m = min(size(P,2), numel(F));  P = P(:,1:m); F = F(1:m);
    n = min(size(P,1), numel(dom)); P = P(1:n,:); dom = dom(1:n);

    % convert F to Hz if needed
    fps = fps_by_fish(f);
    if ~isnan(fps)
        mxF = max(F);
        if mxF <= pi + 1e-6,      F = F * (fps/(2*pi));
        elseif mxF <= 0.5 + 1e-6, F = F * fps;
        end
    end

    % normalize power (robust area-normalization)
    switch lower(normMode)
        case 'area'
            idxF = find(F >= edges(1) & F <= edges(end));
            if numel(idxF) >= 2
                Pband = P(:, idxF);
                areas = trapz(F(idxF), Pband, 2);
                P(:, idxF) = Pband ./ max(eps, areas);
            elseif numel(idxF) == 1         % single-bin fallback
                df_eff = min([diff(F)'; edges(end)-edges(1)]);
                P(:, idxF) = P(:, idxF) ./ max(eps, P(:, idxF) * df_eff);
            else
                areasAll = trapz(F, P, 2);
                P = P ./ max(eps, areasAll);
            end
        case 'max'
            P = P ./ max(eps, max(P, [], 2));
        case 'none'
            % do nothing
    end

    % nearest freq-bin to each dom
    [~, idxNearest] = min(abs(dom - F), [], 2);
    idxNearest = idxNearest(:);
    peakP = arrayfun(@(i) P(i, idxNearest(i)), 1:n).';

    keep = isfinite(dom) & isfinite(peakP) & dom>=edges(1) & dom<=edges(end);
    if ~any(keep), continue; end
    b = discretize(dom(keep), edges);

    allBinIdx = [allBinIdx; b(:)];
    allPeakP  = [allPeakP;  peakP(keep)];

    row = row + 1;
    perFishMedP(row,:) = accumarray(b(:), peakP(keep), [nBins,1], @median, NaN).';
end
perFishMedP = perFishMedP(1:row,:);

if isempty(allBinIdx)
    medPooled = nan(1,nBins);
else
    medPooled = accumarray(allBinIdx, allPeakP, [nBins,1], @median, NaN).';
end

figure('Color','w'); hold on;
for r = 1:size(perFishMedP,1)
    scatter(centers, perFishMedP(r,:), 20, col, 'filled', ...
            'MarkerFaceAlpha',0.55, 'MarkerEdgeAlpha',0.55);
end
% same vertical bin lines as Fig 2
% arrayfun(@(x) xline(x,':k','LineWidth',0.75), edges);

xlabel('Dominant frequency (Hz)');
switch lower(normMode)
    case 'area', ylabel('Power density at dominant frequency (1/Hz)');
    case 'max',  ylabel('Relative power at dominant frequency (a.u., max=1)');
    otherwise,   ylabel('Power at dominant frequency (a.u.)');
end
title(sprintf('Power vs dominant frequency (%s-normalized)', normMode));
grid off; box off; xlim([edges(1) edges(end)]);
yAll = [perFishMedP(:); medPooled(:)];
if all(isnan(yAll)), ylim([0 1]); else, ylim([0, prctile(yAll(~isnan(yAll)), 99)]); end
%% ================== 5) Q-FACTOR LME ==================
Q_all=[]; f0_all=[]; fish_all=[]; type_all=strings(0,1);
fishList = find(useFish(:))';
for f = fishList
    dom = domFreqs_by_fish{f};     if isempty(dom), continue; end
    dom = dom(:);
    P   = powerSpecs_by_fish{f};   if isempty(P),   continue; end
    F   = freqBins_by_fish{f};     if isempty(F),   continue; end
    F   = F(:).';

    if size(P,2) ~= numel(F) && size(P,1)==numel(F), P = P.'; end
    m = min(size(P,2), numel(F));  P = P(:,1:m); F = F(1:m);
    n = min(size(P,1), numel(dom)); P = P(1:n,:); dom = dom(1:n);

    fps = fps_by_fish(f);
    if ~isnan(fps)
        mxF = max(F);
        if mxF <= pi + 1e-6,      F = F * (fps/(2*pi));
        elseif mxF <= 0.5 + 1e-6, F = F * fps;
        end
    end

    inBand = dom>=edges(1) & dom<=edges(end);
    if ~any(inBand), continue; end
    dom = dom(inBand);
    P   = P(inBand,:);

    Q  = nan(numel(dom),1);
    f0 = nan(numel(dom),1);
    for i = 1:numel(dom)
        f0i = dom(i);
        if ~isfinite(f0i), continue; end
        Si = P(i,:); mx = max(Si);
        if ~isfinite(mx) || mx<=0, continue; end
        Si = Si ./ mx;
        [~, k0] = min(abs(F - f0i));
        % left crossing
        kL = k0;  while kL>1 && Si(kL) > 0.5, kL = kL-1; end
        fL = NaN; if kL < k0, fL = interp1([Si(kL) Si(kL+1)], [F(kL) F(kL+1)], 0.5, 'linear','extrap'); end
        % right crossing
        kR = k0;  while kR<numel(F) && Si(kR) > 0.5, kR = kR+1; end
        fR = NaN; if kR > k0 && kR<=numel(F), fR = interp1([Si(kR-1) Si(kR)], [F(kR-1) F(kR)], 0.5, 'linear','extrap'); end
        if isfinite(fL) && isfinite(fR) && fR>fL
            bw   = fR - fL;
            Q(i) = f0i / bw;
            f0(i)= f0i;
        end
    end
    keep = isfinite(Q) & isfinite(f0) & Q>0;
    Q_all    = [Q_all;  Q(keep)];
    f0_all   = [f0_all; f0(keep)];
    fish_all = [fish_all; repmat(string(fish(f).id), sum(keep), 1)];
    type_all = [type_all; repmat(type_by_fish(f), sum(keep), 1)];
end

T = table( categorical(fish_all), f0_all, Q_all, categorical(type_all), ...
           'VariableNames', {'fish','f0','Q','type'});
T.logQ  = log(T.Q);
T.logf0 = log(T.f0);
T = T(all(isfinite(T{:,{'logQ','logf0'}}),2), :);

lme0 = fitlme(T, 'logQ ~ 1 + (1|fish)');
lme1 = fitlme(T, 'logQ ~ logf0 + (1|fish)');
coef = lme1.Coefficients; ix = strcmp(coef.Name,'logf0');
fprintf('Mixed-effects (logQ ~ logf0 + (1|fish)): slope = %.3f ± %.3f, p = %.3g\n', ...
        coef.Estimate(ix), coef.SE(ix), coef.pValue(ix));
cmp = compare(lme0,lme1);
fprintf('  LRT vs null: chi2(1)=%.3f, p=%.3g\n', cmp{2,'LRStat'}, cmp{2,'pValue'});

%% ================== HELPERS ==================
function fps = try_get_fps_from_metadata(mfile)
    fps = NaN;
    try
        S = load(mfile);
        fps = sniff_fps_in_struct(S);
    catch
        fps = NaN;
    end
end

function fps = sniff_fps_in_struct(S)
    fps = NaN;
    keys = lower(string(fieldnames(S)));
    for k = 1:numel(keys)
        val = S.(keys(k));
        if isstruct(val)
            fps = sniff_fps_in_struct(val);
            if ~isnan(fps), return; end
        elseif isnumeric(val) && isscalar(val)
            % common names
            if any(contains(keys(k), ["fps","fs","framerate","frame_rate","volume_rate","volumeRate","volume_rate_hz"]))
                fps = double(val); return;
            end
        end
    end
end

function fps = try_get_fps_from_ini(iniPath)
    fps = NaN;
    try
        txt = fileread(iniPath);
        % your regex
        tok = regexp(txt, 'volume\.rate\.\(in\.Hz\)\s*=\s*([0-9.]+)', 'tokens');
        if ~isempty(tok)
            fps = str2double(tok{1}{1});
            return;
        end
        % some scanners store as "frame_rate ="
        tok = regexp(txt, 'frame[_ ]?rate\s*=\s*([0-9.]+)', 'tokens');
        if ~isempty(tok)
            fps = str2double(tok{1}{1});
            return;
        end
    catch
        fps = NaN;
    end
end

function [dff, chosenPath] = find_dff_matrix(S)
    % 1) direct common paths first (fast path)
    tryPaths = { ...
        'results.DV_DFFmovwindow',' trace', ...
        'results.DFF', ...
        'results.dff', ...
        'results.dff_red', ...
        'results.zdff', ...
        'dff', 'DFF', 'DV_DFFmovwindow', 'zdff' ...
    };
    for i = 1:numel(tryPaths)
        [ok, A] = get_by_path(S, tryPaths{i});
        if ok && isnumeric(A) && ismatrix(A) && ~isscalar(A)
            dff = A; chosenPath = tryPaths{i}; return;
        end
    end

    % 2) heuristic search: collect numeric 2-D arrays and score
    cands = collect_numeric_2d(S, '');
    if isempty(cands)
        dff = []; chosenPath = '';
        return;
    end

    % --- make names a cell array of char, then lowercase safely ---
    namesChar  = cellfun(@char, {cands.path}, 'UniformOutput', false);
    namesLower = cellfun(@lower, namesChar, 'UniformOutput', false);

    nameScore = cellfun(@(s) ~isempty(regexp(s,'(dff|df[_/]?f|dv_dff|zdff)','once')), namesLower);
    nElts     = cellfun(@numel, {cands.val});                % prefer the biggest
    [~, idx]  = sortrows([double(nameScore(:)) nElts(:)], [-1 -2]);
    pick      = cands(idx(1));

    dff        = pick.val;
    chosenPath = pick.path;    % char
end

function list = collect_numeric_2d(S, base)
    list = struct('path',{},'val',{},'size',{});
    fn = fieldnames(S);
    for i = 1:numel(fn)
        f = fn{i};
        p = f; if strlength(base)>0, p = base + "." + f; end
        val = S.(f);
        if isnumeric(val) && ismatrix(val) && ~isscalar(val)
            list(end+1) = struct('path',string(p),'val',val,'size',size(val)); %#ok<AGROW>
        elseif isstruct(val)
            sub = collect_numeric_2d(val, string(p));
            list = [list sub]; %#ok<AGROW>
        end
    end
end

function [ok, val] = get_by_path(S, pathStr)
    ok = false; val = [];
    if isempty(pathStr), return; end
    parts = strsplit(pathStr,'.');
    cur = S;
    for i = 1:numel(parts)
        if isfield(cur, parts{i})
            cur = cur.(parts{i});
        else
            return;
        end
    end
    ok = true; val = cur;
end
