%% ======================= FIGURE 6: Trial-to-trial variability (REAL vs RANDOM; SIG-ONLY) =======================
baseFolder = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';
resultsMat = 'dffs_repact_respcells_randcontrol.mat';

% ---------- global config ----------
minTrialsOK = 3;                % >= this many valid trials to score a neuron
useSpearman = false;            % Pearson (false) or Spearman (true)
SAVE_PDF    = true;

% significance settings
alpha       = 0.05;             % significance level
mccMethod   = 'fdr';            % 'fdr' (BH), 'bonferroni', or 'none'
minSigPairs = 1;                % require at least this many sig pairs to compute rMean

% Scenarios to iterate: label + trial indices (1-based in the 2nd dim)
scenarios = struct( ...
  'name', {'Light','Tap','Light+Tap'}, ...
  'short',{'light','tap','lt'}, ...
  'trials',{1:8, 9:16, 17:24} );

fprintf('\n[Pipeline] Scanning for "%s" under:\n  %s\n', resultsMat, baseFolder);

% find dataset folders
allFolders = strsplit(genpath(baseFolder), pathsep);
allFolders = allFolders(~cellfun(@isempty, allFolders));
validFolders = {};
for i = 1:numel(allFolders)
    if exist(fullfile(allFolders{i}, resultsMat), 'file') == 2
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
    end
end
fprintf('[Pipeline] Found %d dataset folder(s).\n', numel(validFolders));

% ---------- POPULATION ACCUMULATORS (REAL vs RANDOM; per scenario) ----------
% All valid neurons:
popReal = struct('light',[],'tap',[],'lt',[]);
popRand = struct('light',[],'tap',[],'lt',[]);
popRealFish = struct('light',strings(0,1),'tap',strings(0,1),'lt',strings(0,1));
popRandFish = struct('light',strings(0,1),'tap',strings(0,1),'lt',strings(0,1));
% Top-10% per fish:
popTopReal = struct('light',[],'tap',[],'lt',[]);
popTopRand = struct('light',[],'tap',[],'lt',[]);
popTopRealFish = struct('light',strings(0,1),'tap',strings(0,1),'lt',strings(0,1));
popTopRandFish = struct('light',strings(0,1),'tap',strings(0,1),'lt',strings(0,1));

cell_sep = repmat('-',1,70);

% ========================= PER-FOLDER PROCESSING =========================
for f = 1:numel(validFolders)
    dataDir   = validFolders{f};
    resultsFp = fullfile(dataDir, resultsMat);
    [~, fishLabel] = fileparts(dataDir);   % use folder name as fish label

    fprintf('\n%s\n[Dataset %d/%d] %s\n', cell_sep, f, numel(validFolders), dataDir);

    S = load(resultsFp);
    hasReal = isfield(S,'dff_fra_trial_cell');
    hasRand = isfield(S,'dff_fra_trial_cell_r');

    if ~(hasReal || hasRand)
        warning('  Skipping: neither "dff_fra_trial_cell" nor "dff_fra_trial_cell_r" found.');
        continue
    end

    if hasReal
        dff_real = S.dff_fra_trial_cell;       % [frames x trials x neurons]
        if isempty(dff_real) || ndims(dff_real) ~= 3
            warning('  "dff_fra_trial_cell" has unexpected size. Ignoring REAL for this fish.');
            hasReal = false;
        end
    end

    if hasRand
        dff_rand = S.dff_fra_trial_cell_r;     % [frames x trials x neurons]
        if isempty(dff_rand) || ndims(dff_rand) ~= 3
            warning('  "dff_fra_trial_cell_r" has unexpected size. Ignoring RANDOM for this fish.');
            hasRand = false;
        end
    end

    % ---------- process EACH SCENARIO for this fish ----------
    for s = 1:numel(scenarios)
        label = scenarios(s).name;    % "Light", "Tap", "Light+Tap"
        tag   = scenarios(s).short;   % "light","tap","lt"
        trIdx = scenarios(s).trials;

        % ===== REAL =====
        if hasReal
            rMean_real = sig_rmean_per_neuron(dff_real, trIdx, minTrialsOK, useSpearman, alpha, mccMethod, minSigPairs);
            vr = rMean_real(~isnan(rMean_real));
            % all-valid accumulation
            popReal.(tag)     = [popReal.(tag); vr]; %#ok<AGROW>
            popRealFish.(tag) = [popRealFish.(tag); repmat(string(fishLabel), numel(vr),1)]; %#ok<AGROW>
            % top-10% per fish
            if ~isempty(vr)
                K = max(1, round(0.10 * numel(vr)));
                [~, ord] = sort(vr, 'descend');
                topVals = vr(ord(1:K));
                popTopReal.(tag)     = [popTopReal.(tag); topVals]; %#ok<AGROW>
                popTopRealFish.(tag) = [popTopRealFish.(tag); repmat(string(fishLabel), numel(topVals),1)]; %#ok<AGROW>
            end
            fprintf('  [%s | REAL] collected %d neurons (top10%%=%d)\n', label, numel(vr), ...
                    ~isempty(vr) * max(1, round(0.10*numel(vr))));
        end

        % ===== RANDOM =====
        if hasRand
            rMean_rand = sig_rmean_per_neuron(dff_rand, trIdx, minTrialsOK, useSpearman, alpha, mccMethod, minSigPairs);
            vz = rMean_rand(~isnan(rMean_rand));
            % all-valid accumulation
            popRand.(tag)     = [popRand.(tag); vz]; %#ok<AGROW>
            popRandFish.(tag) = [popRandFish.(tag); repmat(string(fishLabel), numel(vz),1)]; %#ok<AGROW>
            % top-10% per fish
            if ~isempty(vz)
                K = max(1, round(0.10 * numel(vz)));
                [~, ord] = sort(vz, 'descend');
                topVals = vz(ord(1:K));
                popTopRand.(tag)     = [popTopRand.(tag); topVals]; %#ok<AGROW>
                popTopRandFish.(tag) = [popTopRandFish.(tag); repmat(string(fishLabel), numel(topVals),1)]; %#ok<AGROW>
            end
            fprintf('  [%s | RANDOM] collected %d neurons (top10%%=%d)\n', label, numel(vz), ...
                    ~isempty(vz) * max(1, round(0.10*numel(vz))));
        end
    end
end

fprintf('\n[Pipeline] Finished per-fish processing. Now making POPULATION plots…\n');

%% ======================= POPULATION: REAL vs RANDOM (3 scenarios; all neurons; SIG-ONLY) =======================
outPop = fullfile(baseFolder, 'fig6_population_real_vs_random_SIG');
if ~exist(outPop,'dir'), mkdir(outPop); end

make_real_vs_random_plot(popReal.light, popRand.light, 'Light',    fullfile(outPop,'light_vs_random.png'),    SAVE_PDF);
make_real_vs_random_plot(popReal.tap,   popRand.tap,   'Tap',      fullfile(outPop,'tap_vs_random.png'),      SAVE_PDF);
make_real_vs_random_plot(popReal.lt,    popRand.lt,    'Light+Tap',fullfile(outPop,'lt_vs_random.png'),       SAVE_PDF);
%% ======================= POPULATION: REAL vs RANDOM (Top-10% per fish; SIG-ONLY) =======================
% Per-scenario figures
make_real_vs_random_plot(popTopReal.light, popTopRand.light, 'Light (Top10%/fish)',    fullfile(outPop,'light_vs_random_TOP10.png'),    SAVE_PDF);
make_real_vs_random_plot(popTopReal.tap,   popTopRand.tap,   'Tap (Top10%/fish)',      fullfile(outPop,'tap_vs_random_TOP10.png'),      SAVE_PDF);
make_real_vs_random_plot(popTopReal.lt,    popTopRand.lt,    'Light+Tap (Top10%/fish)',fullfile(outPop,'lt_vs_random_TOP10.png'),       SAVE_PDF);

% Combined 3×2 panel (side-by-side groups) — CONSISTENT COLORS + WIDER SPACING
figAll = figure('Color','w','Name','Population rMean_{sig}: REAL vs RANDOM (3 scenarios)');
hold on; grid on; box on

C.real   = [0 0.4470 0.7410];
C.random = [0.8500 0.3250 0.0980];
baseX = [1 3 5];  d = 0.25;

[hR1,hZ1] = plot_group(baseX(1)-d, popReal.light,  C.real,   'real');
[~,  hZ2] = plot_group(baseX(1)+d, popRand.light,  C.random, 'random');
[hR2,~]   = plot_group(baseX(2)-d, popReal.tap,    C.real,   'real');
[~,  ~]   = plot_group(baseX(2)+d, popRand.tap,    C.random, 'random');
[hR3,~]   = plot_group(baseX(3)-d, popReal.lt,     C.real,   'real');
[~,  ~]   = plot_group(baseX(3)+d, popRand.lt,     C.random, 'random');

xline(2,'k:','LineWidth',0.5); xline(4,'k:','LineWidth',0.5); % separators
xlim([0 6]); ylim([-1 1]);
xticks(baseX); xticklabels({'Light','Tap','Light+Tap'});
ylabel('rMean_{sig} (trial-to-trial)');
title('Population reliability (sig-only): REAL vs RANDOM');
legend([hR1 hZ1], {'real','random'}, 'Location','southoutside','Orientation','horizontal');
exportgraphics(figAll, fullfile(outPop,'real_vs_random_3scenarios.png'), 'Resolution', 200);
if SAVE_PDF, exportgraphics(figAll, fullfile(outPop,'real_vs_random_3scenarios.pdf'), 'ContentType','vector'); end

% Combined 3-scenario TOP10% panel — CONSISTENT COLORS + WIDER SPACING
figAllTop = figure('Color','w','Name','Population rMean_{sig} TOP10: REAL vs RANDOM (3 scenarios)');
hold on; grid on; box on

C.real   = [0 0.4470 0.7410];
C.random = [0.8500 0.3250 0.0980];
baseX = [1 3 5];  d = 0.25;

[hTR1,hTZ1] = plot_group(baseX(1)-d, popTopReal.light,  C.real,   'real');
[~,    ~]   = plot_group(baseX(1)+d, popTopRand.light,  C.random, 'random');
[hTR2,~]    = plot_group(baseX(2)-d, popTopReal.tap,    C.real,   'real');
[~,    ~]   = plot_group(baseX(2)+d, popTopRand.tap,    C.random, 'random');
[hTR3,~]    = plot_group(baseX(3)-d, popTopReal.lt,     C.real,   'real');
[~,    ~]   = plot_group(baseX(3)+d, popTopRand.lt,     C.random, 'random');

xline(2,'k:','LineWidth',0.5); xline(4,'k:','LineWidth',0.5);
xlim([0 6]); ylim([-1 1]);
xticks(baseX); xticklabels({'Light','Tap','Light+Tap'});
ylabel('rMean_{sig} (Top10% per fish)');
title('Population TOP10 (sig-only): REAL vs RANDOM');
legend([hTR1 hTZ1], {'real','random'}, 'Location','southoutside','Orientation','horizontal');
exportgraphics(figAllTop, fullfile(outPop,'real_vs_random_3scenarios_TOP10.png'), 'Resolution', 200);
if SAVE_PDF
    exportgraphics(figAllTop, fullfile(outPop,'real_vs_random_3scenarios_TOP10.pdf'), 'ContentType','vector');
end

% Long-form CSV for TOP10%
scenTop  = [repmat("light",numel(popTopReal.light),1); repmat("tap",numel(popTopReal.tap),1); repmat("lt",numel(popTopReal.lt),1); ...
            repmat("light",numel(popTopRand.light),1); repmat("tap",numel(popTopRand.tap),1); repmat("lt",numel(popTopRand.lt),1)];
trainTop = [repmat("real", numel(popTopReal.light)+numel(popTopReal.tap)+numel(popTopReal.lt),1); ...
            repmat("random", numel(popTopRand.light)+numel(popTopRand.tap)+numel(popTopRand.lt),1)];
rTopAll  = [popTopReal.light; popTopReal.tap; popTopReal.lt; popTopRand.light; popTopRand.tap; popTopRand.lt];

TpopTop = table(scenTop, trainTop, rTopAll, 'VariableNames', {'scenario','train','rMean_sig_top10'});
writetable(TpopTop, fullfile(outPop,'population_rMeanSIG_long_real_vs_random_TOP10.csv'));

%% ============================== HELPERS (put these at the very end) ==============================
function rMean = sig_rmean_per_neuron(dff_trials, trIdx, minTrialsOK, useSpearman, alpha, mccMethod, minSigPairs)
% dff_trials: [frames x trials x neurons]
% trIdx     : trial indices to use (1-based)
% returns   : rMean (sig-only) per neuron (NaN if not enough sig pairs)
    if isempty(dff_trials) || ndims(dff_trials) ~= 3
        rMean = nan(0,1); return
    end
    nAvailTrials = size(dff_trials,2);
    trIdx = trIdx(trIdx <= nAvailTrials);            % clamp to available trials
    nNeurons = size(dff_trials,3);
    rMean    = nan(nNeurons,1);

    for n = 1:nNeurons
        T = squeeze(dff_trials(:,trIdx,n)).';        % [trials x frames]
        if isempty(T), continue; end
        good = all(isfinite(T),2) & std(T,0,2) > eps; % non-NaN & non-flat
        T    = T(good,:);
        nT   = size(T,1);
        if nT < minTrialsOK, continue; end

        if useSpearman
            [R,P] = corr(T','Type','Spearman','rows','pairwise');
        else
            [R,P] = corr(T','rows','pairwise');      % Pearson (two-tailed)
        end

        msk   = triu(true(nT),1);
        vR    = R(msk);
        vP    = P(msk);

        sigVec = mcc_mask_vec(vP, alpha, mccMethod);
        if sum(sigVec) >= minSigPairs
            vRsig    = max(min(vR(sigVec), 0.999999), -0.999999);
            rMean(n) = tanh(mean(atanh(vRsig), 'omitnan'));   % Fisher-mean on sig pairs
        else
            rMean(n) = NaN;
        end
    end
end

function sig = mcc_mask_vec(p, alpha, method)
% p: vector (or column) of p-values; returns logical vector of significant entries.
    p = p(:);
    switch lower(method)
        case 'none'
            sig = p < alpha;
        case 'bonferroni'
            sig = p < (alpha / max(numel(p),1));
        case 'fdr'  % Benjamini–Hochberg
            sig = fdr_bh_mask(p, alpha);
        otherwise
            error('Unknown mcc method: %s', method);
    end
end

function sig = fdr_bh_mask(p, alpha)
% Benjamini–Hochberg FDR (two-sided). Returns logical vector (same size as p).
    [ps, idx] = sort(p(:));
    m = numel(ps);
    if m == 0
        sig = false(0,1); return
    end
    thresh = (1:m)'/m * alpha;
    pass = ps <= thresh;
    k = find(pass, 1, 'last');
    keep = false(m,1);
    if ~isempty(k), keep(1:k) = true; end
    sig = false(m,1);
    sig(idx) = keep;
end

function make_real_vs_random_plot(realVals, randVals, label, outPng, SAVE_PDF)
    C.real   = [0 0.4470 0.7410];     % blue
    C.random = [0.8500 0.3250 0.0980];% orange

    fig = figure('Color','w','Name',sprintf('%s: REAL vs RANDOM (rMean_{sig})', label));
    hold on; grid on; box on
    d = 0.20;

    nR = numel(realVals); nZ = numel(randVals);
    if ~isempty(realVals)
        hR = boxchart(ones(nR,1)*(1-d), realVals, ...
            'BoxFaceAlpha',0.25, 'BoxFaceColor',C.real, 'MarkerStyle','none', ...
            'DisplayName','real');
        scatter((1-d)+0.045*randn(nR,1), realVals, 10, 'filled', ...
            'MarkerFaceColor',C.real, 'MarkerFaceAlpha',0.35, 'MarkerEdgeColor','none');
    end
    if ~isempty(randVals)
        hZ = boxchart(ones(nZ,1)*(1+d), randVals, ...
            'BoxFaceAlpha',0.25, 'BoxFaceColor',C.random, 'MarkerStyle','none', ...
            'DisplayName','random');
        scatter((1+d)+0.045*randn(nZ,1), randVals, 10, 'filled', ...
            'MarkerFaceColor',C.random, 'MarkerFaceAlpha',0.35, 'MarkerEdgeColor','none');
    end

    xlim([0.5 1.5]); ylim([-1 1]); xticks(1); xticklabels({label});
    ylabel('rMean_{sig} (trial-to-trial)');
    title(sprintf('%s: REAL (n=%d) vs RANDOM (n=%d)', label, nR, nZ));
    if exist('hR','var') && exist('hZ','var')
        legend([hR hZ], {'real','random'}, 'Location','southoutside','Orientation','horizontal');
    end
    exportgraphics(fig, outPng, 'Resolution', 200);
    if SAVE_PDF, exportgraphics(fig, strrep(outPng,'.png','.pdf'), 'ContentType','vector'); end
    close(fig);
end


function [hBox,hSc] = plot_group(xpos, vals, color, label)
% Draw one group at xpos with a given color; returns handles for legend
    if isempty(vals), hBox = gobjects(1); hSc = gobjects(1); return; end
    hBox = boxchart(ones(numel(vals),1)*xpos, vals, ...
        'BoxFaceAlpha',0.25, 'BoxFaceColor',color, 'MarkerStyle','none', ...
        'DisplayName', label);
    hSc = scatter(xpos + 0.045*randn(numel(vals),1), vals, 10, 'filled', ...
        'MarkerFaceColor',color, 'MarkerFaceAlpha',0.35, 'MarkerEdgeColor','none', ...
        'DisplayName', label);
end



function topVals = get_top_vals(rMean, topFrac)
% Return the top-K (fraction) rMean values for one fish in one scenario
    valid = find(~isnan(rMean));
    if isempty(valid)
        topVals = [];
        return
    end
    K = max(1, round(topFrac * numel(valid)));
    [~, ord] = sort(rMean(valid), 'descend');
    topVals  = rMean(valid(ord(1:K)));
end
