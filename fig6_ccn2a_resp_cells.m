%% ======================= FIGURE 6: Trial-to-trial variability (3 scenarios) =======================
baseFolder = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';
resultsMat = 'dffs_repact_respcells_randcontrol.mat';

% ---------- global config ----------
minTrialsOK = 3;          % >= this many valid trials to score a neuron
useSpearman = false;      % Pearson (false) or Spearman (true)
topFrac     = 0.10;       % keep top 10% per fish (per scenario)
CLIMS       = [-0.5 0.5]; % heatmap color limits
SAVE_PDF    = true;
alpha        = 0.05;          % significance level
mccMethod    = 'fdr';         % 'fdr' (BH), 'bonferroni', or 'none'
minSigPairs  = 1;             % require at least this many sig pairs to compute rMean


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

% ---------- POPULATION ACCUMULATORS (one per scenario) ----------
pop = struct('light',[],'tap',[],'lt',[]);
popFish = struct('light',strings(0,1),'tap',strings(0,1),'lt',strings(0,1));
% NEW: top-10% accumulators (per fish, per scenario)
popTop     = struct('light',[],'tap',[],'lt',[]);
popTopFish = struct('light',strings(0,1),'tap',strings(0,1),'lt',strings(0,1));

cell_sep = repmat('-',1,70);

% ========================= PER-FOLDER PROCESSING =========================
for f = 1:numel(validFolders)
    dataDir   = validFolders{f};
    resultsFp = fullfile(dataDir, resultsMat);
    [~, fishLabel] = fileparts(dataDir);   % use folder name as fish label

    fprintf('\n%s\n[Dataset %d/%d] %s\n', cell_sep, f, numel(validFolders), dataDir);

    S = load(resultsFp);
    if ~isfield(S,'dff_fra_trial_cell')
        warning('  Skipping: "dff_fra_trial_cell" not found.');
        continue
    end
    dff_fra_trial_cell = S.dff_fra_trial_cell;   % [frames x trials x neurons]
    if isempty(dff_fra_trial_cell) || ndims(dff_fra_trial_cell) ~= 3
        warning('  Skipping: "dff_fra_trial_cell" has unexpected size.');
        continue
    end

    % Optional vertical marker: look for "duration" in the file
    haveDuration = isfield(S,'duration');
    if haveDuration, duration = S.duration; end %#ok<NASGU>
    getXcut = @(nFrames) ( haveDuration * min((S.duration*2)+1, nFrames) );

    nNeurons = size(dff_fra_trial_cell,3);
    nAvailTrials = size(dff_fra_trial_cell,2);

    % ---------- process EACH SCENARIO for this fish ----------
    for s = 1:numel(scenarios)
        label = scenarios(s).name;    % "Light", "Tap", "Light+Tap"
        tag   = scenarios(s).short;   % "light","tap","lt"

        % Trial indices for this scenario (respect actual available trials)
        trIdx = scenarios(s).trials;
        trIdx = trIdx(trIdx <= nAvailTrials);     % clamp if fewer trials exist
        if numel(trIdx) < minTrialsOK
            fprintf('  [%s] Not enough trials (%d). Skipping.\n', label, numel(trIdx));
            continue
        end

        % ---- per-neuron reliability for this scenario ----
        Rstore  = cell(nNeurons,1);
        rMean   = nan(nNeurons,1);
        rMedian = nan(nNeurons,1);
        rLOO    = nan(nNeurons,1);
        nTrials = nan(nNeurons,1);

        for n = 1:nNeurons
            T = squeeze(dff_fra_trial_cell(:,trIdx,n)).';    % [trials x frames]
            good = all(isfinite(T),2) & std(T,0,2) > eps;    % non-NaN & non-flat
            T    = T(good,:);
            nT   = size(T,1);
            nTrials(n) = nT;
            if nT < minTrialsOK, continue; end

            if useSpearman
                R = corr(T','Type','Spearman','rows','pairwise');
            else
                R = corr(T','rows','pairwise');              % Pearson
            end
            Rstore{n} = R;

            msk = triu(true(nT),1);
            v   = R(msk);
            v   = max(min(v,0.999999),-0.999999);            % safe clamp

            rMean(n)   = tanh(mean(atanh(v), 'omitnan'));    % Fisher-mean
            rMedian(n) = median(v, 'omitnan');

            loo = nan(nT,1);
            for iTrial = 1:nT
                mu_others = mean(T(setdiff(1:nT,iTrial),:), 1);
                loo(iTrial) = corr(T(iTrial,:)', mu_others(:), 'rows','pairwise');
            end
            rLOO(n) = mean(loo, 'omitnan');
        end

        % ---- select top 10% (this fish, this scenario) ----
        valid = find(~isnan(rMean));
        if isempty(valid)
            fprintf('  [%s] No valid neurons.\n', label);
            continue
        end
        K = max(1, round(topFrac * numel(valid)));
        [~, ord] = sort(rMean(valid), 'descend');
        topIdx   = valid(ord(1:K));
        fprintf('  [%s] Valid neurons: %d | keeping top %d (%.1f%%). Median r̄(top)=%.3f\n', ...
                label, numel(valid), K, 100*K/numel(valid), median(rMean(topIdx),'omitnan'));

        % ---- outputs (scenario-tagged) ----
        outRoot = fullfile(dataDir, 'fig6_outputs', tag);
        if ~exist(outRoot,'dir'), mkdir(outRoot); end

        % correlation matrices (top-10%)
        ncols = min(4, K); nrows = ceil(K/ncols);
        figCorr = figure('Color','w','Name',sprintf('%s: Top 10%% trial correlations', label));
        tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');
        for k = 1:K
            n = topIdx(k); R = Rstore{n};
            nexttile; imagesc(R, CLIMS); axis square
            title(sprintf('%s | #%d  r̄=%.2f  LOO=%.2f  nT=%d', label, n, rMean(n), rLOO(n), nTrials(n)), 'FontSize',8);
            set(gca, 'XTick',1:nTrials(n), 'YTick',1:nTrials(n));
        end
        try, colormap(coolwarm); catch, colormap(parula); end
        cb = colorbar; cb.Label.String = 'trial-to-trial r';

        exportgraphics(figCorr, fullfile(outRoot, sprintf('top10_corr_mats_%s.png', tag)), 'Resolution', 200);
        if SAVE_PDF
            exportgraphics(figCorr, fullfile(outRoot, sprintf('top10_corr_mats_%s.pdf', tag)), 'ContentType','vector');
        end
        close(figCorr);

        % stacked trials for each top neuron (scenario trials only)
        outDirStacks = fullfile(outRoot, sprintf('top10_stacked_trials_%s', tag));
        if ~exist(outDirStacks,'dir'), mkdir(outDirStacks); end
        pdfPath = fullfile(outDirStacks, sprintf('top10_stacked_trials_%s.pdf', tag));
        if SAVE_PDF && exist(pdfPath,'file'), delete(pdfPath); end

        for kk = 1:numel(topIdx)
            n = topIdx(kk);
            temp = squeeze(dff_fra_trial_cell(:,trIdx,n)).';   % scenario trials only
            good = all(isfinite(temp),2) & std(temp,0,2) > eps;
            temp = temp(good,:);
            nT   = size(temp,1);
            if nT < 1, continue; end
            xcut = getXcut(size(temp,2));

            figS = figure('Color','w','Position',[100 80 700 220+80*nT], ...
                          'Name', sprintf('%s | Neuron %d', label, n));
            tiledlayout(nT,1,'TileSpacing','compact','Padding','compact');
            ax = gobjects(nT,1);
            for tr = 1:nT
                ax(tr) = nexttile;
                plot(temp(tr,:), 'LineWidth', 1); hold on
                if xcut > 0, xline(xcut, '--', 'LineWidth', 1); end
                ylabel(sprintf('T%02d', tr)); box off
                if tr==1
                    title(sprintf('%s | neuron %d  |  r\\_bar=%.2f  |  LOO=%.2f  |  nT=%d', ...
                          label, n, rMean(n), rLOO(n), nT));
                end
                if tr < nT, set(gca,'XTickLabel',[]); end
            end
            linkaxes(ax,'x'); xlabel('Frame');

            pngPath = fullfile(outDirStacks, sprintf('neuron_%04d_stacked_%s.png', n, tag));
            exportgraphics(figS, pngPath, 'Resolution', 200);
            if SAVE_PDF
                exportgraphics(figS, pdfPath, 'ContentType','vector', 'Append', true);
            end
            close(figS);
        end

        % metrics table (per fish, per scenario)
        keepers = false(nNeurons,1); keepers(topIdx) = true;
        Tmetrics = table((1:nNeurons).', nTrials, rMean, rMedian, rLOO, keepers, ...
                         'VariableNames', {'neuron','nTrials','rMean','rMedian','rLOO','top10'});
        writetable(Tmetrics, fullfile(outRoot, sprintf('trial_reliability_metrics_%s.csv', tag)));

        % accumulate population (per scenario)
        rvals = rMean(valid);
        pop.(tag)     = [pop.(tag); rvals];                     %#ok<AGROW>
        popFish.(tag) = [popFish.(tag); repmat(string(fishLabel), numel(rvals),1)]; %#ok<AGROW>
        % NEW: accumulate population (TOP-10% per fish in this scenario)
        topVals = rMean(topIdx);
        popTop.(tag)     = [popTop.(tag); topVals];
        popTopFish.(tag) = [popTopFish.(tag); repmat(string(fishLabel), numel(topVals),1)];

        fprintf('  [%s] Outputs saved in: %s\n', label, outRoot);
    end
end

fprintf('\n[Pipeline] Finished per-fish processing. Now making POPULATION plots…\n');

%% ======================= POPULATION: 3-scenario box + scatter =======================
outPop = fullfile(baseFolder, 'fig6_population_outputs');
if ~exist(outPop,'dir'), mkdir(outPop); end

% Prepare data vectors (can be empty if a scenario had no valid trials anywhere)
L  = pop.light;  TL = numel(L);
T  = pop.tap;    TT = numel(T);
LT = pop.lt;     TLT= numel(LT);

if TL+TT+TLT == 0
    warning('No rMean values collected. Nothing to plot.');
else
    figPop = figure('Color','w','Name','Population rMean (Light vs Tap vs Light+Tap)');
    hold on; grid on; box on

    % Boxcharts (positions 1,2,3)
    if ~isempty(L),  boxchart(ones(TL,1)*1,  L,  'BoxFaceAlpha',0.25,'MarkerStyle','none'); end
    if ~isempty(T),  boxchart(ones(TT,1)*2,  T,  'BoxFaceAlpha',0.25,'MarkerStyle','none'); end
    if ~isempty(LT), boxchart(ones(TLT,1)*3, LT, 'BoxFaceAlpha',0.25,'MarkerStyle','none'); end

    % Jittered scatters
    if ~isempty(L),  scatter(1  + 0.06*randn(TL,1),  L, 10, 'filled','MarkerFaceAlpha',0.35); end
    if ~isempty(T),  scatter(2  + 0.06*randn(TT,1),  T, 10, 'filled','MarkerFaceAlpha',0.35); end
    if ~isempty(LT), scatter(3  + 0.06*randn(TLT,1), LT,10, 'filled','MarkerFaceAlpha',0.35); end

    xlim([0.5 3.5]); ylim([-1 1]);
    xticks([1 2 3]); xticklabels({'Light','Tap','Light+Tap'});
    ylabel('rMean (trial-to-trial)');
    title(sprintf('Population reliability | Light=%d, Tap=%d, LT=%d neurons', TL,TT,TLT));

    exportgraphics(figPop, fullfile(outPop,'population_rMean_box_scatter_3scenarios.png'), 'Resolution', 200);
    if SAVE_PDF
        exportgraphics(figPop, fullfile(outPop,'population_rMean_box_scatter_3scenarios.pdf'), 'ContentType','vector');
    end
    % close(figPop);

    % Long-form CSV (scenario, fish, rMean)
    scen = [repmat("light",TL,1); repmat("tap",TT,1); repmat("lt",TLT,1)];
    fish = [popFish.light; popFish.tap; popFish.lt];
    rAll = [L; T; LT];
    Tpop = table(scen, fish, rAll, 'VariableNames', {'scenario','fish','rMean'});
    writetable(Tpop, fullfile(outPop,'population_rMean_long_3scenarios.csv'));

    fprintf('[Population] Saved 3-scenario plot + CSV in: %s\n', outPop);
end

fprintf('\n[Pipeline] Done.\n');

%% ======================= POPULATION: Top-10% (per fish) =======================
Ltop  = popTop.light;   TLtop  = numel(Ltop);
Ttop  = popTop.tap;     TTtop  = numel(Ttop);
LTtop = popTop.lt;      TLTtop = numel(LTtop);

if TLtop + TTtop + TLTtop == 0
    warning('No top-10%% values collected. Skipping top-10%% population plot.');
else
    figTop = figure('Color','w','Name','Population rMean (Top 10% per fish)');
    hold on; grid on; box on

    % Boxcharts at positions 1,2,3
    if ~isempty(Ltop),  boxchart(ones(TLtop,1)*1,  Ltop,  'BoxFaceAlpha',0.25,'MarkerStyle','none'); end
    if ~isempty(Ttop),  boxchart(ones(TTtop,1)*2,  Ttop,  'BoxFaceAlpha',0.25,'MarkerStyle','none'); end
    if ~isempty(LTtop), boxchart(ones(TLTtop,1)*3, LTtop, 'BoxFaceAlpha',0.25,'MarkerStyle','none'); end

    % Jittered scatters
    if ~isempty(Ltop),  scatter(1  + 0.06*randn(TLtop,1),  Ltop, 10, 'filled','MarkerFaceAlpha',0.35); end
    if ~isempty(Ttop),  scatter(2  + 0.06*randn(TTtop,1),  Ttop, 10, 'filled','MarkerFaceAlpha',0.35); end
    if ~isempty(LTtop), scatter(3  + 0.06*randn(TLTtop,1), LTtop,10, 'filled','MarkerFaceAlpha',0.35); end

    xlim([0.5 3.5]); ylim([-1 1]);
    xticks([1 2 3]); xticklabels({'Light (top10%)','Tap (top10%)','Light+Tap (top10%)'});
    ylabel('rMean (trial-to-trial)');
    title(sprintf('Population (Top 10%%/fish) | Light=%d, Tap=%d, LT=%d neurons', TLtop, TTtop, TLTtop));

    outPop = fullfile(baseFolder, 'fig6_population_outputs');
    if ~exist(outPop,'dir'), mkdir(outPop); end
    exportgraphics(figTop, fullfile(outPop,'population_rMean_box_scatter_3scenarios_TOP10.png'), 'Resolution', 200);
    if SAVE_PDF
        exportgraphics(figTop, fullfile(outPop,'population_rMean_box_scatter_3scenarios_TOP10.pdf'), 'ContentType','vector');
    end
    % close(figTop);

    % Long-form CSV for top-10%
    scenTop = [repmat("light",TLtop,1); repmat("tap",TTtop,1); repmat("lt",TLTtop,1)];
    fishTop = [popTopFish.light; popTopFish.tap; popTopFish.lt];
    rTopAll = [Ltop; Ttop; LTtop];
    TpopTop = table(scenTop, fishTop, rTopAll, ...
                    'VariableNames', {'scenario','fish','rMean_top10'});
    writetable(TpopTop, fullfile(outPop,'population_rMean_long_3scenarios_TOP10.csv'));
end
