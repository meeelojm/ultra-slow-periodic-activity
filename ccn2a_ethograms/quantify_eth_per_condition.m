%% ===========================================
%  Ethogram: High vs Mid vs Low (f010, f011, f012)
%  Variables: choose any fields present in ethogram_*.mat
%  ===========================================

clear; clc;

% -------- Base directories to search (first existing wins) --------
base_dirs = { ...
  "\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low", ...
  "Z:\temp\Emilian temp\ccn2a\high2low", ...
  "Z:\Emilian temp\ccn2a\high2low" ...
};

% -------- Animals & Conditions --------
fish_ids = {'f010'};
conds    = {'high','mid','low'};     % folder names ('highact' etc. also supported)

% -------- Variables to analyze (edit to taste) --------
% mode = 'ang'  -> take |d/dt| before metrics (tail angle)
% mode = 'raw'  -> use the signal as-is (vir speed, eye/head/mouth/operculum signals)
vars_to_do = [ ...
   struct('field','tai_ang_uni_frames', 'mode','ang', 'label','tail angle'), ...
   struct('field','vir_spe_uni_frames', 'mode','raw', 'label','virtual speed'), ...
   struct('field','eye_sig',            'mode','raw', 'label','eye signal'), ...
   struct('field','hea_sig',            'mode','raw', 'label','head signal'), ...
   struct('field','mou_sig',            'mode','raw', 'label','mouth signal'), ...
   struct('field','ope_sig',            'mode','raw', 'label','operculum signal') ...
];

% -------- Analysis parameters --------
params = struct;
params.smooth_sec     = 0.10;   % smoothing before derivative
params.thr_k          = 4;      % activity threshold = median + k*MAD
params.min_bout_sec   = 0.20;   % minimum bout duration
params.alpha          = 0.05;   % star cutoff
params.threshold_mode = 'pooled'; % use pooled per-fish threshold across stages (recommended)

%% ========== Load all files & build per-variable matrices ==========
DATA = struct();   % DATA.(field).(cond) = [T x N] matrix; DATA.t_common.(field) = t
for v = 1:numel(vars_to_do)
    vf   = vars_to_do(v).field;
    mode = vars_to_do(v).mode;

    % --- Load & keep per fish per condition time+signal ---
    tCell  = cell(numel(fish_ids), numel(conds));
    yCell  = cell(numel(fish_ids), numel(conds));
    for fi = 1:numel(fish_ids)
        for ci = 1:numel(conds)
            fpath = find_etho_path(base_dirs, fish_ids{fi}, conds{ci});
            [t, y] = load_etho_var(fpath, vf);
            tCell{fi,ci} = t(:);
            yCell{fi,ci} = y(:);
        end
    end

    % --- Build common time base (intersection across all fish & conds) ---
    t0 = -inf; t1 = inf; dts = [];
    for fi = 1:numel(fish_ids)
        for ci = 1:numel(conds)
            ti = tCell{fi,ci};
            if isempty(ti), error('Empty time vector in %s/%s for %s',fish_ids{fi},conds{ci},vf); end
            t0 = max(t0, ti(1));
            t1 = min(t1, ti(end));
            dts(end+1) = median(diff(ti));
        end
    end
    assert(t1>t0, 'No overlapping time across fish/conditions for %s', vf);
    dt_common = median(dts);
    t_common  = (t0:dt_common:t1).';

    % --- Resample to common time & assemble matrices (columns=fish) ---
    for ci = 1:numel(conds)
        M = nan(numel(t_common), numel(fish_ids));
        for fi = 1:numel(fish_ids)
            M(:,fi) = resamp_to(tCell{fi,ci}, yCell{fi,ci}, t_common);
        end
        DATA.(vf).(conds{ci}) = M;
    end
    DATA.t_common.(vf) = t_common;
    DATA.mode.(vf)     = mode;
end

%% ========== Quantify & plot for each variable ==========
OUT = struct();
for v = 1:numel(vars_to_do)
    vf   = vars_to_do(v).field;
    mode = vars_to_do(v).mode;
    t    = DATA.t_common.(vf);
    H    = DATA.(vf).high;
    M    = DATA.(vf).mid;
    L    = DATA.(vf).low;

    OUT.(vf) = quantify_by_condition(t, H, M, L, params, mode, vars_to_do(v).label);
end

%% ===================== Helpers ==========================
function fpath = find_etho_path(base_dirs, fish, cond)
    % Try 'condact' then 'cond', pick the newest ethogram_*.mat
    subs = {cond+"act", cond};
    for b = 1:numel(base_dirs)
        for s = 1:numel(subs)
            pat = fullfile(base_dirs{b}, fish, subs{s}, "suite2p", "ethogram_*.mat");
            d = dir(pat);
            if ~isempty(d)
                [~,ix] = max([d.datenum]);
                fpath = string(fullfile(d(ix).folder, d(ix).name));
                return;
            end
        end
    end
    error('No ethogram_*.mat for %s/%s', fish, cond);
end

function [t, y] = load_etho_var(fpath, field)
% Robust loader for ethogram variables.
% Returns a time vector t that matches the signal length.
    S = load(fpath);                         % load everything (some files miss fra_rat)
    assert(isfield(S, field), 'Field %s not found in %s', field, fpath);

    y = S.(field)(:);
    n = numel(y);

    % ---- preferred companion time names that must match length n ----
    pref = {};
    if endsWith(field, '_uni_frames')
        % try field-specific time vectors first (if present in your files)
        pref = {'fra_tim_uni_frames','tim_uni_frames'};
    end
    % generic candidates (will only accept if length matches n)
    cand = [pref, {'fra_tim_uni','fra_tim','time','t'}];

    % 1) try to find a time vector in the file with the SAME length as y
    for k = 1:numel(cand)
        nm = cand{k};
        if isfield(S, nm) && numel(S.(nm)) == n
            t = S.(nm)(:);
            return
        end
    end

    % 2) otherwise, build t from fps and signal length
    fs = [];
    if isfield(S,'fra_rat_get') && ~isempty(S.fra_rat_get), fs = double(S.fra_rat_get); end
    if isempty(fs) && isfield(S,'fra_rat') && ~isempty(S.fra_rat), fs = double(S.fra_rat); end
    if isempty(fs)
        % last resort: try to infer fs from any time vector present
        an = {'fra_tim_uni','fra_tim','time','t'};
        for k = 1:numel(an)
            if isfield(S,an{k}) && numel(S.(an{k}))>1
                tt = S.(an{k})(:);
                dt = median(diff(tt(~isnan(tt))));
                if isfinite(dt) && dt>0, fs = 1/dt; break; end
            end
        end
    end
    assert(~isempty(fs) && isfinite(fs) && fs>0, ...
        'No fps/time information found to build t for %s', fpath);

    t = (0:n-1)'/fs;                          % matches y by construction
end


function yq = resamp_to(t, y, t_common)
    [tu, ia] = unique(t, 'stable');
    yu = y(ia);
    yq = interp1(tu, yu, t_common, 'linear', 'extrap');  % extrap OK: we already clip to intersection
end


function OUT = quantify_by_condition(t, high, mid, low, params, mode, nice_label)
    if nargin < 7, nice_label = ''; end
    fs    = 1/median(diff(t));
    winSm = max(1, round(params.smooth_sec*fs));

    conds = {'High','Mid','Low'};
    X{1} = high; X{2} = mid; X{3} = low;
    [T,N] = size(high);
    assert(all(cellfun(@(A) all(size(A)==[T N]), X)), 'Matrix sizes mismatch');

    % Build “speed-like” signals to threshold:
    SPD = cell(3,1);
    for c = 1:3
        x = X{c};
        if strcmpi(mode,'ang')           % tail angle -> |d/dt|
            x  = movmean(x, winSm, 1, 'omitnan');
            dt = 1/fs;
            [dx,~] = gradient(x, dt, 1);
            SPD{c} = abs(dx);
        else                              % raw intensity / speed
            SPD{c} = movmean(abs(x), winSm, 1, 'omitnan');
        end
    end

    % Per-fish thresholds (pooled across stages = recommended)
    thr = nan(1,N);
    if strcmpi(params.threshold_mode,'pooled')
        for i=1:N
            pool = [SPD{1}(:,i); SPD{2}(:,i); SPD{3}(:,i)];
            medS = median(pool, 'omitnan');
            madS = mad(pool, 1, 1);
            madS = max(madS, 1e-10);
            thr(i) = medS + params.thr_k*madS;
        end
    end

    % Metrics per condition
    M = struct;
    for c = 1:3
        spd = SPD{c};
        % make activity mask
        if strcmpi(params.threshold_mode,'pooled')
            th = thr;
        else
            medS = median(spd,1,'omitnan'); madS = mad(spd,1,1); madS = max(madS,1e-10);
            th = medS + params.thr_k*madS;
        end
        act = spd > th;
        minBout = round(params.min_bout_sec*fs);
        for i=1:N, act(:,i) = bwareaopen(act(:,i), minBout); end

        % Basic metrics (all on spd)
        M(c).name              = conds{c};
        M(c).mean_signal       = colfun(@(z) mean(z, 'omitnan'), spd);
        M(c).p95_signal        = colfun(@(z) prctile(z,95), spd);
        M(c).active_frac       = colfun(@(z) mean(z, 'omitnan'), act);
        M(c).bout_rate_per_min = bout_rate(act, fs);
        M(c).integral          = colfun(@(z) trapz(t,z), spd);  % total “vigor”
    end

    % Assemble per-fish tables
    OUT.per_fish = struct();
    fields = fieldnames(M(1)); fields(strcmp(fields,'name')) = [];
    for f = 1:numel(fields)
        fn = fields{f};
        Tab = table((1:N).','VariableNames',{'fish'});
        for c=1:3
            Tab.(sprintf('%s_%s', fn, lower(M(c).name))) = M(c).(fn);
        end
        OUT.per_fish.(fn) = Tab;
    end

    % Plot/stats for key metrics
    key = {'mean_signal','p95_signal','active_frac','bout_rate_per_min','integral'};
    for i=1:numel(key)
        fn = key{i};
        A  = [M(1).(fn), M(2).(fn), M(3).(fn)];     % N x 3
        [pF, pHolm, med3] = rm_stats(A, conds);
        OUT.stats.(fn) = struct('friedman_p',pF,'holm_pairwise_p',pHolm,'medians',med3);
        ttl = sprintf('%s (%s)', strrep(fn,'_',' '), nice_label);
        paired_plot(A, conds, ttl, pF, med3, pHolm, params.alpha);
    end
end

function v = colfun(f, X)
    v = nan(size(X,2),1);
    for i=1:size(X,2), xi = X(:,i); v(i) = f(xi); end
end

function r = bout_rate(act, fs)
    N = size(act,2); r = nan(N,1);
    for i=1:N
        lab = bwlabel(act(:,i));
        nB  = max(lab);
        r(i) = nB / (numel(act(:,i))/fs/60);  % bouts per minute
    end
end

function [pF, pHolm, med3] = rm_stats(A, conds)
    [pF,~,~] = friedman(A,1,'off');
    pairs = [1 2; 1 3; 2 3]; % H-M, H-L, M-L
    p = nan(3,1);
    for k=1:3
        [p(k),~] = signrank(A(:,pairs(k,1)), A(:,pairs(k,2)));
    end
    [ps,ord] = sort(p); m = numel(p);
    holm = ps .* (m - (1:m)' + 1);
    adj  = holm; for i=2:m, adj(i)=max(adj(i-1), holm(i)); end; adj=min(adj,1);
    pHolm = nan(size(p)); pHolm(ord) = adj;
    med3  = median(A,1,'omitnan');
    fprintf('%s: Friedman p=%.3g | medians %s=%.3g, %s=%.3g, %s=%.3g | Holm p (H-M,H-L,M-L) = [%g %g %g]\n', ...
        inputname(1), pF, conds{1}, med3(1), conds{2}, med3(2), conds{3}, pHolm(1), pHolm(2), pHolm(3));
end

function paired_plot(A, conds, metric_name, pF, med3, pHolm, alpha)
    if nargin < 7 || isempty(pHolm), pHolm = [NaN NaN NaN]; end
    if nargin < 8 || isempty(alpha), alpha = 0.05; end

    figure('Color','w','Name',metric_name); hold on; grid on; box on;
    N = size(A,1); x = 1:3;
    cols = [0.20 0.70 0.20; 1.00 0.55 0.10; 0.00 0.45 1.00]; % High, Mid, Low

    % IQR backdrop + medians
    for g=1:3
        q = quantile(A(:,g), [0.25 0.5 0.75]);
        w = 0.24;
        patch([x(g)-w x(g)+w x(g)+w x(g)-w], [q(1) q(1) q(3) q(3)], cols(g,:), ...
              'FaceAlpha',0.15,'EdgeColor',[.6 .6 .6]);
        plot([x(g)-w x(g)+w],[q(2) q(2)],'k-','LineWidth',2);
    end

    % Paired lines + points
    for i=1:N, plot(x, A(i,:), '-', 'Color', [.7 .7 .7]); end
    for g=1:3, scatter(g*ones(N,1), A(:,g), 40, 'filled', ...
            'MarkerFaceColor', cols(g,:), 'MarkerEdgeColor','k'); end

    xlim([0.5 3.5]); set(gca,'XTick',x,'XTickLabel',conds);
    ylabel(metric_name);
    title(sprintf('Friedman p=%.3g | medians: %s=%.3g, %s=%.3g, %s=%.3g', ...
        pF, conds{1}, med3(1), conds{2}, med3(2), conds{3}));

    % Significance bars (Holm-corrected)
    pairs = [1 2; 1 3; 2 3]; % H-M, H-L, M-L
    allY = A(:); allY = allY(isfinite(allY));
    if isempty(allY), allY = 0; end
    yMax  = max(allY); yMin = min(allY); ySpan = max(yMax-yMin, eps);
    base  = yMax + 0.06*ySpan; step = 0.08*ySpan;

    topY = [];
    for k=1:3
        p = pHolm(k);
        if ~isnan(p) && p < alpha
            x1 = pairs(k,1); x2 = pairs(k,2);
            y  = base + (numel(topY))*step;
            plot([x1 x1 x2 x2], [y y+0.015*ySpan y+0.015*ySpan y], 'k-', 'LineWidth',1.2);
            text((x1+x2)/2, y+0.02*ySpan, pstars(p), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom','FontWeight','bold','Color','k');
            topY(end+1) = y + 0.08*ySpan; %#ok<AGROW>
        end
    end
    if ~isempty(topY), ylim([yMin-0.05*ySpan, max(topY)+0.10*ySpan]); end
end

function s = pstars(p)
    if p < 1e-3, s = '***';
    elseif p < 1e-2, s = '**';
    elseif p < 5e-2, s = '*';
    else, s = 'n.s.'; % not drawn unless you change the logic
    end
end
