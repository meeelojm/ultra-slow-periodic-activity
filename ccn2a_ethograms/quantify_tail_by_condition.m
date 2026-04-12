function OUT = quantify_tail_by_condition(t, high, mid, low, params)
% Quantify tail kinematics across High/Mid/Low conditions.
% t      : Tx1 time (s), uniform (use your fra_tim_uni)
% high   : T x N matrix of tai_ang_uni for High (rows=time, cols=fish)
% mid    : T x N matrix for Mid  (same size as high)
% low    : T x N matrix for Low  (same size as high)
% params : struct (all optional)
%   .fs            sampling rate (Hz). If omitted, computed from t.
%   .smooth_sec    smoothing window for angle (default 0.10 s)
%   .thr_k         activity threshold: median(speed)+k*MAD (default k=4)
%   .min_bout_sec  minimum bout duration (default 0.10 s)
%   .calib_b       3x1 coefficients to map vigor->[b0 b1 b2] to VR speed (optional)
%   .alpha         p-value cutoff for significance bars (default 0.05)
%
% Returns OUT with tables per fish and summary stats/plots.

if nargin < 5, params = struct; end
if ~isfield(params,'fs'),           params.fs = 1/median(diff(t)); end
if ~isfield(params,'smooth_sec'),   params.smooth_sec = 0.10;      end
if ~isfield(params,'thr_k'),        params.thr_k = 4;              end
if ~isfield(params,'min_bout_sec'), params.min_bout_sec = 0.10;    end
if ~isfield(params,'alpha'),        params.alpha = 0.05;           end

fs = params.fs;
winSm = max(1, round(params.smooth_sec*fs));

% ---- all conditions in a loop ----
conds = {'High','Mid','Low'};
X{1} = high; X{2} = mid; X{3} = low;
N = size(high,2); T = numel(t);
assert(all(cellfun(@(A) size(A,1)==T && size(A,2)==N, X)), 'All matrices must be T x N');

M = struct;  % metrics
for c = 1:3
    ang = X{c};
    % smooth angle a bit before deriv
    ang_s = movmean(ang, winSm, 1, 'omitnan');

    % signed angular velocity; unsigned speed (time derivative along rows)
    dt = 1/fs;
    [dAng_dt, ~] = gradient(ang_s, dt, 1);
    vel = dAng_dt;                         % T x N, signed angular velocity
    spd = abs(vel);                        % unsigned speed

    % robust activity threshold per fish
    medS = median(spd,1,'omitnan');
    madS = mad(spd,1,1);                   % median absolute deviation
    thr  = medS + params.thr_k*madS;       % 1 x N

    % clean tiny single-sample blips; enforce min bout duration
    act  = spd > thr;                      % T x N
    minBout = round(params.min_bout_sec*fs);
    for i = 1:size(act,2)
        act(:,i) = bwareaopen(act(:,i), minBout);  % truly per-column
    end

    % metrics
    M(c).name = conds{c};
    M(c).mean_speed        = colfun(@(x) mean(x, 'omitnan'), spd);     % N x 1
    M(c).p95_speed         = colfun(@(x) prctile(x,95), spd);          % N x 1
    M(c).active_frac       = colfun(@(x) mean(x, 'omitnan'), act);     % N x 1
    M(c).bout_rate_per_min = bout_rate(act, fs);                       % N x 1
    M(c).vigor_integral    = colfun(@(x) trapz(t,x), spd);             % N x 1

    % beat rate from peaks in angle (avoid double counting)
    M(c).beat_rate_hz      = beat_rate_from_angle(ang_s, t);

    % optional: predicted virtual distance from tail (needs calib)
    if isfield(params,'calib_b') && numel(params.calib_b)==3
        vigor = movmean(spd, winSm, 1, 'omitnan');     % smoother drive
        yhat  = params.calib_b(1) + params.calib_b(2)*vigor + params.calib_b(3)*(vigor.^2);
        yhat  = max(yhat,0);                           % clip negatives
        M(c).pred_dist = colfun(@(x) trapz(t,x), yhat);
    end
end

% ---- assemble tables ----
metNames = fieldnames(M(1)); metNames = metNames(~strcmp(metNames,'name'));
OUT.per_fish = struct();
for k = 1:numel(metNames)
    fn = metNames{k};
    Ttab = table((1:N).', 'VariableNames', {'fish'});
    for c = 1:3
        v = M(c).(fn); Ttab.(sprintf('%s_%s', fn, lower(M(c).name))) = v;
    end
    OUT.per_fish.(fn) = Ttab;
end

% ---- stats & plots for selected metrics ----
key = {'mean_speed','p95_speed','beat_rate_hz','active_frac','bout_rate_per_min','vigor_integral'};
for i = 1:numel(key)
    fn = key{i};
    A = [M(1).(fn), M(2).(fn), M(3).(fn)];   % N x 3
    [pF, pHolm, med3] = rm_stats(A, conds);
    OUT.stats.(fn) = struct('friedman_p',pF,'holm_pairwise_p',pHolm,'medians',med3);
    paired_plot(A, conds, fn, pF, med3, pHolm, params.alpha);
end

% optional: predicted distance plot
if isfield(M(1),'pred_dist')
    A = [M(1).pred_dist, M(2).pred_dist, M(3).pred_dist];
    [pF, pHolm, med3] = rm_stats(A, conds);
    OUT.stats.pred_dist = struct('friedman_p',pF,'holm_pairwise_p',pHolm,'medians',med3);
    paired_plot(A, conds, 'pred_dist', pF, med3, pHolm, params.alpha);
end
end

% ---------- helpers ----------

function v = colfun(f, X)
    v = nan(size(X,2),1);
    for i=1:size(X,2), xi = X(:,i); v(i) = f(xi); end
end

function r = bout_rate(act, fs)
    % bouts per minute based on boolean activity (T x N)
    N = size(act,2); r = nan(N,1);
    for i=1:N
        a = act(:,i);
        lab = bwlabel(a);
        nB  = max(lab);
        durMin = numel(a)/fs/60;
        r(i) = nB / max(durMin, eps);
    end
end

function fr = beat_rate_from_angle(ang, t)
    fs = 1/median(diff(t));
    N = size(ang,2); fr = nan(N,1);
    for i=1:N
        a = ang(:,i);
        if all(~isfinite(a)), fr(i)=NaN; continue; end
        a = movmean(a, max(1,round(0.05*fs)));              % gentle smooth
        prom = 0.5*mad(a,1);                                % robust prominence
        [~,locs] = findpeaks(a,'MinPeakProminence',prom,'MinPeakDistance',round(0.15*fs));
        if numel(locs) < 2
            fr(i) = NaN;
        else
            Tsec = (t(locs(end))-t(locs(1))) / (numel(locs)-1);
            fr(i) = 1 / Tsec;
        end
    end
end

function [pF, pHolm, med3] = rm_stats(A, conds)
    % A: N x 3 (High, Mid, Low) values per animal
    N = size(A,1);
    [pF,~,~] = friedman(A,1,'off');
    % pairwise signrank with Holm correction
    pairs = [1 2; 1 3; 2 3];   % (H-M, H-L, M-L) in this order
    p = nan(3,1);
    for k=1:3
        [p(k),~] = signrank(A(:,pairs(k,1)), A(:,pairs(k,2)));
    end
    % Holm step-down
    [ps,ord] = sort(p); m = numel(p);
    holm = ps .* (m - (1:m)' + 1);
    adj  = holm;
    for i=2:m, adj(i)=max(adj(i-1), holm(i)); end
    adj = min(adj,1);
    pHolm = nan(size(p)); pHolm(ord) = adj;  % back to (H-M, H-L, M-L)
    med3 = median(A,1,'omitnan');
    fprintf('%s: Friedman p=%.3g | medians %s=%.3g, %s=%.3g, %s=%.3g | Holm p (H-M,H-L,M-L) = [%g %g %g]\n', ...
        inputname(1), pF, conds{1}, med3(1), conds{2}, med3(2), conds{3}, pHolm(1), pHolm(2), pHolm(3));
end

function paired_plot(A, conds, metric_name, pF, med3, pHolm, alpha)
    if nargin < 7 || isempty(pHolm), pHolm = [NaN NaN NaN]; end
    if nargin < 8 || isempty(alpha), alpha = 0.05; end

    figure('Color','w','Name',metric_name); hold on; grid on; box on;
    N = size(A,1); x = 1:3;
    cols = [0.20 0.70 0.20; 1.00 0.55 0.10; 0.00 0.45 1.00]; % High, Mid, Low

    % backdrops (IQR + median)
    for g=1:3
        q = quantile(A(:,g), [0.25 0.5 0.75]);
        w = 0.24;
        patch([x(g)-w x(g)+w x(g)+w x(g)-w], [q(1) q(1) q(3) q(3)], cols(g,:), ...
              'FaceAlpha',0.15,'EdgeColor',[.6 .6 .6]);
        plot([x(g)-w x(g)+w],[q(2) q(2)],'k-','LineWidth',2);
    end

    % paired lines + points
    for i=1:N
        plot(x, A(i,:), '-', 'Color', [.7 .7 .7], 'LineWidth', 1);
    end
    for g=1:3
        scatter(g*ones(N,1), A(:,g), 40, 'filled', 'MarkerFaceColor', cols(g,:), 'MarkerEdgeColor','k');
    end

    xlim([0.5 3.5]); set(gca,'XTick',x,'XTickLabel',conds);
    ylabel(strrep(metric_name,'_',' '));
    title(sprintf('Friedman p=%.3g | medians: %s=%.3g, %s=%.3g, %s=%.3g', ...
          pF, conds{1}, med3(1), conds{2}, med3(2), conds{3}));

    % ---------- significance bars (Holm-corrected) ----------
    % order in pHolm: [H-M, H-L, M-L] = pairs [1 2; 1 3; 2 3]
    pairs = [1 2; 1 3; 2 3];
    yData = A(:); yData = yData(isfinite(yData));
    if isempty(yData), yData = 0; end
    yMax  = max(yData);
    yMin  = min(yData);
    ySpan = max(yMax - yMin, eps);
    base  = yMax + 0.06*ySpan;     % where the first bar sits
    step  = 0.08*ySpan;            % vertical separation between bars
    hUsed = base;

    for k = 1:3
        p = pHolm(k);
        if ~isnan(p) && p < alpha
            x1 = pairs(k,1); x2 = pairs(k,2);
            y  = hUsed + (k-1)*step;
            plot([x1 x1 x2 x2], [y y+0.015*ySpan y+0.015*ySpan y], 'k-', 'LineWidth', 1.2);
            text((x1+x2)/2, y+0.02*ySpan, pstars(p), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', 'FontWeight','bold', 'Color','k');
            hTop = y + 0.08*ySpan;
        end
    end
    if exist('hTop','var')
        ylim([yMin - 0.05*ySpan, hTop + 0.10*ySpan]);
    end
end

function s = pstars(p)
    if p < 1e-3, s = '***';
    elseif p < 1e-2, s = '**';
    elseif p < 5e-2, s = '*';
    else, s = 'n.s.';  % not drawn unless you want to show it
    end
end
