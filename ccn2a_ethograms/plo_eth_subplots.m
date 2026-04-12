function plo_eth_subplots(sta_tim, end_tim, ...
    dru_ons, sti_ons, tai_bea, hea_bea, ope_bea, mou_bea, eye_bea, ...
    tai_sig, hea_sig, ope_sig, mou_sig, eye_sig, ...
    rat_tim, rat_tim_par, fra_tim, tim_rob, rat_rob, spo_bou_ons, ...
    col_tai, col_hea, col_ope, col_mou, col_eye, ...
    col_tai_sig, col_hea_sig, col_ope_sig, col_mou_sig, col_eye_sig, ...
    fon_siz, col_rob, col_dru, col_sti, col_spo_bou_ons, tit_rec)

% Figure & layout
f = figure('Color','w','Name',tit_rec);
t = tiledlayout(f, 5, 1, 'TileSpacing','compact','Padding','compact');
title(t, tit_rec, 'FontWeight','bold');

% unpack rates (columns: tail, heart, operc, mouth, eye)
tai_rat = rat_tim_par(:,1);
hea_rat = rat_tim_par(:,2);
ope_rat = rat_tim_par(:,3);
mou_rat = rat_tim_par(:,4);
eye_rat = rat_tim_par(:,5);

% Row 1: Tail
nexttile;
plot_row(gca, 'Tail', fra_tim, tai_sig, col_tai_sig, rat_tim, tai_rat, col_tai, tai_bea);
draw_vlines(gca, dru_ons, col_dru, '--', .5);
draw_vlines(gca, sti_ons, col_sti, ':',  .5);
draw_vlines(gca, spo_bou_ons, col_spo_bou_ons, '--', .5);
xlim([sta_tim end_tim]);

% Row 2: Heart (also robust HR)
ax = nexttile;
plot_row(ax, 'Heart', fra_tim, hea_sig, col_hea_sig, rat_tim, hea_rat, col_hea, hea_bea);
% robust HR on right axis as extra line (scaled to that axis)
if ~isempty(tim_rob) && ~isempty(rat_rob)
    yyaxis(ax,'right'); hold on;
    plot(tim_rob, rat_rob, 'Color', col_rob, 'LineWidth', 1);
end
draw_vlines(ax, dru_ons, col_dru, '--', .5);
draw_vlines(ax, sti_ons, col_sti, ':',  .5);
draw_vlines(ax, spo_bou_ons, col_spo_bou_ons, '--', .5);
xlim([sta_tim end_tim]);

% Row 3: Operculum
nexttile;
plot_row(gca, 'Operculum', fra_tim, ope_sig, col_ope_sig, rat_tim, ope_rat, col_ope, ope_bea);
draw_vlines(gca, dru_ons, col_dru, '--', .5);
draw_vlines(gca, sti_ons, col_sti, ':',  .5);
draw_vlines(gca, spo_bou_ons, col_spo_bou_ons, '--', .5);
xlim([sta_tim end_tim]);

% Row 4: Mouth
nexttile;
plot_row(gca, 'Mouth', fra_tim, mou_sig, col_mou_sig, rat_tim, mou_rat, col_mou, mou_bea);
draw_vlines(gca, dru_ons, col_dru, '--', .5);
draw_vlines(gca, sti_ons, col_sti, ':',  .5);
draw_vlines(gca, spo_bou_ons, col_spo_bou_ons, '--', .5);
xlim([sta_tim end_tim]);

% Row 5: Eye
ax5 = nexttile;
plot_row(ax5, 'Eye', fra_tim, eye_sig, col_eye_sig, rat_tim, eye_rat, col_eye, eye_bea);
draw_vlines(ax5, dru_ons, col_dru, '--', .5);
draw_vlines(ax5, sti_ons, col_sti, ':',  .5);
draw_vlines(ax5, spo_bou_ons, col_spo_bou_ons, '--', .5);
xlim([sta_tim end_tim]);

% X-label only on bottom
xlabel(ax5, 'Time (s)');   % if you prefer minutes: xlabel(ax5,'Time (min)'); and divide all x by 60 above

% prettier axes
ax = findall(f,'Type','axes');
set(ax, 'Layer','top','FontSize',fon_siz,'Box','off','TickDir','out','TickLength',[.008 .008]);

linkaxes(ax,'x'); % pan/zoom sync
end

% ---------- helpers ----------
function plot_row(ax, name, fra_tim, sig, col_sig, rat_tim, rat, col_rate, events)
% Left: signal (if provided)
axes(ax); cla(ax); hold(ax,'on');
yyaxis(ax,'left');
if ~isempty(sig) && ~isempty(fra_tim)
    plot(ax, fra_tim(:), normalize_01(sig(:)), 'Color', col_sig, 'LineWidth', 1);
end
ylabel(ax, name);

% Right: rate
yyaxis(ax,'right');
if ~isempty(rat) && ~isempty(rat_tim)
    plot(ax, rat_tim(:), rat(:), 'Color', col_rate, 'LineWidth', 1);
end
% event ticks along baseline (right axis)
if ~isempty(events)
    y0 = min_ylim(ax);
    yyaxis(ax,'right');  % ensure correct axis
    scatter(ax, events(:), repmat(y0, numel(events),1), 12, col_rate, 'filled', 'MarkerEdgeColor','none');
end
grid(ax,'off');
end

function draw_vlines(ax, xs, col, sty, lw)
if nargin<4||isempty(sty), sty='--'; end
if nargin<5||isempty(lw),  lw=.5;     end
if isempty(xs), return; end
yl = ylim(ax); hold(ax,'on');
for k=1:numel(xs)
    line(ax, [xs(k) xs(k)], yl, 'Color', col, 'LineStyle', sty, 'LineWidth', lw);
end
end

function y = normalize_01(x)
x = x(:);
mn = min(x); mx = max(x);
if ~isfinite(mn) || ~isfinite(mx) || mx<=mn
    y = zeros(size(x));
else
    y = (x - mn) / (mx - mn);
end
end

function y0 = min_ylim(ax)
yl = ylim(ax);
y0 = yl(1) + 0.03*(yl(2)-yl(1));
end
