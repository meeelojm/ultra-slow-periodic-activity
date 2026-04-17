function plot_itpc_psth_like(itpcCurves, tSec, sortMode, recordingLabel)

    if nargin < 3 || isempty(sortMode)
        sortMode = 'reset';
    end
    if nargin < 4
        recordingLabel = '';
    end

    % --- define windows for sorting ---
    preMask  = tSec >= -40 & tSec <= -5;
    postMask = tSec >= 5   & tSec <= 40;

    % --- sorting variable ---
    switch lower(sortMode)
        case 'reset'
            sortVal = mean(itpcCurves(:,postMask),2,'omitnan') - ...
                      mean(itpcCurves(:,preMask),2,'omitnan');
        case 'post'
            sortVal = mean(itpcCurves(:,postMask),2,'omitnan');
        case 'peak'
            sortVal = max(itpcCurves(:,postMask),[],2);
        otherwise
            sortVal = mean(itpcCurves(:,postMask),2,'omitnan') - ...
                      mean(itpcCurves(:,preMask),2,'omitnan');
    end

    [~, ord] = sort(sortVal, 'descend');
    itpcSorted = itpcCurves(ord,:);

    % --- mean and SEM ---
    m = mean(itpcCurves, 1, 'omitnan');
    s = std(itpcCurves, 0, 1, 'omitnan') ./ sqrt(size(itpcCurves,1));

    % --- figure ---
    figure('Color','w', 'Position', [100 100 500 700]);
    tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

    % =========================
    % Top: mean ITPC trace
    % =========================
    nexttile; hold on;

    fill([tSec(:); flipud(tSec(:))], ...
         [m(:)-s(:); flipud(m(:)+s(:))], ...
         [0.8 0.8 0.8], ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.6);

    plot(tSec, m, 'k', 'LineWidth', 2);
    xline(0, '--k', 'LineWidth', 1);

    ylabel('ITPC');
    title('Bout-triggered phase concentration');
    xlim([min(tSec) max(tSec)]);
    box off;

    % =========================
    % Bottom: neuron heatmap
    % =========================
    nexttile;

    imagesc(tSec, 1:size(itpcSorted,1), itpcSorted);
    set(gca, 'YDir', 'normal');
    xline(0, '--k', 'LineWidth', 1);

    xlabel('Time from bout onset (s)');
    ylabel('Neuron #');
    title(sprintf('Neurons sorted by %s ITPC', sortMode));

    colormap(gca, hot);
    cb = colorbar;
    ylabel(cb, 'ITPC');

    xlim([min(tSec) max(tSec)]);

    if ~isempty(recordingLabel)
        sgtitle(recordingLabel, 'Interpreter', 'none');
    end
end