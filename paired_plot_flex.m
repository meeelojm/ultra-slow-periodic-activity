function paired_plot_flex(A, conds, metric_name, testinfo)
% A: N×G (G=2 or 3). Draws quantile boxes + paired lines + points.
% testinfo: struct with fields .test (char), .p (double), .stat (struct or [])
    arguments
        A double
        conds cell
        metric_name char
        testinfo struct = struct('test','','p',NaN,'stat',[])
    end
    [N,G] = size(A);
    assert(G==2 || G==3, 'This helper supports 2 or 3 conditions.');
    if numel(conds) ~= G, error('conds must have %d labels', G); end

    % Colors (High, Mid, Low palette trimmed to G)
    baseCols = [0.20 0.70 0.20; 1.00 0.55 0.10; 0.00 0.45 1.00];
    cols = baseCols(1:G,:);

    figure('Color','w','Name',metric_name); hold on; grid on; box on;
    x = 1:G; w = 0.24;

    % Backdrop quantile boxes per condition
    for g = 1:G
        q = quantile(A(:,g), [0.25 0.5 0.75], 1);
        patch([x(g)-w x(g)+w x(g)+w x(g)-w], [q(1) q(1) q(3) q(3)], cols(g,:), ...
              'FaceAlpha',0.15,'EdgeColor',[.6 .6 .6]);
        plot([x(g)-w x(g)+w],[q(2) q(2)],'k-','LineWidth',2);
    end

    % Paired lines (one per fish)
    for i=1:N
        plot(x, A(i,:), '-', 'Color', [.7 .7 .7], 'LineWidth', 1);
    end

    % Points
    for g=1:G
        scatter(g*ones(N,1), A(:,g), 40, 'filled', ...
                'MarkerFaceColor', cols(g,:), 'MarkerEdgeColor','k');
    end

    xlim([0.5 G+0.5]); set(gca,'XTick',x,'XTickLabel',conds);
    ylabel(strrep(metric_name,'_',' '));

    % Title with stats
    meds = median(A,1,'omitnan');
    if strcmpi(testinfo.test,'signrank') && isfield(testinfo,'p')
        if isfield(testinfo,'stat') && isstruct(testinfo.stat) && isfield(testinfo.stat,'signedrank')
            statStr = sprintf('signedrank=%d', testinfo.stat.signedrank);
        else
            statStr = '';
        end
        title(sprintf('Wilcoxon signed-rank p=%.3g %s | medians: %s=%.3g, %s=%.3g', ...
              testinfo.p, statStr, conds{1}, meds(1), conds{2}, meds(2)));
    else
        if G==3
            title(sprintf('Medians: %s=%.3g, %s=%.3g, %s=%.3g', ...
                  conds{1}, meds(1), conds{2}, meds(2), conds{3}, meds(3)));
        else
            title(sprintf('Medians: %s=%.3g, %s=%.3g', ...
                  conds{1}, meds(1), conds{2}, meds(2)));
        end
    end
end
