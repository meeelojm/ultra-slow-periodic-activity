function draw_scalebar_xy(ax, len_um, um_per_unit)
    % L-shaped scalebar: horizontal → right, vertical → up (on screen)
    if nargin<3 || isempty(um_per_unit), um_per_unit = 1; end
    L = len_um / um_per_unit;                 % length in data units

    xlm  = get(ax,'XLim');
    ylm  = get(ax,'YLim');
    padX = 0.04 * diff(xlm);
    padY = 0.06 * diff(ylm);

    % Corner at bottom-left (data coords). With YDir 'reverse', smaller y is up.
    x0 = xlm(1) + padX;        % corner x
    y0 = ylm(2) - padY;        % baseline y

    hold(ax,'on');
    % horizontal leg: to the right
    plot(ax, [x0, x0+L], [y0, y0], 'k-', 'LineWidth', 2);
    % vertical leg: upward on screen (decrease y)
    plot(ax, [x0, x0], [y0, y0 - L], 'k-', 'LineWidth', 2);

    % offsets (tweak if needed)
    dx = 0.012*diff(xlm);
    dy = 0.010*diff(ylm);
    
    % horizontal label: just BELOW the horizontal leg, starting at the corner
    text(ax, x0, y0 + dy, sprintf('%d microns', len_um), ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'Color','k','FontSize',9);
    
    % vertical label: LEFT of the vertical leg, aligned to its bottom end
    text(ax, x0 - dx, y0, sprintf('%d microns', len_um), ...
        'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'Color','k','FontSize',9);


end
