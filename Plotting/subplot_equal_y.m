function subplot_equal_y
% sets all subplots in figure to use the same y axis

plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
miny        = cellfun(@min,get(plotaxes,'Ylim'));

maxy        = max(maxy);
miny        = min(miny);

set(plotaxes,'YLim',[miny maxy]); % add margin
