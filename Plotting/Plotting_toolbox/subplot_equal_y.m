function subplot_equal_y(ymax)
% function subplot_equal_y
% Sets all subplots in figure to use the same y axis; does nothing if there
% is only a single plot in the figure so far.

plotaxes    = get(gcf,'Children');

y_limits    = get(plotaxes,'Ylim');

if iscell(y_limits)
    maxy        = cellfun(@max,y_limits);
    miny        = cellfun(@min,y_limits);
    
    maxy        = max(maxy);
    miny        = min(miny);
    
    if nargin < 1
        set(plotaxes,'YLim',[miny maxy]); % add margin
    elseif nargin == 1
        set(plotaxes,'YLim',[miny ymax]); % add margin
    end
end