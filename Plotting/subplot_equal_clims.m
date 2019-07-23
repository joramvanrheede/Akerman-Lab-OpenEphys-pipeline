function subplot_equal_clims
% function subplot_equal_clims
% Sets all subplots in figure to use the same color scale; does nothing if 
% there is only a single axes object in the figure.

% Get the 'children' of the figure object (any object in the figure)
plotchildren    = get(gcf,'Children');

% Make sure to select only those objects that are axes (the subplots)
obj_type        = get(plotchildren,'Type');
q_axes          = strcmp(obj_type,'axes');
plotaxes        = plotchildren(q_axes);

% get the CLim property
color_lims      = get(plotaxes,'Clim');

% Conditional - if there is only one axes object in the figure, color_lims
% is a 2-element vector
if iscell(color_lims)
    
    % Get max and minimum values of the color scales for each axes object
    max_clim    = cellfun(@max,color_lims);
    min_clim    = cellfun(@min,color_lims);
    
    % Get the overal maximum and minimum between axes objects
    max_clim  	= max(max_clim);
    min_clim   	= min(min_clim);
    
    % Set all axes to overall maximum and minimum color range
    set(plotaxes,'CLim',[min_clim max_clim]);
end