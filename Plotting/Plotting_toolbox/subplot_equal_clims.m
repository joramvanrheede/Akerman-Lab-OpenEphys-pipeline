function subplot_equal_clims(colour_range)
% function SUBPLOT_EQUAL_CLIMS(COLOUR_RANGE)
% 
% SUBPLOT_EQUAL_CLIMS by itself sets all subplots in figure to use the same 
% color scale, equal to the maximum colour range between all plots. Does nothing 
% if there is only a single axes object in the figure.
% 
% If COLOUR_RANGE is specified as [cmin cmax], all subplots will be set to 
% this specified COLOUR RANGE.

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
    if nargin < 1
        % Get max and minimum values of the color scales for each axes object
        max_clim    = cellfun(@max,color_lims);
        min_clim    = cellfun(@min,color_lims);
        
        % Get the overal maximum and minimum between axes objects
        max_clim  	= max(max_clim);
        min_clim   	= min(min_clim);
    else
        min_clim    = colour_range(1);
        max_clim    = colour_range(2);
    end
    % Set all axes to overall maximum and minimum color range
    set(plotaxes,'CLim',[min_clim max_clim]);
else
    if nargin == 1
        set(plotaxes,'CLim',colour_range);
    end
end