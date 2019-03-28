function shaded_region(xvals,colour,alpha_val)
% function shaded_region(xvals,colour,alpha_val)
% 
% Makes shaded region in a plot between specified x values.
% Uses area plot function to generate plot and will shade a region spanning
% the entire y axis (but not more).
%
% Plots in current axes.
% 
% XVALS is a 2-element vector [xmin xmax].
% 
% Other inputs are optional - output defaults to 'cyan' and alpha of 0.7.
% 
% COLOUR sets the colour, you can use matlab colour string specifications 
% ('r' 'g' 'b' 'c', etc) or  [r g b] triplets.
% 
% ALPHA_VAL sets transparency (1 = opaque, 0 = fully transparent).
%
% Joram van Rheede 21/03/2019

% Set sensible defaults
if nargin < 2
    colour      = 'c'; % cyan, because it's an underrated colour
end

if nargin < 3
    alpha_val   = 0.7; % for clearly visible shaded region
end

% Check if 'hold on' is true for current graph
if ~ishold
    hold on
    set_hold_off    = true;
else
    set_hold_of     = false;
end
    
% Get the values for the top and bottom of y-axis which will determine the
% extent of the shaded region along the y-axis
yvals   = ylim;
ymin    = yvals(1);
ymax    = yvals(2);

xmin    = xvals(1);
xmax    = xvals(2);

% Use area plot to create shaded region of defined colour and alpha level
area([xmin xmin xmax xmax],[ymin ymax ymax ymin],ymin,'FaceColor',colour,'FaceAlpha',alpha_val,'EdgeAlpha',0)

% If hold was 'off' originally, reset hold to original hold status
if set_hold_off
    hold off
end
