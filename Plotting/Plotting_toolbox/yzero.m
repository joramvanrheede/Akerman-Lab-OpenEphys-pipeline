function yzero
% Sets the bottom of the y axis range to 0

y_lims = ylim;
ylim([0 y_lims(2)]);
