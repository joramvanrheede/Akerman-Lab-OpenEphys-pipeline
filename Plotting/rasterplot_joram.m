function line_handles = rasterplot_joram(spike_data,split_dim,spacing)
% function line_handles = rasterplot_joram(spike_data,split_dim,spacing)
% m*n*o matrix of spike time data
% takes one input dimension to split rasterplot into rows
% m = number of channels
% n = number of trials
% o = number of spikes
% you can specify 'split_dim' to choose which dimension to use to split
% rows of the rasterplot
% it is assumed that empty spike time values are NaNs
% 
% Joram van Rheede - Jan 2019

if nargin < 2
    split_dim = 1;
end

if nargin < 3
    spacing = 0.2;
end

n_rows  = size(spike_data, split_dim);

for a = 1:n_rows
    if split_dim == 1
        trial_spikes      	= spike_data(a,:,:);
    elseif split_dim == 2
        trial_spikes      	= spike_data(:,a,:);
    end
    xvals               = trial_spikes(:);
    %xvals               = xvals(~isnan(xvals));
    
    x_raster            = NaN(length(xvals)*3,1);
    x_raster(1:3:end)   = xvals;
    x_raster(2:3:end)   = xvals;
    
    y_raster            = NaN(length(xvals)*3,1);
    y_raster(1:3:end)   = a-spacing/2;
    y_raster(2:3:end)   = a-(1-spacing/2);
    
    line_handles(a)     = line(x_raster(:),y_raster(:),'Color',[0 0 0]);
    
end

axis tight
axis ij


