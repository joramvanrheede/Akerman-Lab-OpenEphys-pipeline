function image_handle = spike_density_plot(spike_data,split_dim,bins)
% function line_handles = spike_density_plot(spike_data,split_dim)
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

%% use histcounts to make a histogram for each row of the plot, stack them together, then use imagesc to plot...

if nargin < 2
    split_dim = 1;
end

if nargin < 3
    spacing = 0.2;
end

n_rows  = size(spike_data, split_dim);

count_data = zeros(n_rows, length(bins)-1);
for a = 1:n_rows
    if split_dim == 1
        trial_spikes      	= spike_data(a,:,:);
    elseif split_dim == 2
        trial_spikes      	= spike_data(:,a,:);
    end
    
    count_data(a,:)     = histcounts(trial_spikes,bins,'Normalization','probability');
    
end

image_handle = imagesc(count_data);
colormap hot
axis tight
axis ij


