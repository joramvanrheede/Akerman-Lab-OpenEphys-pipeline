function [image_handle, count_data] = spike_density_plot(spike_data,split_dim,bins)
% FUNCTION [IMAGE_HANDLE, COUNT_DATA] = SPIKE_DENSITY_PLOT(SPIKE_DATA,SPLIT_DIM,BINS)
% 
% Make heatmap of spike density over time by trial or by channel; 
% effectively this turns a set of PSTHs into a visually pleasing image that
% may be easier to interpret than a rasterplot or superimposed PSTHs.
% 
% SPIKE_DATA: M*N*O matrix of spike times. In current Akerman lab
% 'ephys_data' format, the dimensions represent:
% M = number of channels
% N = number of trials
% O = number of spikes
% 
% NaNs are ignored (it is assumed that NaNs represent empty slots in the 
% SPIKE_DATA matrix.
% 
%
% SPLIT_DIM: Determines which dimension of SPIKE_DATA to display as individual 
% rows, and which to average over. E.g. for SPLIT_DIM == 1, output is a
% row for each channel and averages over N trials. SPLIT_DIM == 2 outputs a
% row for each individual trial and averages over M channels.
% (the only options are 1 or 2)
% 
% 
% BINS: specifies the time bins in which SPIKE_DATA points are pooled to 
% generate the histogram that determines the colour in the colour map.
% e.g. for BINS == [0:0.001:1] the function outputs a spike density plot
% of spikes falling between 0 and 1s, with 1ms time bins.
% 
% 
% IMAGE_HANDLE: A handle to the image plotted with MATLAB's IMAGESC function.
% 
% COUNT_DATA: The histogram count data that underlies the heatmap in the
% image. Use imagesc(COUNT_DATA) to recreate the heatmap.
% 
% Joram van Rheede - March 2019

% Default: split by channel, average over trials
if nargin < 2 || isempty(split_dim)
    split_dim = 1;
end

if nargin < 3
    % generate binedges using all spike data, divide into 100 bins which is
    % a good amount for the matlab figure size
    [temp_counts, bins] = histcounts(spike_data(:),100);
end

% determine size of the heatmap
n_rows  = size(spike_data, split_dim);
n_bins  = length(bins);

% initialise count_data that will make up the heatmap image
count_data = zeros(n_rows, n_bins-1); 

% loop over each row
for a = 1:n_rows
    
    % Get the relevant spike data
    if split_dim == 1
        trial_spikes      	= spike_data(a,:,:);
    elseif split_dim == 2
        trial_spikes      	= spike_data(:,a,:);
    end
    
    % Make histogram for this row
    count_data(a,:)     = histcounts(trial_spikes,bins);
    
end

% Display the count data as an image
image_handle    = imagesc(count_data);

% Do plot aesthetics BEFORE changing tick labels (font size change will change the number of ticks and mess up labeling)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')

% X axis tick labels will currently label the pixel nr in the image
x_tick_labels   = xticklabels;

% Substitute the current X axis tick labels with the relevant bin so that
% they reflect the input data (spike times) rather than the image pixel nr.
new_x_tick_labels   = cell(size(x_tick_labels));
for i = 1:length(x_tick_labels)
    x_tick_ind              = str2num(x_tick_labels{i});
    x_tick_val              = bins(x_tick_ind+1);
    new_x_tick_labels{i}    = num2str(x_tick_val);
end

% Replace the tick labels and do some other aesthetic stuff
set(gca,'TickDir','out','Box','off','XTickLabel',new_x_tick_labels)

% Because your spike data are hot stuff:
colormap hot

% invert y axis so that figure can be read from top to bottom (i.e. top row
% is the first / most superficial channel, or the first trial / sweep)
axis ij 


