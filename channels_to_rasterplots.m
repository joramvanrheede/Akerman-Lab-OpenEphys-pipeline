function fighandles = channels_to_rasterplots(channelfile,chans,split_plots,split_figures,trialrange,x_ax_lims,condition_name,condition_units)
% channels_to_rasterplots(channelfile,chans,split_plots,split_figures,trialrange,x_ax_lims)
% Figures with rasterplots separated by conditions
% channelfile:      Full file name of 'channels' data struct
% chans:            Which channels to include (e.g. [1:4])
% split_plots:      Which conditions to split plots by?
% split_figures:    Which conditions to split figures by?
% trialrange:       Which trials to show (e.g. [1 10])
% x_ax_lims:        x axes limits in seconds
% 
% Example values for Frequency experiment:
%
% channelfile    	= '/Users/Joram/Dropbox/Akerman Postdoc/Data/Extracted data/Dual Stim/Frequency/2018_04_24/2018_04_24-4-Frequency.mat';
%
% chans           	= [1:16];
% split_plots     	= [4 1];
% split_figures   	= [6];
% trialrange      	= [1 12];
% x_ax_lims       	= [0 12];

trialnrs            = trialrange(1):trialrange(end);

load(channelfile)

%% Get spikes into easily addressable, 'plotSpikeRaster'-compatible format
spike_cell = {};
for a = 1:length(chans)
    chan = chans(a);
    for b = 1:length(channels(chan).conditions)
        for c = 1:length(channels(chan).conditions(b).episodes)
            spike_cell(a,b,c) = {channels(chan).conditions(b).episodes(c).spikes'};
        end
    end
end

qempty              = cellfun(@isempty,spike_cell);
spike_cell(qempty)  = deal({[0]});

%% Retrieve conditions 
condition_mat           = cell2mat({channels(1).conditions.timings}');

% which condition to split figures by?
split_figure_mat        = condition_mat(:,split_figures);

% work out which conditions go in which figure
[split_fig_rows, indxa, cond_fig_inds] = unique(split_figure_mat,'rows');

fighandles = [];
for a = 1:length(split_fig_rows)
    
    % Create new figure and set it up 
    fighandles(a) = figure;
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.2 .1 .6 .9])
    set(gcf,'Color',[1 1 1])
    
    % select appropriate subset of conditions
    qselect         = condition_mat(:,split_figures) == split_fig_rows(a);
    thiscondmat     = condition_mat(qselect,:);
    
    % split this subset of conditions by split_plots
    split_plot_mat  = thiscondmat(:,split_plots);
    [split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');
    
    % number of rows and columns in subplots
    n_cols      = length(unique(split_plot_rows(:,2)));
    n_rows      = length(unique(split_plot_rows(:,1)));
    
    % use 'tight_subplot' function to ensure more effective plotting area
    % in multiple subplots
    [axh axpos] = tight_subplot(n_rows,n_cols,0.02,0.02,0.02);
    
    for b = 1:length(split_plot_rows)
        
        % get the current row of condition values
        this_row    = split_plot_rows(b,:);
        
        % make boolean to select appropriate condition from spike_cell
        qcond       = condition_mat(:,split_plots(1)) == this_row(1) & condition_mat(:,split_plots(2)) == this_row(2) & qselect;
        
        % select appropriate axes
        axes(axh(b)) 
        
        % fetch appropriate trials from spike_cell
        spike_episodes  = spike_cell(:,qcond,trialnrs);
        
        % remove spurious dimensions
        spike_episodes  = squeeze(spike_episodes);
        
        % reshape spike_episodes such that spike responses are grouped by
        % channel
        spike_episodes  = spike_episodes';
        spike_episodes  = spike_episodes(:);
        
        % create plot
        plotSpikeRaster(spike_episodes,'PlotType','vertline','XLimForCell',x_ax_lims,'VertSpikeHeight',.8);
        
        
        set(gca,'xtick',x_ax_lims(1):x_ax_lims(end),'ytick',[0:max(trialnrs):max(trialnrs)*length(chans)],'BoxStyle','full')
        title([condition_name ' = ' num2str(condition_mat(qcond,split_plots(1))) condition_units])
    end
end


