% Laser_Pulse analysis

data_dir                = '/Volumes/Akermanlab/Joram/Preprocessed data/AVK S1/Laser_pulse/2019_07_17_P2';
qreload                 = 1;
channels                = 1:10;             % Which channels to include in analysis 

% Histogram and spike count bins
hist_binsize            = [0.002];          % bin size in s for histograms and heatmap counts
direct_resp_win         = [0.007 0.025];    % response window in s - RELATIVE TO OPTO ONSET
effect_resp_win      	= [0.030 0.100];    % response window in s - RELATIVE TO OPTO ONSET
spont_win               = [-1 -.1];         % spontaneous window in s - RELATIVE TO OPTO ONSET
plot_win                = [-.2 .4];         % plotting window in s - RELATIVE TO OPTO ONSET
smoothwin               = 10;                % Window to smooth PSTH - affects psth plot only

power_bins              = [4 6 7 8 9 10];   % Used to group all experiments into easily digested subplots aggregating powers between these bin edges

% Remove spikes during suspected artifact window?
kill_artifact           = 1;                % remove artifact or not?
artifact_win            = [-0.002 0.007];   % window for removing artifact (set all spike times in this window to NaN)

% Colormaps:
color_max_percentile    = 0.5;              % Is used to shave a percentile off the data range before determining colour scale max

show_all_expts          = false;            % If 'true', shows a heatmap for every experiment in addition to those grouped using power_bins

% Saving figures
save_figs               = true;
fig_save_dir            = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/AVK/Pulse';


%% 
close all

if qreload
    
    folder_files    = dir([data_dir filesep '*.mat']); % read .mat files from folder
    
    % Pre-allocate
    laser_power     = [];
    laser_psth      = [];
    spike_resp      = [];
    spike_serr      = [];
    effect_resp     = [];
    effect_serr     = [];
    all_count_data  = [];
    max_count_data  = [];
    for a = 1:length(folder_files)
        
        this_file       = folder_files(a).name;
        full_filenm     = fullfile(data_dir, this_file);
        
        disp(['Loading ' this_file]);
        load(full_filenm);  % this should load variables 'ephys_data' and 'parameters'
        
        laser_power(a)  = ephys_data.parameters.LED_power;  % What laser/LED power was used?
        spikes          = ephys_data.conditions(1).spikes(channels,:,:);  % get spikes for target channels
        opto_onset      = ephys_data.conditions(1).LED_onset; % when did opto stimulus come on?
        
        spikes          = spikes - opto_onset; % make spike times relative to opto_onset
        q_artifact      = spikes >= artifact_win(1) & spikes <= artifact_win(2); % 
        
        if kill_artifact
            spikes(q_artifact) = NaN;
        end
        
        figure
        [im_h, spike_counts]    = psth(spikes,hist_binsize,plot_win); % get psth in chose time window
        laser_psth(a,:)         = spike_counts;
        
        
        [spike_resp(a), spike_serr(a)]  = spike_rate_in_win(spikes,direct_resp_win);
        spike_resp_spont                = spike_rate_in_win(spikes,spont_win);
        spike_resp(a)                   = spike_resp(a) - spike_resp_spont;
        
        
        [effect_resp(a), effect_serr(a)] = spike_rate_in_win(spikes,effect_resp_win);
        spike_resp_spont                = spike_rate_in_win(spikes,spont_win);
        spike_resp(a)                   = spike_resp(a) - spike_resp_spont;
        
        % Get spike counts over time by channel
        spike_bins              = [plot_win(1):hist_binsize:plot_win(2)];
        [imh, count_data]       = spike_density_plot(spikes,1,spike_bins);
        close(gcf)
        
        all_count_data          = cat(3,all_count_data,count_data);
        
    end
    
end

uniq_powers     = unique(laser_power);
n_powers        = length(uniq_powers);

%% Loop to make PSTH and get response for each laser power
mean_psth_by_power      = [];
spike_resp_by_power     = [];
spike_serr_by_power     = [];
effect_resp_by_power    = [];
effect_serr_by_power    = [];
for b = 1:n_powers
    this_power                  = uniq_powers(b);
    q_power                     = laser_power == this_power;
    
    these_psths                 = laser_psth(q_power,:);     % Select all PSTHs for this laser power
    
    mean_psth_by_power(b,:)     = smooth(mean(these_psths,1),smoothwin);    % Make mean PSTH trace and smooth for visualisation
    spike_resp_by_power(b)      = mean(spike_resp(q_power));
    spike_serr_by_power(b)      = mean(spike_serr(q_power));
    
    effect_resp_by_power(b)     = mean(effect_resp(q_power));
    effect_serr_by_power(b)     = mean(effect_serr(q_power));
end

% Create a vector of time bins to select relevant responses
psth_bins       = [plot_win(1):hist_binsize:plot_win(2)];

%% Plot PSTH
figure
plot(psth_bins(1:length(mean_psth_by_power)),mean_psth_by_power','LineWidth',2)
title('PSTH of spikes following LASER pulse')
xlabel('Post-pulse Time (s)')
ylabel('Spike count in time bin')
fixplot
xlim(plot_win)
set(gcf,'Units','Normalized','Position',[0.2 0.3 0.6 0.4],'PaperPositionMode','auto')

power_labels = [];
for a = 1:length(uniq_powers)
    power_labels{a} = num2str(uniq_powers(a));
end
legend(power_labels)

%% Plot laser power
figure
errorbar(uniq_powers,spike_resp_by_power,spike_serr_by_power,'k.-','MarkerSize',25,'LineWidth',2)
title('LASER power vs spike rate, EARLY')
xlabel('LASER Power')
ylabel('Delta spike rate (Hz) vs. spontaneous')
fixplot
xlim([min(uniq_powers) - 1, max(uniq_powers) + 1])

%% Plot laser power
figure
errorbar(uniq_powers,effect_resp_by_power,effect_serr_by_power,'k.-','MarkerSize',25,'LineWidth',2)
title('LASER power vs spike rate, LATE')
xlabel('LASER Power')
ylabel('Delta spike rate (Hz) vs. spontaneous')
fixplot
xlim([min(uniq_powers) - 1, max(uniq_powers) + 1])

%% Make figure of heat maps by power range

[sort_laser_power, sortinds] = sort(laser_power);   % Sort all the opto powers
all_count_data = all_count_data(:,:,sortinds);      % Sort the count data based on opto power

figure
n_power_figs            = length(power_bins) - 1;
ax_h                    = tight_subplot(n_power_figs,1,0.05,0.05,0.05);
for i = 1:n_power_figs
    q_power_win         = sort_laser_power > power_bins(i) & sort_laser_power <= power_bins(i+1); 
    mean_count_data     = mean(all_count_data(:,:,q_power_win),3);
    
    axes(ax_h(i))
    imagesc(mean_count_data)
    colormap hot
    
    % set a sensible maximum value for the color scale; shave off top percentile which can shift the whole colormap
    max_count_data      = robust_max(mean_count_data,color_max_percentile,'all');
    
    % if there are literally no spikes outside of this percentile, set the range to be between 0 and 1
    if max_count_data < 1 || isnan(max_count_data)
        max_count_data = 1;
    end
    
    set(gca,'CLim',[0 max_count_data])
    
    title(['Opto power '  num2str(power_bins(i)) ' - ' num2str(power_bins(i+1))])
    
    % X axis tick labels will currently label the pixel nr in the image
    x_tick_labels   = xticklabels;
    
    % Substitute the current X axis tick labels with the relevant bin so that
    % they reflect the input data (spike times) rather than the image pixel nr.
    new_x_tick_labels   = cell(size(x_tick_labels));
    for j = 1:length(x_tick_labels)
        x_tick_ind              = str2num(x_tick_labels{j});
        x_tick_val              = spike_bins(x_tick_ind+1);
        new_x_tick_labels{j}    = num2str(x_tick_val);
    end
    
    % Replace the tick labels and do some other aesthetic stuff
    set(gca,'TickDir','out','Box','off','XTickLabel',new_x_tick_labels)
    
    colorbar
    
end
% Set all axes in figure to same colour range:
subplot_equal_clims

set(gcf,'Units','Normalized','Position',[0.3 0.1 0.4 0.8],'PaperPositionMode','auto')


%% Make figure with heat map for each experiment
if show_all_expts
    % Make figure with subplots that extend all the way to the edge
    figure
    ax_h    = tight_subplot(length(sort_laser_power),1,0,0,[.1 0]);
    
    % Loop over the laser powers in ascending order
    for a = 1:length(sort_laser_power)
        
        % Select relevant spike count data
        this_count_data 	= all_count_data(:,:,a);
        
        % set a sensible maximum value for the color scale; shave off top percentile which can shift the whole colormap
        max_count_data      = robust_max(this_count_data,color_max_percentile,'all');
        
        % if there are literally no spikes outside of this percentile, set the range to be between 0 and 1
        if max_count_data < 1
            max_count_data = 1;
        end
        
        % Make 'spike density' heatmap of each of the individual experiments
        axes(ax_h(a))
        imagesc(this_count_data)
        xlabel(num2str(sort_laser_power(a)))
        colormap hot
        set(gca,'CLim',[0 max_count_data])
    end
    
    % Set subplots to the same colour range
    subplot_equal_clims
end


%%

min_chan                    = num2str(min(channels));
max_chan                    = num2str(max(channels));
[filepath, expt_name, ext]  = fileparts(data_dir);

if save_figs
    print(1,fullfile(fig_save_dir,['Pulse ' expt_name ' ' min_chan '-' max_chan  ' PSTHs']),'-dpng')
    print(2,fullfile(fig_save_dir,['Pulse ' expt_name ' ' min_chan '-' max_chan  ' DIRECT power plot']),'-dpng')
    print(3,fullfile(fig_save_dir,['Pulse ' expt_name ' ' min_chan '-' max_chan  ' EFFECT power plot']),'-dpng')
    print(4,fullfile(fig_save_dir,['Pulse ' expt_name ' ' min_chan '-' max_chan  ' density plots']),'-dpng')

end


