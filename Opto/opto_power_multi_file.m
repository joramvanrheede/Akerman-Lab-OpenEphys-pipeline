% Laser_Pulse analysis

data_dir        = '/Volumes/Akermanlab/Joram/Preprocessed data/AVK POM/Laser_pulse/2019_07_19_P1';

qreload         = 1;

% Do LFP?

hist_binsize    = [0.005];
rate_window     = [0.020];

resp_win        = [0.007 0.150]; %
spont_win       = [0 .9];
plot_win        = [-.5 1.5];
win_range       = [-1 1.5];
smoothwin       = 1;




color_scaling   = 0.05;

channels        = 1:8;
trials          = 'all';

kill_artifact   = 1;
artifact_win    = [-0.001 0.006];

color_max_percentile    = 0.1;
subplot_margin          = 0;

%% 
close all


if qreload
    
    folder_files    = dir([data_dir filesep '*.mat']); % read .mat files from folder
    
    % Pre-allocate
    laser_power     = [];
    laser_psth      = [];
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
        [im_h, spike_counts]    = psth(spikes,hist_binsize,win_range); % get psth in chose time window
        laser_psth(a,:)         = spike_counts;
        
        % Get spike counts over time by channel
        [imh, count_data]       = spike_density_plot(spikes,1,[plot_win(1):hist_binsize:plot_win(2)]);
        close(gcf)
        
        all_count_data          = cat(3,all_count_data,count_data);
        
    end
    
end

[sort_laser_power, sortinds] = sort(laser_power);   % Sort all the opto powers
all_count_data = all_count_data(:,:,sortinds);      % Sort the count data based on opto power

% Make figure with subplots that extend all the way to the edge
figure
ax_h    = tight_subplot(length(sort_laser_power),1,subplot_margin,subplot_margin,[.1 subplot_margin]);

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

% Create a vector of time bins to select relevant responses
psth_bins       = [win_range(1):hist_binsize:win_range(2)];
q_respwin       = psth_bins >= resp_win(1) & psth_bins <= resp_win(2);
q_spontwin   	= psth_bins >= spont_win(1) & psth_bins <= spont_win(2);

uniq_powers     = unique(laser_power);
n_powers        = length(uniq_powers);


%% Loop to make PSTH for each laser power
mean_psth_by_power  = [];
spike_sum_by_power  = [];
spont_sum_by_power  = [];
for b = 1:n_powers
    this_power                  = uniq_powers(b);
    q_power                     = laser_power == this_power;
    
    these_psths                 = laser_psth(q_power,:);     % Select all PSTHs for this laser power
    
    mean_psth_by_power(b,:)     = smooth(mean(these_psths,1),smoothwin);    % Make mean PSTH trace and smooth for visualisation
    spike_sum_by_power(b)       = mean(mean(these_psths(:,q_respwin),1));   % Get spike rate in target response window
    spont_sum_by_power(b)       = mean(mean(these_psths(:,q_spontwin),1));  % Get spike rate in spontaneous window
end

% plot all psths
figure
plot(psth_bins(1:length(mean_psth_by_power)),mean_psth_by_power')
title('PSTH of spikes following LASER pulse')
xlabel('Post-pulse Time (s)')
ylabel('Spike count in time bin')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')
xlim(plot_win)


figure
plot(uniq_powers,spike_sum_by_power - spont_sum_by_power,'k.-','MarkerSize',20)
title('LASER Power setting vs binned spiking response')
xlabel('LASER Power setting')
ylabel('Spiking response')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')
