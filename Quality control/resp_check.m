
% check_whisk_resp_baseline

% % 03/04
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-ON/2019_04_03/2019_04_03-17-Laser_power.mat' - YES, Post = REDUCED; 24%

% % 04/04
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/CHRONOS-ON/2019_04_04/2019_04_04-1-Laser_power.mat' - no - weird file, wrong led / whisk times
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-ON/2019_04_04/2019_04_04-2-Laser_power.mat' - YES, 813%
% 
% % 07/04
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-ON/2019_04_07/2019_04_07-1-Laser_power.mat' - YES



% 03/05
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-OFF/2019_05_03/2019_05_03-3-Laser_power.mat' - YES
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-OFF/2019_05_03/2019_05_03-10-Laser_power.mat' - YES

% 06/05 - 1
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-OFF/2019_05_06/2019_05_06-3-Laser_power.mat' - YES, most convincing profile...
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-OFF/2019_05_06/2019_05_06-9-Laser_power.mat' - nearly

% 07/05 - 1
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-OFF/2019_05_07_1/2019_05_07-5-Laser_power.mat' - YES
% '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-OFF/2019_05_07_1/2019_05_07-11-Laser_power.mat' - nearly


data_file           = '/Volumes/Akermanlab/Joram/Preprocessed data/MJB laser power/Laser_power/CHRONOS-ON/2019_04_03/2019_04_03-17-Laser_power.mat';
q_reload            = 1;

%% Visualisation settings
psth_win            = [-.05 .1];
bin_size            = [0.001];

stim_type           = 'whisk'; % 'whisk' or 'opto'

%% Responsiveness settings
spike_resp_win    	= [0.005 0.08]; % where to measure response (relative to stim onset)
spike_spont_win  	= [-1 -.1]; % where to measure spontaneous activity (relative to stim onset)

n_trials_for_drift  = 10; % number of trials on each side of recording to use for drift assessment

%% Misc
LFP_res             = 1000;
channels            = 1:32;

clim_max_percentile = 0.01;

%% Load and preprocess data
close all
load(data_file);

%% Unpack data

spikes              = ephys_data.conditions(1).spikes(channels,:,:);
LFPs                = ephys_data.conditions(1).LFP_trace;

switch stim_type
    case 'whisk'
        stim_onset      = ephys_data.conditions(1).whisk_onset;
    case 'opto'
        stim_onset      = ephys_data.conditions(1).LED_onset;
end

spikes              = spikes - stim_onset; % set whisk onset time to 0
LFP_timestamps      = ([1:size(LFPs,3)]/LFP_res) - stim_onset;

%% Visualisation

figure
set(gcf,'Units','normalized','Position',[0.05 .1 .9 .3])
subplot(1,5,1)
psth(spikes,bin_size,psth_win)
fixplot

subplot(1,5,2)
[im_h, spike_count] = spike_density_plot(spikes,1,psth_win(1):bin_size:psth_win(2));
set(gca,'CLim',[0 robust_max(spike_count(:),clim_max_percentile)])


subplot(1,5,3)
raster_plot(spikes,2)
xlim(psth_win)
fixplot

subplot(1,5,4)
plot_LFP_traces(LFPs,1,LFP_timestamps)
xlim(psth_win)
fixplot

subplot(1,5,5)
plot_LFP_traces(LFPs,1,LFP_timestamps,0.2,0.6)
xlim(psth_win)
fixplot

%% Responsiveness assessment

resp_spikes         = spikes >= spike_resp_win(1) & spikes < spike_resp_win(2);
spont_spikes        = spikes >= spike_spont_win(1) & spikes < spike_spont_win(2);

resp_spikes         = sum(resp_spikes,3);
spont_spikes        = sum(spont_spikes,3);

resp_win_size       = spike_resp_win(2) - spike_resp_win(1);
spont_win_size      = spike_spont_win(2) - spike_spont_win(1);

win_ratio           = resp_win_size / spont_win_size;

corr_spont_spikes   = spont_spikes .* win_ratio;

spont_by_trial      = mean(corr_spont_spikes,1);
mean_spont_spikes   = mean(spont_by_trial);
spont_serr_spikes   = serr(spont_by_trial);

resp_by_trial       = mean(resp_spikes,1);
mean_resp_spikes    = mean(resp_by_trial);
resp_serr_spikes    = serr(resp_by_trial);

% for drift assessment:

early_spike_resps   = resp_by_trial(1:n_trials_for_drift);
late_spike_resps    = resp_by_trial(((end-n_trials_for_drift)+1):end);

median_early        = median(early_spike_resps);
iqr_early           = iqr(early_spike_resps);
median_late         = median(late_spike_resps);

median_diff         = abs(median_early - median_late);
nonparam_zscore     = median_diff / iqr_early;

% for spont vs resp beeswarmplot:
spikes_by_trial     = [spont_by_trial(:); resp_by_trial(:)];
grouping_var        = ones(size(spikes_by_trial));
grouping_var((length(spont_by_trial)+1):end)    = 2;
labels              = {'Spontaneous' 'Evoked'};

[h,p]               = ttest2(resp_by_trial,spont_by_trial);
resp_ratio          = mean(resp_by_trial) / mean(spont_by_trial);
resp_perc_spont     = round(resp_ratio*100);
resp_p          	= round(p,3);

% for baseline drift beeswarmplot
drift_spike_resps   = [early_spike_resps(:) late_spike_resps(:)];
drift_groups        = ones(size(drift_spike_resps));
drift_groups((length(early_spike_resps)+1):end)     = 2;
drift_labels        = {'Early' 'Late'};

%% Responsiveness visualisation

figure
set(gcf,'Units','normalized','Position',[0.05 .1 .9 .3])

subplot(1,4,1)
beeswarmplot(spikes_by_trial,grouping_var,labels)
fixplot
title(['Resp = ' num2str(resp_perc_spont) ' % spont, p = ' num2str(resp_p)])
ylabel('Spike count')
yzero

subplot(1,4,2)
beeswarmplot(drift_spike_resps,drift_groups,drift_labels)
fixplot
yzero
title(['Baseline z-score = ' num2str(nonparam_zscore)])
ylabel('Spike count')

subplot(1,2,2)
plot(resp_by_trial - mean_spont_spikes,'k.','MarkerSize',20)
fixplot
%axis tight
% yzero
xlabel('Trial number')
ylabel('Spike count')
title('Spiking response over time')

