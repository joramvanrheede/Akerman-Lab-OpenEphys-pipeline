

%% Pre stim examples: 
% '/Volumes/Akermanlab/Joram/Cortex_20x100Hz_Plasticity/Test_1_pre_0/2019_02_04/2019_02_04-3-Test_1_pre_0.mat'

%% RWS 100Hz examples: 
% '/Volumes/Akermanlab/Joram/Cortex_20x100Hz_Plasticity/RWS_1/2019_02_12/2019_02_12-2-RWS_1.mat'

%% RWS 8Hz examples: 
% '/Volumes/Akermanlab/Joram/Plasticity_8Hz_PW/RWS_1/2019_03_09/2019_03_09-4-RWS_1.mat'
% '/Volumes/Akermanlab/Joram/Plasticity_8Hz_PW/RWS_1/2019_03_10/2019_03_10-3-RWS_1.mat'

%% VPM examples:
% '/Volumes/Akermanlab/Joram/Thalamus/RWS_1/2019_01_30/2019_01_30-2-RWS_1.mat'

%% PoM examples: 
% '/Volumes/Akermanlab/Joram/Thalamus/RWS_1/2019_01_26/2019_01_26-6-RWS_1.mat'

%% RJB PV-Cre
% '/Volumes/Akermanlab/Joram/RJB/Laser_pulse/2019_01_24/2019_01_24-1-Laser_pulse.mat'
% '/Volumes/Akermanlab/Joram/RJB/Drive/2019_01_24/2019_01_24-2-Drive.mat'
% '/Volumes/Akermanlab/Joram/RJB/Test_post_1/2019_02_11/2019_02_11-2-Test_post_1.mat'
% '/Volumes/Akermanlab/Joram/RJB/Test_post_seizure/2019_02_11/2019_02_11-5-Test_post_seizure.mat'

%% User set variables:

data_file       = '/Volumes/Akermanlab/Joram/Cortex_20x100Hz_Plasticity/Test_1_pre_0/2019_02_04/2019_02_04-3-Test_1_pre_0.mat';

condition_nr    = 1; % Which condition nr to make plots for? Use 1 for data with only a single condition

split_rows      = 'trials'; % 'trials' or 'channels'. 'Trials' has one row for each trial and averages over channels, 'Channels' does vice versa

channels        = [1:32]; % Which channels to include
trials          = [1:30]; % Which trial numbers to include

psth_win        = [-0.1 0.2]; % sets x-axis values for all plots
psth_offset     = 1; % set this time to zero (e.g. stimulus time)
binsize         = [0.001]; % bin size for psth

fig_title       = 'Whisker stimulation';


%% Shaded regions in graphs?
do_shade1       = true;
shade1_colour   = 'r';
shade1_alpha    = 0.5;
shade1_xvals    = [0 0.005];

do_shade2       = false;
shade2_colour   = 'r';
shade2_alpha    = 0.7;
shade2_xvals    = [0 0.005];

%% 
RWS_continuous 	= false; % for continuous RWS where data are stored as a single very long trial; this triggers the generation of a number of 'sweeps' to show multiple repeats of a continuously repeating stimulus on the same row
n_reps          = 8; % for repeating stimuli, how many repeats to show in one 'sweep'
do_offset_corr  = false; % if there is an accumulating offset between repeats this allows for correction
offset_val      = 0.0002; % value for offset



%% Saving options
save_fig    = true;
save_dir    = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/RJB/';
fig_format  = '-depsc'; % '-depsc' for vector graphics or '-dpng' for png image file

%% Axis labels; Y-axis for raster is given by split_rows

x_ax_label  = 'Time (s)';

%% Code execution starts here

close all

load(data_file); % This loads 'ephys_data' struct

spikes      = ephys_data.conditions(condition_nr).spikes; % get spikes from the relevant condition

%% new for continuous 8Hz RWS - generate a number of fake 'sweeps' to show repetitive nature of stimulus
if RWS_continuous
    
    if do_offset_corr
        for a = 1:size(spikes,2)
            spikes(:,a,:) = spikes(:,a,:) - offset_val*a; % offset correction
        end
    end
    
    new_spikes  = [];
    for a = 1:size(spikes,2)/n_reps
        spike_reps = [];
        for b = 1:n_reps
            spike_reps = cat(3,spike_reps,spikes(:,(a-1)*b+b,:)+(b-1)*0.125);
        end
        new_spikes = [new_spikes spike_reps];
    end
    spikes = new_spikes;
end

%% Rasterplot 

switch split_rows
    case 'trials'
        split_dim = 2;
        y_ax_label  = 'Trials';
    case 'channels'
        split_dim = 1;
        y_ax_label  = 'Channels';
end

%%

figure(1)

rasterplot_joram(spikes(channels,trials,:)-psth_offset,split_dim);
xlim(psth_win)
title([fig_title ' raster plot'])
ylabel(y_ax_label)
xlabel(x_ax_label)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold')

if do_shade1
    shaded_region(shade1_xvals, shade1_colour, shade1_alpha);
end
if do_shade2
    shaded_region(shade2_xvals, shade2_colour, shade2_alpha);
end

if save_fig
    print(gcf,[save_dir 'Rasterplot_' fig_title],fig_format)
end

%% PSTH

figure(2)

target_spikes 	= spikes(channels,trials,:); % get relevant spikes
psth_bins     	= (psth_win(1):binsize:psth_win(2)); % set bin edges for psth counts
psth_counts    	= histc(target_spikes(:)-psth_offset,psth_bins); % count spikes within bins

bar(psth_bins,psth_counts,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]) % plot the spike counts vs time bins

% plot aesthetics:
xlim(psth_win);
title([fig_title ' PSTH'])
xlabel('Post-stimulus time (s)')
ylabel('Spike count')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold')

if do_shade1
    shaded_region(shade1_xvals, shade1_colour, shade1_alpha);
end

if do_shade2
    shaded_region(shade2_xvals, shade2_colour, shade2_alpha);
end

if save_fig
    print(gcf,[save_dir 'PSTH_' fig_title],fig_format)
end

%% LFP average across trials / channels, aligned with rasterplot idea

figure(3)

LFP_traces      = ephys_data.conditions(condition_nr).LFP_trace(channels,trials,:);
LFP_timepoints  = (1:size(LFP_traces,3))/1000 - psth_offset;
q_LFP           = LFP_timepoints > psth_win(1) & LFP_timepoints < psth_win(2);

switch split_rows
    case 'trials'
        mean_LFP_traces     = squeeze(mean(LFP_traces,1));
    case 'channels'
        mean_LFP_traces     = squeeze(mean(LFP_traces,2));
end

mean_LFP_traces = notch_filt(mean_LFP_traces',1000,50)'; % remove 50Hz noise

plot_LFPs_by_channel(mean_LFP_traces(:,q_LFP),LFP_timepoints(q_LFP),.5)

title([fig_title ' LFP by ' split_rows])
axis tight
xlabel('Post-stimulus time (s)')
ylabel(['LFP ' split_rows])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold')

if do_shade1
    shaded_region(shade1_xvals, shade1_colour, shade1_alpha);
end

if do_shade2
    shaded_region(shade2_xvals, shade2_colour, shade2_alpha);
end

if save_fig
    print(gcf,[save_dir 'LFP_' fig_title],fig_format)
end


%% LFP

figure(4)

LFP_traces      = ephys_data.conditions(condition_nr).LFP_trace(channels,trials,:);
LFP_timepoints  = (1:size(LFP_traces,3))/1000 - psth_offset;
q_LFP           = LFP_timepoints > psth_win(1) & LFP_timepoints < psth_win(2);

switch split_rows
    case 'Trials'
        mean_LFP_traces     = squeeze(mean(LFP_traces,1));
    case 'Channels'
        mean_LFP_traces     = squeeze(mean(LFP_traces,2));
end

mean_LFP_traces     = notch_filt(mean_LFP_traces',1000,50)'; % remove 50Hz noise

overall_LFP_trace   = mean(mean_LFP_traces);

plot(LFP_timepoints(q_LFP),overall_LFP_trace(q_LFP),'k-','LineWidth',2)
title([fig_title ' LFP'])
axis tight
xlabel('Post-stimulus time (s)')
ylabel('LFP (uV)')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold')

if do_shade1
    shaded_region(shade1_xvals, shade1_colour, shade1_alpha);
end

if do_shade2
    shaded_region(shade2_xvals, shade2_colour, shade2_alpha);
end

if save_fig
    print(gcf,[save_dir 'LFP_' fig_title],fig_format)
end


