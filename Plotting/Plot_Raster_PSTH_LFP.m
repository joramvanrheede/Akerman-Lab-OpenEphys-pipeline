% Plot_Raster_PSTH_LFP
% 
% This script will generate a number of plots for inspecting your data.
% For a single condition in ephys_data, it will generate a raster plot (by
% layer or by trial), a post-stimulus time histogram, an LFP plot (again, by 
% layer or by trial) and an overall LFP trace.
% 
% You need to set a number of variables at the top of the script before it
% will run.
% 
% Suggested use: When you want to make figures for a particular experiment,
% copy this script to an ad hoc folder and change the variables as needed.
% 

%% User set variables:

%% Change these variables for every experiment:
data_file       = 'Full/Path/To/Processed/Data_file.mat'; % Full path to data file, preprocessed into trial-synchronised spike time data using ephys_metadata_reader
fig_title       = 'Experiment title here'; % This sets the title of the figures - NOTE: also used for the filename of saved figures

% This determines the selection of your data:
condition_nr    = 1;        % Which condition nr to make plots for? Use 1 for data with only a single condition
channels        = [1:32];   % Which channels to include
trials          = [1:10];   % Which trial numbers to include. Use a single trial to plot a 'canonical' raster plot

% How much of the data to show and what time bins to use?
psth_win        = [-.5 1];  % sets x-axis values for all plots
psth_offset     = 1;        % set this time to zero (e.g. stimulus time)
binsize         = [0.005];  % bin size for psth


%% Change the following variables as necessary:

split_rows      = 'channels'; % 'trials' or 'channels'. 'Trials' has one row for each trial and averages over channels, 'Channels' does vice versa

% Saving options for output figures (Note - will overwrite files of the same name if script is re-run) 
save_fig        = false;
save_dir        = 'Path/To/Figure/Saving/Directory';
fig_format      = '-depsc'; % '-depsc' for vector graphics or '-dpng' for png image file


% Shaded regions in graphs? Currently 2 shades are supported but it should be obvious
% from this script how to add additional shaded regions if required:

do_shade1       = false;    % 
shade1_colour   = 'r';
shade1_alpha    = 0.5;
shade1_xvals    = [0 0.005]; % time window for shaded region 1

do_shade2       = false;
shade2_colour   = 'r';
shade2_alpha    = 0.7;
shade2_xvals    = [0 0.005]; % time window for shaded region 2


% Only for continuous RWS where data are stored as a single very long trial; 
% this triggers the generation of a number of 'sweeps' to show multiple 
% repeats of a continuously repeating stimulus on the same row
RWS_continuous 	= false;    % leave to false unless the above is true
n_reps          = 8;        % for repeating stimuli, how many repeats to show in one 'sweep'
do_offset_corr  = false;    % if there is an accumulating offset between repeats this allows for correction
offset_val      = 0.0002;   % value for offset (correcting for drift caused by bug in PulsePal protocol)

% Axis labels; Y-axis for raster is given by split_rows
x_ax_label  = 'Time (s)';

%% Code execution starts here

close all

load(data_file); % This loads 'ephys_data' struct

spikes  = ephys_data.conditions(condition_nr).spikes; % get spikes from the relevant condition

%% new for continuous 8Hz RWS - generate a number of fake 'sweeps' to show repetitive nature of stimulus
if RWS_continuous
    
    if do_offset_corr
        for a = 1:size(spikes,2)
            spikes(:,a,:) = spikes(:,a,:) - offset_val*a; % offset correction for repeated whisks within same trial (due to bug in PulsePal protocol)
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

%% Rasterplot

figure(1)

raster_plot(spikes(channels,trials,:)-psth_offset,split_dim);
xlim(psth_win)
title([fig_title ' raster plot'])
ylabel(y_ax_label)
xlabel(x_ax_label)
fixplot

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

psth(target_spikes - psth_offset,binsize,psth_win)

% plot aesthetics:
xlim(psth_win);
title([fig_title ' PSTH'])
fixplot

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


