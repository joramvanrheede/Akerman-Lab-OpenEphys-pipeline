function drive_data = drive_plots_by_unit(ephys_data, units, resp_win, psth_bins, artifact_win)
% function drive_data = drive_plots_by_unit(EPHYS_DATA),
% function drive_data = drive_plots_by_unit(EPHYS_DATA, UNITS, RESP_WIN, PSTH_BINS, ARTIFACT_WIN)
% 
% Visualise data and extract key measures from 'drive' experiment, in which
% a whisker stimulus happens either during (condition 1) or separate from
% a continuous optogenetic stimulus. Geared towards single unit analysis.
% 
% This function will plot the following:
% 
% Figure 1: A raster plot for both conditions:
% Figure 2: A PSTH for both conditions:
% Figure 3: A spike density plot for both conditions:
% Figure 4: A raster plot, PSTH and spike density plot of the opto response
% Figure 5: A comparison of binned spike rate, peak spike rate and peak spike
%           time for both conditions including p-values
% 
% REQUIRED INPUTS:
% 
% EPHYS_DATA: EPHYS_DATA structure as generated by preprocess_multiunit, containing
% spike data and metadata for a 'Drive' type experiment.
% 
% OPTIONAL:
% 
% UNITS:        Specify which channels, e.g. [1:8]; Default is ':' (all). 
% RESP_WIN:     Window for assessing spiking response; Default is [0.007 0.030].
% PSTH_BINS:    Bins for PSTH and spike density plots; Default is [-0.1:0.001:0.3]
% ARTIFACT_WIN: To remove potential artifacts, spike times in this window are 
%               set to NaN. Default is [-0.001 0.006].
% 
% OUTPUT: a data structure DRIVE_DATA with fields:
%

% Default to all channels / units
if nargin < 2 || isempty(units)
    units        = ':';
end

% Default to resp win from 6ms (after any artifacts) to 30ms (should capture
% most of the direct stimulus-driven activity
if nargin < 3 || isempty(resp_win)
    resp_win        = [0.006 0.030];
end

% Default PSTH range; 300ms post stimulus should capture even long-tailed responses
if nargin < 4 || isempty(psth_bins)
    psth_bins       = [-0.1:0.001:0.3];
end

% Set any spikes during this window to NaN; -0.001 to 0.006 is where any piezo artifacts
% may occur
if nargin < 5
    artifact_win    = [-0.001 0.006];
end

% Hardcoded for now:
cond_names          = {'Whisk + opto' 'Whisk only'};    % Condition nicknames
rate_kernel_win     = [0.01];	% 3SD gaussian kernel to aggregate spikes for peak spike rate
clim_perc           = 0.5;      % Percentile for setting max color limit;
opto_psth_margin    = 0.3;      % margin (in seconds) around opto stimulation to include in opto PSTHs

%% Running code starts here

cond_data             	= ephys_data.conditions;

% Find the time window in the opto response
opto_whisk_win_offset   = cond_data(1).whisk_onset - cond_data(1).LED_onset; % First condition should be opto and whisk
opto_resp_win           = resp_win + opto_whisk_win_offset;

opto_duration           = cond_data(2).LED_duration;
bin_size                = mean(diff(psth_bins));
opto_psth_bins          = [(-opto_psth_margin):bin_size:(opto_duration+opto_psth_margin)];

plot_lims               = [psth_bins(1) psth_bins(end)];
opto_plot_lims          = [opto_psth_bins(1) opto_psth_bins(end)];

% Pre-create figures of the right size
raster_fig_h            = figure;
set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
psth_fig_h              = figure;
set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
density_fig_h           = figure;
set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
opto_fig_h              = figure;
set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
        
for a = 1:length(cond_data)
    
    spikes                          = cond_data(a).spikes(units,:,:) - cond_data(a).whisk_onset;
    
    q_artifact                      = spikes > artifact_win(1) & spikes < artifact_win(2);
    spikes(q_artifact)              = NaN;
    
    % Quantify spikes
    all_spike_rates(:,:,a)          = spike_rates_individual(spikes, resp_win);
    [spike_probs(:,a), n_hits(:,a), n_trials(:,a)]      = spike_prob_by_channel(spikes, resp_win);
    [first_spike_times(:,a), first_spike_jitters(:,a)]  = first_spike_by_channel(spikes, resp_win);
    
    [peak_spike_rates(:,a), peak_spike_times(:,a)]         	= peak_ROF_by_channel(spikes, resp_win, rate_kernel_win);
    
    first_spikes_all(:,:,a)         = first_spike_individual(spikes, resp_win);
    
    % Quantify opto spikes
    opto_spikes                     = cond_data(a).spikes(units,:,:) - cond_data(a).LED_onset;
    all_opto_spike_rates(:,:,a)     = spike_rates_individual(opto_spikes, opto_resp_win);
    
    % Quantify spikes in spontaneous window
    spont_spikes                    = cond_data(a).spikes(units,:,:);
    all_spont_spike_rates(:,:,a)	= spike_rates_individual(spont_spikes, opto_resp_win);
    
    
    %% Make plots
    figure(raster_fig_h)
    
    % Raster plot
    subplot(1,2,a)
    raster_plot(spikes,1)
    xlim(plot_lims)
    ylabel('Unit number')
    xlabel('Time (s)')
    fixplot
    title(cond_names{a})
    
    % PSTH
    figure(psth_fig_h)
    subplot(1,2,a)
    [psth_handle psth_counts(:,a)] = psth(spikes,psth_bins);
    title(cond_names{a})
    fixplot
    
    % Spike density plot
    figure(density_fig_h)
    subplot(1,2,a)
    [image_handle, spike_density_counts(:,:,a)] = spike_density_plot(spikes,1,psth_bins);
    ylabel('Unit number')
    xlabel('Time (s)')
    color_lims = [0 robust_max(spike_density_counts(:),clim_perc)];
    set(gca,'CLim',color_lims)
    colorbar
    title([cond_names{a} ' spike rate'])
    
    if a == 2 % Control condition; make plots for opto only as well
        
        figure(opto_fig_h)
        
        % Opto raster
        subplot(1,3,1)
        raster_plot(opto_spikes,1)
        xlim(opto_plot_lims)
        ylabel('Unit number')
        xlabel('Time (s)')
        fixplot
        title('Opto only raster plot')
        
        % Opto PSTH
        subplot(1,3,2)
        [opto_psth_handle, opto_psth_counts] = psth(opto_spikes,opto_psth_bins);
        title('Opto only PSTH')
        fixplot
        
        % Opto spike density plot
        subplot(1,3,3)
        [image_handle, opto_spike_density_counts] = spike_density_plot(opto_spikes,1,opto_psth_bins);
        ylabel('Unit number')
        xlabel('Time (s)')
        color_lims = [0 robust_max(opto_spike_density_counts(:),clim_perc)];
        set(gca,'CLim',color_lims)
        colorbar
        title(['Opto only spike density'])
        
    end  
end

% Make sure that y-axes are the same in PSTH subplots in the same figure
figure(psth_fig_h)
subplot_equal_y

% Make sure that colour scaling is the same in spike density subplots in the same figure 
figure(density_fig_h)
subplot_equal_clims

%% Corrections for opto background and for spontaneous

% correct for opto background
all_spike_rates(:,:,1)    = all_spike_rates(:,:,1) - mean(all_opto_spike_rates(:,:,2),2);

% correct for spontaneous
all_spike_rates(:,:,2)    = all_spike_rates(:,:,2) - mean(all_spont_spike_rates(:,:,2),2);


%% Plotting:

% Create figure of the right size
mean_plots_h            = figure;
set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')

%% Beeswarmplot spike rates

lineplot_spike_rates    = squeeze(mean(all_spike_rates,2));

subplot(1,3,1)
scatter(lineplot_spike_rates(:,2),lineplot_spike_rates(:,1))
hold on
xlim(ylim)
plot(xlim,xlim);
fixplot
xlabel('Spike rate control')
ylabel('Spike rate opto')




%% Beeswarmplot peak spike rate



subplot(1,3,2)
scatter(peak_spike_rates(:,2),peak_spike_rates(:,1))
hold on
xlim(ylim)
plot(xlim,xlim);
fixplot
xlabel('Peak rate control')
ylabel('Peak rate opto')

% [peak_rate_h,peak_rate_p]               = ttest2(peak_spike_rates(:,1)',peak_spike_rates(:,2)');


%% Beeswarmplot peak spike time


subplot(1,3,3)
set(gcf,'PaperPositionMode','auto')
scatter(peak_spike_times(:,2),peak_spike_times(:,1))
hold on
xlim(ylim)
plot(xlim,xlim);
fixplot
xlabel('Peak time control')
ylabel('Peak time opto')

% [first_time_h,first_time_p]               = ttest2(first_spike_times(:,1)',first_spike_times(:,2)');

%% 
figure
set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')

subplot(1,2,1)
set(gcf,'PaperPositionMode','auto')
scatter(spike_probs(:,2),spike_probs(:,1))
hold on
xlim(ylim)
plot(xlim,xlim);
fixplot
xlabel('Spike probability control')
ylabel('Spike probability opto')


subplot(1,2,2)
set(gcf,'PaperPositionMode','auto')
scatter(first_spike_times(:,2),first_spike_times(:,1))
hold on
xlim(ylim)
plot(xlim,xlim);
fixplot
xlabel('First spike control')
ylabel('First spike opto')

%%

[spike_rate_h,spike_rate_p]                 = ttest2(all_spike_rates(:,:,1)',all_spike_rates(:,:,2)');

[first_spike_h, first_spike_p]              = ttest2(first_spikes_all(:,:,1)',first_spikes_all(:,:,2)');

keyboard
spike_prob_hits     = n_hits;
spike_prob_misses   = n_trials - spike_prob_hits;

for a = 1:size(spike_prob_hits,1)
    control_row             = [spike_prob_hits(a,2) spike_prob_misses(a,2)];
    test_row                = [spike_prob_hits(a,1) spike_prob_misses(a,1)];
    
    fisher_tab              = int32([control_row; test_row]);
    
    [h, spike_prob_p(a)]    = fishertest(fisher_tab);
end


%% Generate output data structure


drive_data.data_folder          = ephys_data.data_folder;
drive_data.channels             = units;
drive_data.resp_win             = resp_win;
drive_data.psth_bins            = psth_bins;
drive_data.opto_psth_bins       = opto_psth_bins;
drive_data.bin_size             = bin_size;
drive_data.artifact_win         = artifact_win;

drive_data.binned_rates         = all_spike_rates;
drive_data.spike_rate_p         = spike_rate_p;

drive_data.spike_probs          = spike_probs;
drive_data.spike_prob_p         = spike_prob_p;

drive_data.first_spike_times    = first_spike_times;
drive_data.first_spikes_all     = first_spikes_all;
drive_data.first_spike_p        = first_spike_p;

drive_data.psths                = psth_counts;
drive_data.spike_density        = spike_density_counts;

drive_data.opto_psths           = opto_psth_counts;
drive_data.opto_spike_density   = opto_spike_density_counts;




