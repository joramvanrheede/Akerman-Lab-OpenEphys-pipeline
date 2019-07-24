

%% Input variables 
data_file           = ['/Volumes/Akermanlab/Joram/Preprocessed data/AVK S1/Timing/2019_05_21/2019_05_21-22-Timing.mat'];

q_reload            = 1;

%% PSTH settings
early_resp_win      = [0.0055 0.020]; % response window for assessing 'early' response (spike count is from this time window)
late_resp_win       = [0.020 0.2];   % response window for assessing 'late' response (spike count is from this time window)
opto_resp_win       = [0 0.03];     % window for assessing opto response

psth_win            = [-0.1 .3];   % window for plotting PSTH
psth_bin_size       = [0.001];  % bin size for PSTH

channels            = [1:32];    % which recording sites / channels to include


%% Only necessary if auto-saving figures
save_figs           = true;
fig_save_dir        = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/AVK/Timing';

early_save_name     = 'Early spikes';
late_save_name      = 'Late spikes';
psth_save_name      = 'PSTHs';
peak_ROF_save_name  = 'Peak ROF';

fig_format          = '-dpng'; % '-dpng' for png bitmap file, '-depsc' for .eps vector graphics file

%% Code execution starts here

close all

n_channels      = length(channels);

if q_reload
    load(data_file)
end

early_spike_count   = [];
late_spike_count    = [];
serr_early_spikes   = [];
serr_late_spikes    = [];
peak_spike_rate     = [];
opto_spike_count    = [];
delta_t             = [];
whisker_nr          = [];
opto_power          = [];
n_trials            = [];
all_psth_counts    	= [];
counter             = 0;
for a = 1:length(ephys_data.conditions)
    
    this_cond               = ephys_data.conditions(a);
    
    % for compatibility with multiple-stimulator files; will do stim 1 only.
    if this_cond.whisk_stimulator ~= 1
        continue
    end
    
    counter                         = counter + 1;
    this_t_whisk                    = this_cond.whisk_onset;
    this_t_opto                     = this_cond.LED_onset;
    this_whisker_nr                 = this_cond.whisk_stimulator;
    this_n_trials                   = this_cond.n_trials;
    
    spikes                          = this_cond.spikes(channels, :, :) - this_t_whisk;
    
    delta_t(counter)                = this_t_opto - this_t_whisk;
    whisker_nr(counter)             = this_whisker_nr;
    opto_power(counter)             = this_cond.LED_power;
    
    [early_spike_count(counter), serr_early_spikes(counter)]    = spike_rate_in_win(spikes, early_resp_win);
    [late_spike_count(counter), serr_late_spikes(counter)]      = spike_rate_in_win(spikes, late_resp_win);
    
    peak_spike_rate(counter)        = peak_ROF_in_win(spikes, early_resp_win);
    
    n_trials(counter)           	= this_n_trials;
    
    % get full post_stimulus_time_histogram
    figure
    [plot_handle, psth_counts, psth_bins]  = psth(spikes, psth_bin_size, psth_win);
    close(gcf)
    
    all_psth_counts                 = [all_psth_counts; psth_counts(:)'];
    
    %% opto resp
    
    spikes                          = (spikes + this_t_whisk) - this_t_opto;
    spikes_in_opto_win              = spikes > opto_resp_win(1) & spikes <= opto_resp_win(2);
    opto_spike_count(counter)       = sum(spikes_in_opto_win(:)) / n_channels / this_n_trials;
end

q_control         	= delta_t > 1;
q_below_0       	= delta_t <= 0; % Select only instances where the opto stimulus happened before the whisk

%% response plotting
figure(1)
errorbar(delta_t(q_below_0),early_spike_count(q_below_0),serr_early_spikes(q_below_0),'k.-','LineWidth',2,'MarkerSize',20)
xlimits     = xlim;
line([xlimits],[early_spike_count(end) early_spike_count(end)],'Color',[1 0 0],'LineWidth',2,'LineStyle',':')
title('Early response window')
ylabel('Mean spike rate')
xlabel('Opto-whisk time delay')
fixplot
yzero

figure(2)
errorbar(delta_t(q_below_0),late_spike_count(q_below_0),serr_late_spikes(q_below_0),'k.-','LineWidth',2,'MarkerSize',20)
xlimits     = xlim;
line([xlimits],[late_spike_count(end) late_spike_count(end)],'Color',[1 0 0],'LineWidth',2,'LineStyle',':')
title('Late response window')
ylabel('Mean spike rate')
xlabel('Opto-whisk time delay')
fixplot
yzero

figure(3)
plot(delta_t(q_below_0),peak_spike_rate(q_below_0),'k.-','LineWidth',2,'MarkerSize',20)
xlimits     = xlim;
line([xlimits],[peak_spike_rate(end) peak_spike_rate(end)],'Color',[1 0 0],'LineWidth',2,'LineStyle',':')
title('Peak spike rate (early window)')
ylabel('Peak spike rate')
xlabel('Opto-whisk time delay')
fixplot
yzero

%% Delta-t PSTH plots

figure;
n_delta_ts          = length(delta_t);
plot_ymax           = [];
for k = 1:n_delta_ts
    
    % What delta_t are we processing
    this_delta_t    = delta_t(k);
    q_delta_t       = delta_t == this_delta_t;

    % Select relevant PSTHs
    this_psth       = mean(all_psth_counts(q_delta_t,:),1);   %
    this_n_trials   = n_trials(q_delta_t);
    
    % Convert to spike rate per channel per trial
    plot_psth       = this_psth / n_channels / this_n_trials / psth_bin_size; 
    
    % Determine maximum response in whisk window (purely for setting Y axes later)
    psth_bins       = psth_win(1):psth_bin_size:psth_win(2);
    q_resp          = psth_bins > early_resp_win(1) & psth_bins < late_resp_win(2);
    plot_ymax(k)  	= max(plot_psth(q_resp));
    
    % Actual plotting
    subplot(n_delta_ts, 1, k)
    plot(psth_bins(1:end-1),plot_psth,'k-','LineWidth',2)
    
    % Shade response regions
    shaded_region(early_resp_win,'g',0.4)
    shaded_region(late_resp_win,'r',0.4)
    
    % Set sensible axis and labels
    xlim(psth_win)
    title(['Opto-whisk delay = ' num2str(-this_delta_t * 1000) 'ms'])
    ylabel('Spike rate (Hz)')
    
end
% set all y axes to equal
subplot_equal_y(robust_max(plot_ymax,20,'all') * 1.2);
set(gcf,'Units','normalized','Position',[.3 0 .4 1])


% Save figures?
if save_figs
    [file_path, file_name, file_ext]   = fileparts(data_file); % unpack data_file_name to generate figure file name
    
    % Generate file names for figures
    early_spike_filenm  = fullfile(fig_save_dir,[file_name ' ' early_save_name]); 
    late_spike_filenm   = fullfile(fig_save_dir,[file_name ' ' late_save_name]); 
    peak_ROF_filenm     = fullfile(fig_save_dir,[file_name ' ' peak_ROF_save_name]); 
    psth_filenm         = fullfile(fig_save_dir,[file_name ' ' psth_save_name]); 
    
    
    % Select figures and save them
    figure(1)
    print(early_spike_filenm,fig_format,'-painters')
    figure(2)
    print(late_spike_filenm,fig_format,'-painters')
    figure(3)
    print(peak_ROF_filenm,fig_format,'-painters')
    figure(4)
    print(psth_filenm,fig_format,'-painters')
end


