function [timing_data] = timing_function(ephys_data, whisk_resp_win, psth_bins)
% function [timing_data] = timing_function(ephys_data, whisk_resp_win, psth_bins)
% 
% Function to analyse timing experiment data; only takes data where there is 
% a single whisker being stimulated and a single laser / LED power
% 
% EPHYS_DATA - Variable (struct) containing preprocessed data with spike 
% times for the different conditions of the timing experiment
% 
% WHISK_RESP_WIN - window for assessing whisker responses, relative to whisker onset time.
% 
% PSTH_BINS - bin edges for count data
% 
% TIMING_DATA
% struct with fields:
% data_folder             = original data folder
% 
% opto_power              = LED or LASER power for this timing file
% delta_t                 = the time differences between opto and whisk
% 
% spont_spike_rate        = spontaneous spike rate by channel
% spont_spike_std         = spontaneous spike standard deviation by channel
%
% spike_rate              = spike rate in response window by channel for each delta t (spontaneous rate by channel is subtracted)
% spike_std               = standard deviation of spike rate in response window, by channel, for each delta t
%
% norm_spike_rate         = spike rate normalised to the control condition (1 = equal to control spike rate), by channel, for each delta t
% delta_spike_rate        = delta spike rate compared to control condition (condition rate - control rate), by channel, for each delta t
% control_spike_rate      = spike rate by channel in the control condition
% 
% opto_spike_rate         = opto spike rate (minus spontaneous) by channel for each condition
% opto_spike_std          = standard deviation of opto spike rate
% norm_opto_rate          = corr_opto_spike_rate ./ corr_opto_spike_rate(:,end);
% 
% peak_spike_rate         = peak spike rate (minus spontaneous rate)
% norm_peak_spike_rate    = peak spike rate normalised to control (control = 1)
% delta_peak_spike_rate   = delta peak spike rate (Hz) compared to control
% 
% peak_spike_time         = time of peak spike rate, by channel, for each delta t
% diff_peak_spike_time    = delta peak spike time (s) compared to control
% 
% first_spike_time        = time of first spike by channel, for each delta t
% first_spike_jitter      = trial-by-trial jitter (std of first spike time) by channel, for each delta t;
% 
% diff_first_spike_time   = first_spike_jitters - first_spike_jitters(:,end);
% norm_first_spike_jitter = first_spike_jitters ./ first_spike_jitters(:,end);
% 
% psth_counts             = psth_counts;
% spike_density_rates     = spike_density_rates;

%% Input variables

%% Hardcoded variables
peak_rate_kernel    = [0.01]; % peak rate gaussian kernel setting
show_psth           = false; % Flag to not show psth plot from psth function

%% Code execution starts here

if isfield(ephys_data.conditions,'whisk_start')
    whisk_onsets         = [ephys_data.conditions.whisk_start];
elseif isfield(ephys_data.conditions,'whisk_onset')
    whisk_onsets         = [ephys_data.conditions.whisk_onset];
end


if isfield(ephys_data.conditions,'opto_onset')
    opto_onsets         = [ephys_data.conditions.opto_onset];
elseif isfield(ephys_data.conditions,'LED_onset')
    opto_onsets         = [ephys_data.conditions.LED_onset];
end

delta_t             = opto_onsets(:) - whisk_onsets(:);

whisk_time          = median(whisk_onsets);
spont_win           = whisk_resp_win - whisk_time; 

uniq_delta_ts       = unique(delta_t);
n_delta_ts          = length(uniq_delta_ts);

q_control         	= uniq_delta_ts >= 1;    % what is the control condition
control_opto_spikes = ephys_data.conditions(q_control).spikes - opto_onsets(q_control); % get spikes relative to opto pulse in control condition


for k = 1:n_delta_ts
    
    % Select appropriate condition
    q_delta_t       = delta_t == uniq_delta_ts(k);
    
    this_cond       = ephys_data.conditions(q_delta_t);
    
    if isfield(this_cond,'whisk_onset')
        this_t_whisk    = this_cond.whisk_onset;
    elseif isfield(this_cond,'whisk_start')
        this_t_whisk    = this_cond.whisk_start;
    end
    % Retrieve spike data and generate spike count & rate measures
    spikes        	= this_cond.spikes - this_t_whisk;
    
    [spont_spike_rate(:,k), spont_spike_std(:,k)]	= spike_rate_by_channel(spikes, spont_win);   
    [spike_rate(:,k), spike_std(:,k)]               = spike_rate_by_channel(spikes, whisk_resp_win);
    [spike_prob(:,k), n_hits, n_trials(k)]          = spike_prob_by_channel(spikes, whisk_resp_win); 
    
    [peak_spike_rates(1:size(spikes,1),k), peak_spike_times(1:size(spikes,1),k)]        = peak_ROF_by_channel(spikes, whisk_resp_win, peak_rate_kernel);
    [first_spike_times(1:size(spikes,1),k), first_spike_jitters(1:size(spikes,1),k)]	= first_spike_by_channel(spikes, whisk_resp_win);
    
    [~, spike_density_rates(:,:,k)]    = spike_density_plot(spikes,1,psth_bins,show_psth); 
    
    %% Opto resp
    % Get spikes from the window of spontaneous activity in the control condition 
    % corresponding to the whisk window in the current condition
    opto_spikes                                     = control_opto_spikes + uniq_delta_ts(k); % set temporal offset
    [opto_spike_rate(:,k), opto_spike_std(:,k)]     = spike_rate_by_channel(opto_spikes, whisk_resp_win);

    %% do something with psth?
    % Get full post_stimulus_time_histogram
    [plot_handle, psth_counts(:,k), psth_bins]      = psth(spikes, psth_bins, show_psth);
    
end

spont_spike_rate                    = mean(spont_spike_rate,2);
spont_spike_std                     = mean(spont_spike_std,2);

corr_spike_rates                    = spike_rate - spont_spike_rate;
corr_peak_rates                     = peak_spike_rates - spont_spike_rate;
corr_opto_spike_rate                = opto_spike_rate - spont_spike_rate;

%% output variables
timing_data.data_folder             = ephys_data.data_folder;

if isfield(ephys_data,'parameters')
    timing_data.opto_power              = ephys_data.parameters.LED_power;
else
    timing_data.opto_power              = NaN;
end

timing_data.delta_t                 = delta_t';
timing_data.n_trials                = n_trials;

timing_data.spont_spike_rate        = spont_spike_rate;
timing_data.spont_spike_std         = spont_spike_std;
timing_data.spike_rate              = corr_spike_rates;
timing_data.spike_std               = spike_std;
timing_data.norm_spike_rate         = corr_spike_rates ./ corr_spike_rates(:,end);
timing_data.delta_spike_rate        = corr_spike_rates - corr_spike_rates(:,end);
timing_data.control_spike_rate      = corr_spike_rates(:,end);

timing_data.spike_probability   	= spike_prob;

timing_data.opto_spike_rate         = corr_opto_spike_rate;
timing_data.opto_spike_std          = opto_spike_std;
timing_data.norm_opto_rate          = corr_opto_spike_rate ./ corr_opto_spike_rate(:,end);

timing_data.peak_spike_rate         = corr_peak_rates;
timing_data.norm_peak_spike_rate    = corr_peak_rates ./ corr_peak_rates(:,end);
timing_data.delta_peak_spike_rate   = corr_peak_rates - corr_peak_rates(:,end);

timing_data.peak_spike_time         = peak_spike_times;
timing_data.diff_peak_spike_time    = peak_spike_times - peak_spike_times(:,end);

timing_data.first_spike_time        = first_spike_times;
timing_data.first_spike_jitter      = first_spike_jitters;

timing_data.diff_first_spike_time   = first_spike_jitters - first_spike_jitters(:,end);
timing_data.norm_first_spike_jitter = first_spike_jitters ./ first_spike_jitters(:,end);

timing_data.psth_counts             = psth_counts;
timing_data.psth_bins               = psth_bins;
timing_data.spike_density_rates     = spike_density_rates;
