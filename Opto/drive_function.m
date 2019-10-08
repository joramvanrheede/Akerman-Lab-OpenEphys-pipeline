function [drive_data] = drive_function(ephys_data, whisk_resp_win, psth_bins)
% function [DRIVE_DATA] = drive_function(EPHYS_DATA, WHISK_RESP_WIN, PSTH_BINS)
% 
% Extract summary data of the spiking response in the 'Drive' experiment
% 
% Inputs:
% 
% EPHYS_DATA: Variable (struct) containing preprocessed data with spike 
% times for the different conditions of the drive experiment
% 
% WHISK_RESP_WIN: Window for assessing whisker responses, relative to whisker onset time.
% 
% PSTH_BINS: Bin edges for binned spike rate data.
% 
% 
% DRIVE_DATA: a struct with fields:
% 
% whisk_resp_win            : Post-stimulus window for assessing whisker-evoked response
% 
% control_whisk_spike_rates	: Whisker-evoked spike rates in control condition, by channel
% opto_whisk_spike_rates   	: Whisker-evoked spike rates in test condition, corrected for control opto response in the same window, by channel
% delta_spike_rates        	: Difference in whisker-evoked spike rate in test condition relative to control, by channel
% 
% whisk_spike_rates       	: Whisker-evoked spike rates per channel, corrected for spontaneous rates by channel, for each condition
% whisk_spike_stds        	: Standard deviation of whisker-evoked spike rates by channel, for each condition
% 
% spont_spike_rates      	: Spontaneous spike rate per channel in a window of equal size to the whisker resp win
% spont_spike_stds        	: Standard deviation of spontaneous spike rate per channel;
% 
% opto_spike_rates         	: Optogenetically evoked spike rate in equivalent whisker response window in control condition;
% opto_spike_stds         	: Standard deviation of opto response in control condition;
% 
% peak_spike_rates       	: Peak whisker-evoked spike rate by channel for each condition
% 
% peak_spike_times        	: Time of peak whisker-evoked spike rate by channel for each condition
% diff_peak_spike_times   	: Difference in peak whisker-evoked spike rate in test condition relative to control, by channel
% 
% first_spike_times        	: Average time of first spike in whisker response window per channel for each condition
% diff_first_spike_times   	: Difference in spike time in opto condition relative to control
% 
% whisk_spike_probabilities	: Estimated spike probability in whisk_resp_win by channel (probability of at least one spike)
% whisk_n_hits            	: Number of trials with at least one spike, per channel;
% n_trials                	: Total number of trials for each condition;
% 
% opto_spike_probabilities	: opto_spike_probs;
% opto_n_hits            	: opto_n_hits;
% 
% psth_bins                 : Bins for assessing the binned spike rates by channel
% spike_density_rates    	: Binned spike rates over time bins by channel for each condition
% 


% hardcoded peak spike rate kernel
rate_kernel_win         = [0.01];           % 3SD gaussian kernel to aggregate spikes for peak spike rate

% isolate conditions data for ease of writing / reading the below code:
cond_data               = ephys_data.conditions;

% Find the time window in the opto response in the test condition
opto_whisk_win_offset   = cond_data(1).whisk_onset - cond_data(1).LED_onset; % First condition should be opto and whisk
opto_resp_win           = cond_data(2).LED_onset + opto_whisk_win_offset + whisk_resp_win;


for a = 1:length(cond_data)
    
    % Get spike times relative to whisk
    spikes                  = cond_data(a).spikes - cond_data(a).whisk_onset; 
    
    % Spontaneous window relative to start of trial
    spont_resp_win          = whisk_resp_win - cond_data(a).whisk_onset;

    % Quantify spikes:
    
    % Whisker-evoked spike rate in target response window
    [whisk_spike_rates(:,a), whisk_spike_stds(:,a)]     = spike_rate_by_channel(spikes, whisk_resp_win);
    
    % Spontaneous spikes (in window of equal size to whisk_resp_win but counting from start of trial)
    [spont_spike_rates(:,a), spont_spike_stds(:,a)]     = spike_rate_by_channel(spikes, spont_resp_win);
    
    % Peak spike rate in response to whisker stimulus
    [peak_spike_rates(:,a), peak_spike_times(:,a)]      = peak_ROF_by_channel(spikes, whisk_resp_win, rate_kernel_win);
    
    % First spike time in response to the whisker stimulus for each channel
    [first_spike_times(:,a), first_spike_jitter(:,a)]  	= first_spike_by_channel(spikes, whisk_resp_win);

    % Spike probability by channel
    [spike_probs(:,a), n_hits(:,a), n_trials(a)]    	= spike_prob_by_channel(spikes, whisk_resp_win);
    
    % get spike rate / density data for all channels but don't make the plot
    [~, spike_density_rates(:,:,a)]                     = spike_density_plot(spikes,1,psth_bins, false);
    
    
    % Condition 2 = control condition (separate whisker and opto stimulus)
    if a == 2
     	% Spike rate in opto-only stimulus window that corresponds to the time of 
        % the whisker stimulus in the test condition
        [opto_spike_rates, opto_spike_stds]             = spike_rate_by_channel(spikes, opto_resp_win);
    
        % Spike probability in opto-only stimulus window that corresponds to the time of 
        % the whisker stimulus in the test condition
        [opto_spike_probs, opto_n_hits, ~]              = spike_prob_by_channel(spikes, opto_resp_win);
    end
    
    
end


%% Make some corrections to response measures

spont_spike_rates               = mean(spont_spike_rates,2);
spont_spike_stds                = mean(spont_spike_stds,2);

% correct both for spontaneous
whisk_spike_rates               = whisk_spike_rates - spont_spike_rates; % subtract mean of spontaneous; take mean of spontaneous over both conditions for more robust measure
opto_spike_rates                = opto_spike_rates - spont_spike_rates;
peak_spike_rates                = peak_spike_rates - spont_spike_rates;

% Time difference measures
diff_peak_spike_times           = peak_spike_times(:,1) - peak_spike_times(:,2);
diff_first_spike_times          = first_spike_times(:,1) - first_spike_times(:,2);


% Correct for the expected number of spikes given the opto stimulation
control_whisk_spike_rates       = whisk_spike_rates(:,2);
opto_whisk_spike_rates          = whisk_spike_rates(:,1) - opto_spike_rates;

% Get difference in spike rate between +opto condition and control for a 
% one-data-point summary of the result of the experiment
delta_spike_rates               = opto_whisk_spike_rates - control_whisk_spike_rates;


%% Create output variable
drive_data.whisk_resp_win           	= whisk_resp_win;

drive_data.control_whisk_spike_rates    = control_whisk_spike_rates;
drive_data.opto_whisk_spike_rates       = opto_whisk_spike_rates;
drive_data.delta_spike_rates            = delta_spike_rates;

drive_data.whisk_spike_rates            = whisk_spike_rates;
drive_data.whisk_spike_stds             = whisk_spike_stds;

drive_data.spont_spike_rates            = spont_spike_rates;
drive_data.spont_spike_stds             = spont_spike_stds;

drive_data.opto_spike_rates             = opto_spike_rates;
drive_data.opto_spike_stds              = opto_spike_stds;

drive_data.peak_spike_rates             = peak_spike_rates;

drive_data.peak_spike_times             = peak_spike_times;
drive_data.diff_peak_spike_times        = diff_peak_spike_times;

drive_data.first_spike_times            = first_spike_times;
drive_data.diff_first_spike_times       = diff_first_spike_times;

drive_data.whisk_spike_probabilities  	= spike_probs;
drive_data.whisk_n_hits               	= n_hits;
drive_data.n_trials                     = n_trials;

drive_data.opto_spike_probabilities     = opto_spike_probs;
drive_data.opto_n_hits                  = opto_n_hits;

drive_data.psth_bins                    = psth_bins;
drive_data.spike_density_rates          = spike_density_rates;

