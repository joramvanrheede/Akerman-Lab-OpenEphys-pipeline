function [drive_multi_data] = drive_multi_function(ephys_data, whisk_resp_win, psth_bins, n_stims)
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
% opto_power                : What was the power of the LED / Laser?
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

if nargin < 4
    n_stims = 4;
end

% isolate conditions data for ease of writing / reading the below code:
cond_data               = ephys_data.conditions;

% Get whisk frequency
whisk_freq              = ephys_data.conditions(1).whisk_frequency;
whisk_interval          = 1/whisk_freq; % interval between whisks

multi_whisk_wins        = NaN(n_stims,2);
for i = 1:n_stims
    multi_whisk_wins(i,:)   = whisk_resp_win + (i-1)*whisk_interval;
end


% Find the time window in the opto response in the test condition
opto_whisk_win_offset   = cond_data(1).whisk_onset - cond_data(1).LED_onset; % First condition should be opto and whisk
opto_resp_wins         	= cond_data(2).LED_onset + opto_whisk_win_offset + multi_whisk_wins;


for a = 1:length(cond_data)
    
    % Get spike times relative to whisk
    spikes                  = cond_data(a).spikes - cond_data(a).whisk_onset;
    
    % Spontaneous window relative to start of trial
    spont_resp_wins        	= multi_whisk_wins - cond_data(a).whisk_onset;
    
    % get spike rate / density data for all channels but don't make the plot
    [~, spike_density_rates(:,:,a)]                     = spike_density_plot(spikes,1,psth_bins, false);
    
    for b = 1:n_stims
        % Quantify spikes:
        
        % Whisker-evoked spike rate in target response window
        [whisk_spike_rates(:,a,b), whisk_spike_stds(:,a,b)]     = spike_rate_by_channel(spikes, multi_whisk_wins(b,:));
        
        % Spontaneous spikes (in window of equal size to whisk_resp_win but counting from start of trial)
        [spont_spike_rates(:,a,b), spont_spike_stds(:,a,b)]     = spike_rate_by_channel(spikes, spont_resp_wins);
        
        % Peak spike rate in response to whisker stimulus
        [peak_spike_rates(:,a,b), peak_spike_times(:,a,b)]      = peak_ROF_by_channel(spikes, multi_whisk_wins(b,:), rate_kernel_win);
        
        % First spike time in response to the whisker stimulus for each channel
        [first_spike_times(:,a,b), first_spike_jitter(:,a,b)]  	= first_spike_by_channel(spikes, multi_whisk_wins(b,:));
        
        % Spike probability by channel
        [spike_probs(:,a,b), n_hits(:,a,b), n_trials(a,b)]    	= spike_prob_by_channel(spikes, multi_whisk_wins(b,:));
        
       
        
        % Condition 2 = control condition (separate whisker and opto stimulus)
        if a == 2
            % Spike rate in opto-only stimulus window that corresponds to the time of
            % the whisker stimulus in the test condition
            [opto_spike_rates(:,1,b), opto_spike_stds]             = spike_rate_by_channel(spikes, opto_resp_wins(b,:));
            
            % Spike probability in opto-only stimulus window that corresponds to the time of
            % the whisker stimulus in the test condition
            [opto_spike_probs(:,1,b), opto_n_hits, ~]              = spike_prob_by_channel(spikes, opto_resp_wins(b,:));
        end
    end
    
end


%% Make some corrections to response measures

spont_spike_rates               = mean(mean(spont_spike_rates,2),3);
spont_spike_stds                = mean(mean(spont_spike_stds,2),3);

opto_spike_rates                = mean(opto_spike_rates,3);

% correct both for spontaneous
opto_spike_rates                = opto_spike_rates - spont_spike_rates;

whisk_spike_rates               = whisk_spike_rates - spont_spike_rates; % subtract mean of spontaneous; take mean of spontaneous over both conditions for more robust measure

% Correct for the expected number of spikes given the opto stimulation
whisk_spike_rates(:,1,:)        = whisk_spike_rates(:,1,:) - opto_spike_rates; % subtract opto component

peak_spike_rates                = peak_spike_rates - spont_spike_rates;
peak_spike_rates(:,1,:)      	= peak_spike_rates(:,1,:) - opto_spike_rates;
delta_peak_spike_rates         	= peak_spike_rates(:,1,:) - peak_spike_rates(:,2,:);

% Time difference measures
diff_peak_spike_times           = peak_spike_times(:,1,:) - peak_spike_times(:,2,:);
diff_first_spike_times          = first_spike_times(:,1,:) - first_spike_times(:,2,:);

control_whisk_spike_rates       = whisk_spike_rates(:,2,:);
opto_whisk_spike_rates          = whisk_spike_rates(:,1,:);

% Get difference in spike rate between +opto condition and control for a 
% one-data-point summary of the result of the experiment
delta_spike_rates               = opto_whisk_spike_rates - control_whisk_spike_rates;


%% Create output variable
drive_multi_data.opto_power             	= ephys_data.parameters.LED_power; % Report the power level

drive_multi_data.whisk_resp_win          	= whisk_resp_win;

drive_multi_data.control_whisk_spike_rates	= squeeze(control_whisk_spike_rates);
drive_multi_data.opto_whisk_spike_rates  	= squeeze(opto_whisk_spike_rates);
drive_multi_data.delta_spike_rates      	= squeeze(delta_spike_rates);

drive_multi_data.whisk_spike_rates      	= whisk_spike_rates;
drive_multi_data.whisk_spike_stds       	= whisk_spike_stds;

drive_multi_data.spont_spike_rates       	= spont_spike_rates;
drive_multi_data.spont_spike_stds       	= spont_spike_stds;

drive_multi_data.opto_spike_rates         	= opto_spike_rates;
drive_multi_data.opto_spike_stds            = opto_spike_stds;

drive_multi_data.peak_spike_rates         	= peak_spike_rates;
drive_multi_data.delta_peak_spike_rates  	= squeeze(delta_peak_spike_rates);

drive_multi_data.peak_spike_times         	= peak_spike_times;
drive_multi_data.diff_peak_spike_times   	= squeeze(diff_peak_spike_times);

drive_multi_data.first_spike_times       	= first_spike_times;
drive_multi_data.diff_first_spike_times 	= diff_first_spike_times;

drive_multi_data.whisk_spike_probabilities	= spike_probs;
drive_multi_data.whisk_n_hits            	= n_hits;
drive_multi_data.n_trials                	= n_trials;

drive_multi_data.opto_spike_probabilities 	= opto_spike_probs;
drive_multi_data.opto_n_hits             	= opto_n_hits;

drive_multi_data.psth_bins               	= psth_bins;
drive_multi_data.spike_density_rates     	= spike_density_rates;

