function [spike_rate, spike_rate_std] = spike_rate_in_win(spikes, time_win)
% function [spike_rate, spike_rate_std] = spike_rate_in_win(SPIKES, TIME_WIN)
% 
% Returns SPIKE_RATE in Hz in a given TIME_WIN. Spike rate is average 
% for all channels and all trials in SPIKES. Also returns SPIKE_RATE_STD, 
% the standard deviation over trials.
% 
% SPIKES: a N_CHANNELS * N_TRIALS * N_SPIKES matrix of spike times, 
% padded with NaNs for empty values. Function assumes spike times are in
% seconds.
% 
% TIME_WIN: a time window, [T1 T2]. Function counts spikes from 
% TSPIKE >= T1 to TSPIKE <= T2

% Count spike times within time window
spikes_in_win                       = spikes >= time_win(1) & spikes <= time_win(2);
spike_counts_by_trial_by_channel    = sum(spikes_in_win,3);
spike_counts_by_trial               = mean(spike_counts_by_trial_by_channel,1)';
mean_spike_count                    = mean(spike_counts_by_trial);
std_spike_count                     = std(spike_counts_by_trial);

% How long is time window?
delta_t             = time_win(2) - time_win(1);

% Convert to rate by dividing by delta_t
spike_rate          = mean_spike_count / delta_t;
spike_rate_std      = std_spike_count / delta_t;
