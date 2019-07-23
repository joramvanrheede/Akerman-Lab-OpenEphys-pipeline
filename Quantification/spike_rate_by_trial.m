function spike_rates = spike_rate_by_trial(spikes, time_win)
% function spike_rates = spike_rate_by_trial(SPIKES, TIME_WIN)
% 
% Returns SPIKE_RATES in Hz in a given TIME_WIN. Spike rate is average by 
% trial for all channels in SPIKES
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

% How long is time window?
delta_t                             = time_win(2) - time_win(1);

% Convert to rate by dividing by delta_t
spike_rates                         = spike_counts_by_trial / delta_t;