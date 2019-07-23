function spike_rate = spike_rate_in_win(spikes, time_win)
% function spike_rate = spike_rate_in_win(SPIKES, TIME_WIN)
% 
% Returns SPIKE_RATE in Hz in a given TIME_WIN. Spike rate is average 
% for all channels and all trials in SPIKES
% 
% SPIKES: a N_CHANNELS * N_TRIALS * N_SPIKES matrix of spike times, 
% padded with NaNs for empty values. Function assumes spike times are in
% seconds.
% 
% TIME_WIN: a time window, [T1 T2]. Function counts spikes from 
% TSPIKE >= T1 to TSPIKE <= T2

% Count spike times within time window
spikes_in_win       = spikes >= time_win(1) & spikes <= time_win(2);
spike_count         = sum(spikes_in_win(:));

% Work out mean spike count over trials and channels
n_channels          = size(spikes,1);
n_trials            = size(spikes,2);
mean_spike_count  	= spike_count / n_channels / n_trials;

% How long is time window?
delta_t             = time_win(2) - time_win(1);

% Convert to rate by dividing by delta_t
spike_rate          = mean_spike_count / delta_t;