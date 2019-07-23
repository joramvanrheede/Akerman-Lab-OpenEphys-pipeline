function spike_count = spike_count_in_win(spikes, time_win)
% function SPIKE_COUNT = spike_count_in_win(SPIKES, TIME_WIN)
% 
% Returns SPIKE_COUNT in in a given TIME_WIN. Spike count is average 
% over all channels and all trials in SPIKES
% 
% SPIKES: a N_CHANNELS * N_TRIALS * N_SPIKES matrix of spike times, 
% padded with NaNs for empty values. Function assumes spike times are in
% seconds.
% 
% TIME_WIN: a time window, [T1 T2]. Function counts spikes from 
% TSPIKE >= T1 to TSPIKE <= T2

% Count spike times within time window
spikes_in_win       = spikes >= time_win(1) & spikes <= time_win(2);
all_spike_count  	= sum(spikes_in_win(:));

% Work out mean spike count per trial per channel
n_channels          = size(spikes,1);
n_trials            = size(spikes,2);
spike_count         = all_spike_count / n_channels / n_trials;
