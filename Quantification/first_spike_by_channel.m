function [first_spike_times, first_spike_jitters] = first_spike_by_channel(spikes, time_win)
% function [first_spike_times, first_spike_jitters] = first_spike_by_channel(SPIKES, TIME_WIN)
% 
% Returns FIRST_SPIKE_TIMES in seconds in a given TIME_WIN. First spike time 
% is average by channel over all trials in SPIKES. FIRST_SPIKE_JITTERS gives
% the standard deviation of the first spike time over all trials.
% 
% SPIKES: a N_CHANNELS * N_TRIALS * N_SPIKES matrix of spike times, 
% padded with NaNs for empty values. Function assumes spike times are in
% seconds.
% 
% TIME_WIN: a time window, [T1 T2]. Function counts spikes from 
% TSPIKE >= T1 to TSPIKE <= T2

% Count spike times within time window
spikes_in_win                       = spikes >= time_win(1) & spikes <= time_win(2); % boolean of which spike times fall within window

spike_times_in_win                  = spikes; % Copy original spike times
spike_times_in_win(~spikes_in_win)  = NaN; % Set spike times out of window to 0

% Find first spike in window by trial by channel
first_spikes_by_trial_by_channel    = nanmin(spike_times_in_win,[],3);

% Find the mean first spike time over all trials, by channel
first_spike_times                   = nanmean(first_spikes_by_trial_by_channel,2);

% Find the jitter (standard deviation) of first spike time over all trials, by channel
first_spike_jitters                 = nanstd(first_spikes_by_trial_by_channel,[],2);
