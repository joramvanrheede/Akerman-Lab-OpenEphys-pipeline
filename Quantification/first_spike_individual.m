function [first_spikes_by_trial_by_channel] = first_spike_individual(spikes, time_win)
% function [first_spikes_by_trial_by_channel] = first_spike_individual(SPIKES, TIME_WIN)
% 
% Returns FIRST_SPIKE_TIMES_BY_TRIAL_BY_CHANNEL in seconds in a given TIME_WIN. First spike time 
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
