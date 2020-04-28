function [spike_probs, n_hits, n_trials]  = spike_prob_by_channel(spikes, time_win)
% function [spike_probs, n_hits, n_trials]  = spike_prob_by_channel(SPIKES, TIME_WIN)
% 
% Returns SPIKE_PROBS, estimated spike probabilities as a number between 
% 0 and 1 in a given TIME_WIN. SPIKE_PROBS is calculated for each channel / 
% unit in SPIKES. Also returns N_HITS, the number of trials with one or more 
% spikes per channel, and N_TRIALS, the overall number of trials. 
% The latter two should be sufficient to do a chi square or Fisher's exact 
% test, or similar, on spike fractions.
% 
% SPIKES: a N_CHANNELS * N_TRIALS * N_SPIKES matrix of spike times, 
% padded with NaNs for empty values. Function assumes spike times are in
% seconds. (N_CHANNELS can also be N_UNITS)
% 
% TIME_WIN: a time window, [T1 T2]. Function counts spikes from 
% TSPIKE >= T1 to TSPIKE <= T2

% Count spike times within time window
spikes_in_win                       = spikes >= time_win(1) & spikes <= time_win(2);

% Find 'hit' trials (with at least 1 spike)
if size(spikes_in_win,3) > 1 % if not working with squeezed data eg n units x m trials x c spikes
    spike_hits_by_trial_by_channel      = sum(spikes_in_win,3) > 0;
    % Take sum over trials to get number of 'hit' trials
    n_hits                              = sum(spike_hits_by_trial_by_channel,2);
    % Take mean over trials to get estimated spike probability
    spike_probs                         = mean(spike_hits_by_trial_by_channel,2);
    % Report number of trials as well, for ease of following up with chi square
    % or Fisher's exact test
    n_trials                            = size(spikes, 2);

else % if working with squeezed data (m trials x c spikes).
    spike_hits_by_trial_by_channel      = sum(spikes_in_win,2) > 0;
    n_hits                              = sum(spike_hits_by_trial_by_channel,1);
    spike_probs                         = mean(spike_hits_by_trial_by_channel,1);
    n_trials                            = size(spikes, 1);
end; % if statement 
end

