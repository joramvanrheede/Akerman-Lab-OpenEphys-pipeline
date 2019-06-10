function q_burst = burst_control(spikes, burst_win, burst_thresh, resp_onset)
% function q_burst = burst_control(spikes, burst_win, burst_thresh, resp_onset)
% 
% spikes: n_channels x n_trials x n_spikes array of spike times
% burst_win: time win relative to resp_onset in which to check for spikes
% burst_thresh: how many spikes (sum over all channels) to count as a burst
% (defaults to 1)
% resp_onset: response onset (defaults to time 0)
% 
% q_burst: a vector of length n_trials which is TRUE for bursts and FALSE
% for non-bursts

if nargin < 4
    resp_onset = 0;
end
if nargin < 3
    burst_thresh = 1;
end

spikes_in_burst_win         = spikes > burst_win(1) & spikes <= burst_win(2);

burst_spikes                = sum(spikes_in_burst_win,3);

burst_spikes_per_trial      = sum(burst_spikes,1);

q_burst                     = burst_spikes_per_trial >= burst_thresh;

