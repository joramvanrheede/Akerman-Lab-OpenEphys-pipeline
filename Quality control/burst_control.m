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

% Default response onset = 0 (assumes 'spikes' has already been zeroed on response onset)
if nargin < 4
    resp_onset = 0;
end

% Default is a single spike disqualifies the response
if nargin < 3
    burst_thresh = 1;
end

% set actual response onset to 0
spikes                      = spikes - resp_onset;

% Create matrix that is 1 for spike times within the window for detecting bursts & 0 otherwise
spikes_in_burst_win         = spikes > burst_win(1) & spikes <= burst_win(2);

% Sum over n_spikes dimension, should give a nr of spikes within window per channel and per trial
burst_spikes                = sum(spikes_in_burst_win,3);

% Now sum over channels 
burst_spikes_per_trial      = sum(burst_spikes,1);

% Use threshold to decide on bursts
q_burst                     = burst_spikes_per_trial >= burst_thresh;

