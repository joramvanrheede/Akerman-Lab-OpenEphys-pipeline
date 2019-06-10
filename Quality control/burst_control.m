function q_burst = burst_control(spikes, burst_win, burst_thresh, burst_chan_thresh, resp_onset)
% function q_burst = burst_control(spikes, burst_win, burst_thresh, resp_onset)
% 
% spikes: n_channels x n_trials x n_spikes array of spike times
% burst_win: time win relative to resp_onset in which to check for spikes
% burst_thresh: how many spikes (single trial, single channel) to count as
% a burst (defaults to 1)
% burst_chan_thresh: number of bursty channels on a single trial that make 
% that trial 'bursty' (defaults to 1)
% resp_onset: response onset (defaults to time 0)
% 
% q_burst: a vector of length n_trials which is TRUE for bursts and FALSE
% for non-bursts

% Default response onset = 0 (assumes 'spikes' has already been zeroed on response onset)
if nargin < 5
    resp_onset = 0;
end

% Default: a single bursty channel makes a bursty trial
if nargin < 4
    burst_chan_thresh = 1;
end

% Default: a single spike disqualifies the response for a channel
if nargin < 3
    burst_thresh = 1;
end

% set actual response onset to 0
spikes                      = spikes - resp_onset;

% Create matrix that is 1 for spike times within the window for detecting bursts & 0 otherwise
spikes_in_burst_win         = spikes > burst_win(1) & spikes <= burst_win(2);

% Sum over n_spikes dimension, should give a nr of spikes within window per channel and per trial
burst_spikes                = sum(spikes_in_burst_win,3);

% PREVIOUS:
% % Now sum over channels 
% burst_spikes_per_trial      = sum(burst_spikes,1);
% 
% % Use threshold to decide on bursts
% q_burst                     = burst_spikes_per_trial >= burst_thresh;

% Determine for each channel and each trial whether there is burstiness going on
is_bursty                   = burst_spikes >= burst_thresh;

% How many bursty channels are there for each trial?
n_bursty_chans          	= sum(is_bursty,1);

% use burst_chan_thresh to decide on bursts
q_burst                     = n_bursty_chans >= burst_chan_thresh;


