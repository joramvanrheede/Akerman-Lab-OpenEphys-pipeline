function [max_LFP_sink, min_profile, LFP_profile] = get_max_sink(LFP_traces,LFP_timestamps,resp_win,chan_lims,smooth_win)
% function  [max_LFP_sink, min_profile, LFP_profile] = get_max_sink(LFP_traces,LFP_timestamps,resp_win,chan_lims,smooth_win)
% 
% Finds the channel with the maximum LFP sink (minimum LFP value).
% 
% INPUTS:
% 
% LFP_TRACES: an n_channels * n_trials * n_samples matrix of LFP data
% 
% LFP_TIMESTAMPS: an 1:n_samples vector of timestamps associated with LFP_data.
% Defaults to [1:size(LFP_data,3)]/1000.
% 
% RESP_WIN: 2-element response window for determining the sink. Defaults to 
% [min(LFP_TIMESTAMPS) max(LFP_TIMESTAMPS)].
%
% CHAN_LIMS: Two-element vector that sets an upper and lower limit (minimum 
% and maximum channel number) on the channels where to look for a sink. E.g. 
% if CHAN_LIMS == [5 16] it will only look for a max sink between these channels.
% Default is all channels.
% 
% SMOOTH_WIN: moving smoothing window across channels to smooth out minor
% fluctuations between channels.
% 
% OUTPUT:
% 
% MAX_SINK_CHAN: Location (index) of the channel with the maximum sink.
% 
% MIN_PROFILE: Profile of the minimum values across channels.
% 
% LFP_PROFILE: an n_channels x n_samples matrix representing the mean LFP for 
% all channels, with the smoothing applied as used for deterining MAX_SINK_CHAN.
% 

if nargin < 2 || isempty(LFP_timestamps)
    LFP_timestamps = [1:size(LFP_traces,3)]/1000;
end

if nargin < 3 || isempty(resp_win)
    resp_win = [min(LFP_timestamps) max(LFP_timestamps)];
end

if nargin < 4 || isempty(chan_lims)
    chan_lims   = [1 size(LFP_traces,1)];
end

if nargin < 5
    smooth_win = 3;
end

% take mean across trials
LFP_means           = squeeze(mean(LFP_traces,2)); 

% use boolean to select only LFP data within resp_win
q_resp_win          = LFP_timestamps >= resp_win(1) & LFP_timestamps < resp_win(2); 
LFP_snippet         = LFP_means(:,q_resp_win); 

% find minimum for each channel
min_profile         = min(LFP_snippet,[],2);

% smooth across channels
min_profile         = smooth(min_profile,smooth_win); 

% find index of channel with minimum value, within chan_lims
max_LFP_sink        = find(min_profile == min(min_profile(chan_lims(1):chan_lims(2)))); 

% Generate smoothed LFP profile
LFP_profile         = smoothdata(LFP_means,1,'movmean',smooth_win); 
