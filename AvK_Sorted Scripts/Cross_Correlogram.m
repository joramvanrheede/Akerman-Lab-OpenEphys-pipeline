function [N,bin_edges] = Cross_Correlogram(ref_spikes,target_spikes,lag_window, bin_size);
%% Computes a Cross Correlogram for two spike trains
% Inputs :
%           NB ref_spikes and target spikes must be same size matrices
%            ref_spikes : n (trials) x m (no of spikes, NaN padded) array 
%                       of spike times in seconds
%           target_spikes : n (trials) x m (no of spikes, NaN padded) array
%                           of spike times in seconds;
%           lag window : time in ms before and after spike to test
%                        correlation
%           bin size : size in ms of each bin. 
% Outputs : 
%           N : array of -lagwindow : bin_size : + lag
%                               window size. with spike count in each bin. 
%           Edges : Bin_edges used
%   Alex von klemperer 2020

bin_edges = 0-lag_window : bin_size : lag_window; % creates bin edges of size bin_size
bin_edges = bin_edges/1000; % converts to time in seconds.

for j = 1 : size(ref_spikes,1) % for each trial
    ref = ref_spikes(j,:); % creates reference spike train for trial
    target = target_spikes(j,:); % creates target spike train for trial
        for k = 1:size(ref_spikes,2) %  for each spike in trial
            r = ref(k); % gets time of reference spike
            t = target-r; % gets timing of each spike in target train relative to refernce spike
            h(k,:) = histcounts(t,bin_edges); % calculates spike counts in each bin
        end;
    h = sum(h,1); % sums for all reference spikes
    i(j,:) = h; % creates temp count for each trial
end;
i = sum(i,1); % sums for all trials
N = i;
end
