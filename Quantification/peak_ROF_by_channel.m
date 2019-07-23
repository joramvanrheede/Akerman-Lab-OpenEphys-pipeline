function [peak_rate_of_fire, peak_time] = peak_ROF_by_channel(spikes,time_win,kernel_size,sample_rate)
% function [PEAK_RATE_OF_FIRE, PEAK_TIME] = peak_ROF_by_channel(SPIKES,TIME_WIN,KERNEL_SIZE, SAMPLE_RATE)
% 
% Returns estimate of peak rate of fire within a time window as well as the 
% time of this peak, by channel.
% Peak rate of fire is estimated by aggregating spike events using a gaussian 
% kernel spanning 3 standard deviations.
% 
% PEAK_RATE_OF_FIRE: estimated peak rate of fire in Hz, by channel
% 
% PEAK_TIME: time of this peak, in seconds, by channel
% 
% 
% SPIKES: an N_CHANNELS * N_TRIALS * N_SPIKES matrix of spike times, padded with
% NaNs for empty spike time slots
%
% TIME_WIN: a time window in seconds, [T1 T2], in which to look for a peak rate of fire
% 
% OPTIONAL ARGUMENTS:
%
% KERNEL_SIZE: the extent of the 3 SD gaussian aggregation window to 
% estimate peak spike rate in seconds (default = 0.02s)
%
% SAMPLE_RATE: the sample rate of the signal (default = 30000Hz)


% Set defaults
if nargin < 3
    kernel_size = 0.02; % 20ms seems reasonable default
end

if nargin < 4
    sample_rate = 30000; % standard acquisition rate of OpenEphys system
end

% set range for looking for peak:
hist_min = time_win(1);
hist_max = time_win(2);

% create a vector of bins at the samplerate and count spikes
spike_rate_bins 	= hist_min:1/sample_rate:hist_max;
spikes_by_channel   = reshape(spikes,[size(spikes,1) size(spikes,2)*size(spikes,3)]);
spike_rate_counts 	= histc(spikes_by_channel,spike_rate_bins,2);

% determine the number of samples for the smoothing kernel from the
% smoothing window (s) and sample rate (Hz)
n_smooth_samples            = kernel_size * sample_rate;

% for convolution, a kernel with an odd number of elements is
% needed; add extra element if it is even
if mod(n_smooth_samples,2) == 0
    n_smooth_samples        = n_smooth_samples + 1;
end

kernel_vect                 = [-3:(6/n_smooth_samples):3];  % let vector run from -3(SDs) to +3(SDs)
gauss_kernel                = normpdf(kernel_vect,0,1);     % create gaussian using the normal probability density function

% ensure that the kernel has a total area under the curve of 1; this
% keeps the units of spike rate essentially the same
gauss_kernel                = gauss_kernel / sum(gauss_kernel);

% use kernel for convolution with the spike counts binned at sample
% rate (using conv2 this can be done in parallel on all 32 channels)
spike_rate_profile       	= conv2(spike_rate_counts,gauss_kernel,'same');

% convert the spike rate profile into an estimate of actual spike rate
% by multiplying by the sample rate (== dividing by time window) and by
% dividing by the number of episodes
n_episodes                      = size(spikes,2);
spike_rate_profile              = spike_rate_profile * sample_rate / n_episodes;

[peak_rate_of_fire, maxinds]    = max(spike_rate_profile, [],2);
peak_time                       = spike_rate_bins(maxinds)';
