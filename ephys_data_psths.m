function [ephys_data] = ephys_data_psths(ephys_data, hist_binsize, rate_window)
% function [ephys_data] = ephys_data_psths(ephys_data, hist_binsize, rate_window)
% This function generates PSTHs, mean LFP traces and a spike rate trace,
% and puts them in the ephys_data struct.


hist_min            = 0;                        % histogram min time in seconds
hist_max            = ephys_data.trial_length;  % histogram max time in seconds

samplerate          = 30000;                    % sample rate

%% 
hist_edges          = [hist_min:hist_binsize:hist_max];

for a = 1:length(ephys_data.conditions)

    these_spikes        = ephys_data.conditions(a).spikes;
    
    trial_psths         = histc(these_spikes,hist_edges,3);
    trial_rate_psths    = trial_psths / hist_binsize;
    expt_psths          = squeeze(mean(trial_psths,2));
    expt_rate_psths     = expt_psths / hist_binsize;
    
    these_LFPs          = ephys_data.conditions(a).LFP_trace;
    expt_LFPs           = squeeze(mean(these_LFPs,2));
    
    %% Using the spike times for all episodes, make a profile
    % of the firing rate, in order to determine the maximum rate
    % of fire and the time of that peak ROF
    
    % create a vector of bins at the samplerate and count spikes
    spike_rate_bins 	= hist_min:1/samplerate:hist_max;
    spikes_by_channel   = reshape(these_spikes,[size(these_spikes,1) size(these_spikes,2)*size(these_spikes,3)]); 
    spike_rate_counts 	= histc(spikes_by_channel,spike_rate_bins,2);
    
    % determine the number of samples for the smoothing kernel from the
    % smoothing window (s) and sample rate (Hz)
    n_smooth_samples            = rate_window * samplerate;
    
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
    n_episodes                  = size(ephys_data.conditions(1).spikes,2);
    spike_rate_profile          = spike_rate_profile * samplerate / n_episodes;
    
    %% put the psth / mean / spike rate data back in the ephys_data struct
    ephys_data.conditions(a).psth_bins          = hist_edges;
    ephys_data.conditions(a).trial_psths        = trial_psths;
    ephys_data.conditions(a).trial_rate_psths   = trial_rate_psths;
    ephys_data.conditions(a).expt_psths         = expt_psths;
    ephys_data.conditions(a).expt_rate_psths    = expt_rate_psths;
    ephys_data.conditions(a).expt_LFPs          = expt_LFPs;
    ephys_data.conditions(a).spike_rate_times   = spike_rate_bins;
    ephys_data.conditions(a).spike_rate_profile = spike_rate_profile;
end



