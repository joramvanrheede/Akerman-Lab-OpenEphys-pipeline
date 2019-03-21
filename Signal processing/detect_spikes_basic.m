function spike_times = detect_spikes_basic(spike_signal,n_SDs,sample_rate)
% function spike_times = detect_spikes_basic(SPIKE_SIGNAL,N_SDS,SAMPLE_RATE)
% 
% Basic spike detection function.
% 
% Takes filtered SPIKE_SIGNAL (vector) and returns spike times based on
% crossings of threshold set by N_SDS (number of standard deviations).
% 
% N_SDS can be negative, in which case 
% 
% SAMPLE_RATE (Hz) is used to convert spike position in the vector into 
% spike times in seconds.
% 
% Joram van Rheede 21/03/2019

SD              	= std(spike_signal); % get standard deviation

% Detect events based on positive or negative peaks?
if n_SDs > 0
    above_threshold     = spike_signal > n_SDs * SD; % determine when signal is above threshold
else
    above_threshold  	= spike_signal < n_SDs * SD; % determine when signal is above threshold
end

threshold_crossings = diff(above_threshold); % will return +1 for threshold crossing samples only

spike_samples       = find(threshold_crossings == 1); % only spike onsets

spike_times         = spike_samples / sample_rate; % convert sample indices to spike times