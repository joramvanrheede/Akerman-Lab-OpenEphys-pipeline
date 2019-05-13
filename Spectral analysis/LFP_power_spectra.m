function [power_spectra, mean_power_spectrum] = LFP_power_spectra(LFP_traces,freq_bins,LFP_signal_freq,spect_win_size,spect_win_overlap)
% [POWER_SPECTRA, MEAN_POWER_SPECTRUM] = LFP_POWER_SPECTRA(LFP_TRACES,FREQ_BINS,LFP_SIGNAL_FREQ,SPECT_WIN_SIZE,SPECT_WIN_OVERLAP)
% 
% Provides power spectra of LFP_TRACES for visualisation purposes.
% 
% POWER_SPECTRA is an M*N*N_COLUMNS*N_ROWS matrix of power spectra, where
% M and N are determined by the first 2 dimensions of LFP_TRACES, N_COLUMNS 
% is determined by N_SAMPLES of LFP_traces and SPECT_WIN_SIZE, and N_ROWS
% determined by FREQ_BINS.
% Power spectra have been transformed by taking log(power+1) and then 
% dividing by the standard deviation to aid visualisation.
% 
% MEAN_POWER_SPECTRUM is an N_COLUMNS*N_ROWS matrix with a mean across M
% and N (e.g. across channels and across trials), for quick visualisation
% of the overall spectrum, e.g.:
% 
% imagesc(MEAN_POWER_SPECTRUM), axis xy
% 
% INPUTS
% 
% LFP_TRACES: an M*N*N_SAMPLES matrix of LFP data (in current Akerman lab 
% data format, M = N_CHANNELS and N = N_TRIALS).
%
% FREQ_BINS: the frequency bins over which to compute the spectrogram.
% 
% LFP_SIGNAL_FREQ: the sample rate of the signal in Hz (currently LFPs are 
% stored at 1000 Hz in Akerman lab data format).
% 
% SPECT_WIN_SIZE: The windows (number of samples) over which to compute the 
% columns of the spectrogram.
% 
% SPECT_WIN_OVERLAP (OPTIONAL): Specifies an overlap window (in samples)
% from one column of the spectrogram to the next.
% 
% Joram van Rheede, Akerman Lab, May 2019

% Get dimensions of input LFP traces
LFP_dims            = size(LFP_traces);

% Set default to no overlap between windows
if nargin < 5
    spect_win_overlap = 0;
end

% 100 samples should allow for detecting a reasonable range of frequencies
if nargin < 4
    spect_win_size  = 100;
end

% 
[s, f, t, temp_power] = spectrogram(squeeze(LFP_traces(1,1,:)),spect_win_size,spect_win_overlap,freq_bins,LFP_signal_freq);
 
spect_dims          = size(temp_power);
 
power_spectra     	= NaN(LFP_dims(1),LFP_dims(2),spect_dims(1),spect_dims(2));
% target_band_power   = NaN(LFP_dims(1),LFP_dims(2));
for a = 1:LFP_dims(1) % loop over channels
    for b = 1:LFP_dims(2)
        this_LFP    = squeeze(LFP_traces(a,b,:));
        [spect, freqs, times, power]    = spectrogram(this_LFP,spect_win_size,spect_win_overlap,freq_bins,LFP_signal_freq);
        
        power = log(power+1);
        power = power / std(power(:));
        
        power_spectra(a,b,:,:)  = power;
        
%         freq_range                      = [(stim_freq*.9) (stim_freq*1.1)];
%         target_band_power(b,c)          = bandpower(this_LFP(power_win(1):power_win(2)),LFP_signal_freq,freq_range);
%         
    end
end

mean_power_spectrum     = squeeze(mean(mean(power_spectra,1),2));

