function [band_power, mean_band_power] = LFP_band_power(LFP_traces,freq_range,signal_freq,time_win)
% [BAND_POWER, MEAN_BAND_POWER] = LFP_band_power(LFP_TRACES,FREQ_RANGE,SIGNAL_FREQ,TIME_WIN)
% 
% Provides power spectra of LFP_TRACES for visualisation purposes.
% 
% BAND_POWER is an M*N matrix of band power values, where M and N are 
% determined by the first 2 dimensions of LFP_TRACES.
% Power spectra have been transformed by taking log(power+1) and then 
% dividing by the standard deviation to aid visualisation.
% 
% MEAN_BAND_POWER is a single value - the mean across M and N (e.g. across 
% channels and across trials), for quick access to the overall mean.
% 
% 
% INPUTS
% 
% LFP_TRACES: an M*N*N_SAMPLES matrix of LFP data (in current Akerman lab 
% data format, M = N_CHANNELS and N = N_TRIALS).
%
% FREQ_RANGE: the frequency band for which to calculate spectral power. Can
% be specified as [min_freq max_freq] or as a single value; if a single
% value is specified the range will be [freq+10%  freq-10%].
% 
% SIGNAL_FREQ: the sample rate of the signal in Hz (currently LFPs are 
% stored at 1000 Hz in Akerman lab data format).
% 
% TIME_WIN: Specifies the time window in seconds in which to calculate band 
% power; if not specified band power will be calculated across all samples.
% 
% Joram van Rheede, Akerman Lab, May 2019

% Get dimensions of input LFP traces
LFP_dims            = size(LFP_traces);

if nargin < 4
    % Use all samples of the signal
    sample_idx  	= 1:LFP_dims(3);
else
    % Convert time window to sample indices
    sample_idx      = (time_win(1)*signal_freq):(time_win(2)*signal_freq);
    sample_idx      = sample_idx + 1; % Avoid starting index at 
    if sample_idx(end) > LFP_dims(3)
        sample_idx(sample_idx > LFP_dims(3)) = [];
    end
end

% If only a target frequency is given, create a 10% margin frequency band around it
if isscalar(freq_range)
    freq_range  = [freq_range*0.9 freq_range*1.1];
end

band_power          = NaN(LFP_dims(1),LFP_dims(2));
for a = 1:LFP_dims(1) % loop over channels
    for b = 1:LFP_dims(2)
        
        this_LFP    = squeeze(LFP_traces(a,b,:));
        
        band_power(a,b)          = bandpower(this_LFP(sample_idx),signal_freq,freq_range);
        
    end
end

mean_band_power     = mean(band_power,'all');

