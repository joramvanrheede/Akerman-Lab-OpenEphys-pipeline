function filtsignal = notch_filt(signal, sample_freq, filt_freq)
% filtsignal = notch_filt(signal, sample_freq, filt_freq)
% Does notch filter of signal at specified frequency (filt_freq), currently with fixed Q-factor
% of 35. Needs sample frequency of signal as input (sample_freq)
% NOTE: OPERATES ALONG THE FIRST NON-SINGLETON DIMENSION (see 'FILTER' function) 

wo                  = filt_freq/(sample_freq/2);       	% 50Hz
bw                  = wo/35;                           	% notch filter with Q-factor of 35 (= quite sharp)
[notch_B,notch_A]   = iirnotch(wo,bw);                 	% Create filter parameters

filtsignal        	= filter(notch_B,notch_A,signal); 	% Apply filter to LFP trace
