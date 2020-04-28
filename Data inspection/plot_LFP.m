% plot_LFP(chan_data)

% where are the data files?
data_folder     	= '/Volumes/Akermanlab/Joram/In_vivo_mouse_data/2019_01_17/CBLK_2019-01-17_11-26-02_6';       % currently set to 'cd' so code should work if script is run from directory containing the data

chan_nr             = 1;        % which channel? 1 = most superficial channel, then incrementing for deeper channels up to e.g. 16 or 32.

% Data resampling
resample_freq    	= 1000;     % Resample trace to this frequency (in Hz)

% Data filtering
LFP_band_max        = 300;      % LFP bandpass filter max
LFP_band_min        = .5;       % LFP bandpass filter min

do_50Hz_filt        = true;     % if true, notch filter for 50Hz noise

% Time within recording to plot
target_time         = 4500;     % Target time (in seconds) - sets the CENTRE of the time window displayed and zoomed in on
zoom_out_time       = 1000;     % Total amount of time to show in zoom out view (seconds)

% Zoom parameters
mid_zoom_factor     = 5;        % zoom factor relative to zoomed out trace
full_zoom_factor    = 5;        % zoom factor relative to mid zoom trace

% Which recording probe type?
probe_type          = '32ch'; 	% '16ch' for 16 channel or '32ch' for 32 channel

q_reload            = true;        % Re-load data or just re-plot with new parameters? Script will fail if 

%% Parameters that probably won't change much

sample_freq     	= 30000; 	% Sampling frequency of the data set

%% Currently supporting channel maps for 16 channel and 32 channel neuronexus linear probes

switch probe_type
    case '16ch'
        % 16ch linear silicon A16 probe neuronexus; corrected February 2019
        chan_map        = [13 29 4 20 15 31 3 19 16 32 1 17 2 18 14 30];
    case '32ch'
        % 32Ch linear A32 probe neuronexus; corrected February 2019
        chan_map        = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];
end

%% Load data
if q_reload
    target_channel                  = chan_map(chan_nr);    % Use the channel map to find appropriate channel number
    
    [chan_data timestamps info]     = load_open_ephys_data([data_folder filesep '100_CH' num2str(target_channel) '.continuous']);
elseif ~exist('chan_data')
    error('Data not loaded yet - set q_reload to true')
end

%% Apply bandpass filter to reveal LFP frequencies
[filt_B,filt_A]     = butter(2,[LFP_band_min LFP_band_max]/sample_freq/2);	% create 2nd order butterworth bandpass filter for LFP
LFP_trace        	= filter(filt_B,filt_A,chan_data);          % apply filter for LFP trace

%% Apply notch filter to remove 50Hz noise

if do_50Hz_filt
    wo                  = 50/(sample_freq/2);                 	% 50Hz
    bw                  = wo/35;                                % notch filter with Q-factor of 35 (= quite sharp)
    [notch_B,notch_A]   = iirnotch(wo,bw);                      % Create filter parameters
    
    LFP_trace        	= filter(notch_B,notch_A,LFP_trace);    % Apply filter to LFP trace
end

%% resampling (e.g. 1kHz more than enough for signal filtered for 300Hz and below)
LFP_trace_resampled = LFP_trace(1:sample_freq/resample_freq:end);                  % Resample to 1000Hz

time_vect           = [1:length(LFP_trace_resampled)]/resample_freq; % Create a vector of time stamps for the LFP data based on resampled 1kHz frequency


%% Work out what part of the LFP trace to use for different zoom level plots

zoom_mid_time       = zoom_out_time / mid_zoom_factor;      % Total amount of time to show in zoom mid view
zoom_in_time        = zoom_mid_time / full_zoom_factor;  	% Total amount of time to show in zoom in view

% Generate vector of samples for the zoomed out trace
zoom_out_t_start 	= target_time - (zoom_out_time / 2);
zoom_out_t_end   	= target_time + (zoom_out_time / 2);
zoom_out_samples    = (zoom_out_t_start * resample_freq):(zoom_out_t_end * resample_freq);

% Generate vector of samples for the mid-zoom trace
zoom_mid_t_start 	= target_time - (zoom_mid_time / 2);
zoom_mid_t_end   	= target_time + (zoom_mid_time / 2);
zoom_mid_samples    = (zoom_mid_t_start * resample_freq):(zoom_mid_t_end * resample_freq);

% Generate vector of samples for the zoomed in trace
zoom_in_t_start 	= target_time - (zoom_in_time / 2);
zoom_in_t_end   	= target_time + (zoom_in_time / 2);
zoom_in_samples     = (zoom_in_t_start * resample_freq):(zoom_in_t_end * resample_freq);


%% Plotting at 3 zoom levels

subplot(3,1,1)
plot(time_vect(zoom_out_samples),LFP_trace_resampled(zoom_out_samples),'k-')
xlabel('Time (s)')
axis tight
title('In vivo dysrhythmia following local 4AP injection','FontSize',24)

subplot(3,1,2)
plot(time_vect(zoom_mid_samples),LFP_trace_resampled(zoom_mid_samples),'k-')
xlabel('Time (s)')
axis tight

subplot(3,1,3)
plot(time_vect(zoom_in_samples),LFP_trace_resampled(zoom_in_samples),'k-')
xlabel('Time (s)')
axis tight

% Enlarge figure
set(gcf,'Units','normal')
set(gcf,'Position',[0 .1 1 .8])

