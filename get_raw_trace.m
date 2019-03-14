function raw_trial_trace = get_raw_trace(ephys_data_file, condition_nr, trial_nr, channel_nr, raw_data_folder)
% get_raw_trace(ephys_data_file, condition_nr, trial_nr, raw_data_folder, channel_number) 
% 
% get raw data trace for a given data file - trial, condition, trial number, etc.
% needs the name of an ephys_data_file to load, and needs the full path to
% the raw data folder where the raw data lives

% constants:
SAMPLEFREQ      = 30000;    % Sampling frequency of the raw data
SPK_MIN_FREQ    = 500;      % Low bound of the bandpass filter for spike trace
SPK_MAX_FREQ    = 6000;     % High bound of the bandpass filter for spike trace

% load preprocessed data file with sync info; this will load variable
% 'ephys_data'
load(ephys_data_file)

% get trial start data
trial_starts = ephys_data.conditions(condition_nr).trial_starts;

% get trial end data
trial_ends = ephys_data.conditions(condition_nr).trial_ends;

% get start and end time for the target trial
trace_start     = trial_starts(trial_nr);
trace_end       = trial_ends(trial_nr);

% get file name for this channel number
channel_file_name = dir(fullfile(raw_data_folder, sprintf('*CH%d.continuous', channel_nr)));

% load data
[channel_trace, timestamps, info]   = load_open_ephys_data([channel_file_name.folder filesep channel_file_name.name]);

% set all timing relative to the start of this data file
timestamps  = timestamps - min(timestamps(:));

% determine what part of the trace falls within the target trial window
trace_win   = timestamps > trace_start & timestamps < trace_end;

% extract the relevant part of the trace
raw_trial_trace     = channel_trace(trace_win);

% filter for spikes
samplefreq          = 30000;
[filt_b,filt_a]  	= butter(2, [500 6000]/(samplefreq/2));
raw_trial_trace    	= filter(filt_b,filt_a,raw_trial_trace); % filter data with butterworth bandpass filter to get spike traces
