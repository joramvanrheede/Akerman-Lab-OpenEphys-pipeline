% Script to quickly visualise the spike events from an in vivo multiunit
% ephys trace (based on simple thresholding)

% which file to use
datafile                    = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data/2018_04_19/ChR2-YFP-ON_2018-04-19_16-52-41_5/100_CH19.continuous';

threshold                   = 8;            % threshold in SDs for detecting spikes
range                       = 30;           % samples are plotted from -range to +range (centered on sample '0', the negative peak of the signal)
qplot                       = true;         % plot traces?
num_traces                  = 'all';        % how many traces to plot (number or 'all')

% Load data
[data, timestamps, info] 	= load_open_ephys_data(datafile);

% Run align_spikes function which finds spikes and aligns them to the peak;
% It plots overlaid spike traces if requested.
aligned_traces              = align_spikes(data,threshold, range, qplot, 'all');


