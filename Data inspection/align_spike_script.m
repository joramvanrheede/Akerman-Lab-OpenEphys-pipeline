% Script to quickly visualise the spike events from an in vivo multiunit
% ephys trace (based on simple thresholding)

% which file to use
% datafile                    = '/Volumes/PS2Akermanlab/MJB/In_Vivo/2019_04_04/CBLK_2019-04-04_18-10-49_10/100_CH16.continuous';
datafile                    = '/Volumes/PS2Akermanlab/MJB/In_Vivo/2019_09_16/CBLK_2019-04-04_18-10-49_10/100_CH31.continuous';
% datafile                    = '/Volumes/PS2Akermanlab/MJB/In_Vivo/2019_04_04/CBLK_2019-04-04_16-59-57_1/100_CH31.continuous' % No opto

threshold                   = 8;            % threshold in SDs for detecting spikes
range                       = 30;           % samples are plotted from -range to +range (centered on sample '0', the negative peak of the signal)
qplot                       = true;         % plot traces?
num_traces                  = 3000;        % how many traces to plot (number or 'all')

qreload                     = 1;

if qreload
    % Load data
    [data, timestamps, info] 	= load_open_ephys_data(datafile);
end

% Run align_spikes function which finds spikes and aligns them to the peak;
% It plots overlaid spike traces if requested.
aligned_traces              = align_spikes(data,threshold, range, qplot, num_traces);

hold on
