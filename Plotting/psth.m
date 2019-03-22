function [psth_handle counts binedges] = psth(spike_times,bin_size,psth_win)
% function [PSTH_HANDLE COUNTS BINEDGES] = psth(SPIKE_TIMES,BIN_SIZE,PSTH_WIN)
% 
% Makes post-stimulus time histogram of (spike) time data in SPIKE_TIMES.
% Each data point in SPIKE_TIMES is a spike time (NaNs are ignored).
% SPIKE_TIMES is vectorised ( 'SPIKE TIMES = SPIKE_TIMES(:)' ) so the input
% data can have any shape or number of dimensions, but it will be treated
% as a 1D vector.
% 
% BIN_SIZE specifies the width of the bins in the PSTH in the same units as
% spike_times (e.g. seconds, milliseconds, etc).
% 
% PSTH_WIN is a 2-element vector [TMIN TMAX] that specifies the time range 
% over which the PSTH is plotted, e.g. if SPIKE_TIMES includes spike times 
% over a 10 second period, but PSTH_WIN is [0 1], only spike times over the
% first second are included in the PSTH.
%
% PSTH_HANDLE is a handle to the bar graph object that is used to plot the
% PSTH - this means that it will have the same properties as any matlab bar
% graph created with the 'bar' function, which can be modified for
% aesthetics by the user after running this script.
% 
% COUNTS provides the spike counts that are used for the histogram.
%
% BINEDGES provides the edges of the bins that are used for the histogram.
% 
% Joram van Rheede 22/03/2019


% Vectorise spike_times into an N spikes x 1 vector
spike_times     = spike_times(:);

% Let input variables determine how to run histcounts:
if ~exist('psth_win','var') && ~exist('bin_size','var')
    % Use histcounts built in best bin size and bin limits set by the range of the data
	[counts binedges]	= histcounts(spike_times);
    psth_win            = [min(binedges) max(binedges)];
    bin_size            = binedges(2) - binedges(1);
elseif exist('bin_size','var') && ~exist('psth_win','var')
    % Use bin limits determined by the range of the data but use user specified bin size
	[counts binedges]	= histcounts(spike_times,'BinWidth',bin_size);
    psth_win            = [min(binedges) max(binedges)];
elseif exist('bin_size','var') && isempty(bin_size) && exist('psth_win','var')
    % Use histcounts built in best bin size but use user specified limits
	[counts binedges]	= histcounts(spike_times,'BinLimits',psth_win);
    psth_win            = [min(binedges) max(binedges)];
    bin_size            = binedges(2) - binedges(1);
else
    % Use user-specified bin size and psth_win
    binedges            = (psth_win(1):bin_size:psth_win(2));
    counts              = histcounts(spike_times,binedges);
end

% Plot the histogram using Matlab's bar function; offset binedges by .5
% binsize so that the edge of the first bar starts at 0, set colour to
% black, and set bar width to 1 so there are no gaps between bars
psth_handle     = bar(binedges(1:end-1)+0.5*bin_size,counts,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',1); % plot the spike counts vs time bins

% Tighten the x-axis around the bar plot
xlim(psth_win)

% Basic labels for axes
xlabel('Time')
ylabel('Spike count')
