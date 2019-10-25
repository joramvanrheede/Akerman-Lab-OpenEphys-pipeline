function [psth_handle, counts, binedges] = psth(spike_times,bin_size,psth_win,show_plot)
% function [PSTH_HANDLE COUNTS BINEDGES] = psth(SPIKE_TIMES,BIN_SIZE,PSTH_WIN,SHOW_PLOT)
% OR:
% function [PSTH_HANDLE COUNTS BINEDGES] = psth(SPIKE_TIMES,BINEDGES,SHOW_PLOT)
% 
% Makes post-stimulus time histogram of (spike) time data in SPIKE_TIMES.
% Each data point in SPIKE_TIMES is a spike time (NaNs are ignored).
% SPIKE_TIMES is vectorised ( 'SPIKE TIMES = SPIKE_TIMES(:)' ) so the input
% data can have any shape or number of dimensions, but it will be treated
% as a 1D vector.
% 
% BIN_SIZE specifies the width of the bins in the PSTH in the same units as
% spike_times (e.g. seconds, milliseconds, etc). If not provided or left 
% empty this will be decided from the density of the data. IF BIN_SIZE is a
% vector it will be interpreted as user-defined BINEDGES instead -->
%
% BINEDGES is a vector that explicitly specifies the edges of the bins used 
% for the PSTH, and when specified, makes PSTH_WIN obsolete. 
% 
% PSTH_WIN is a 2-element vector [TMIN TMAX] that specifies the time range 
% over which the PSTH is plotted, e.g. if SPIKE_TIMES includes spike times 
% over a 10 second period, but PSTH_WIN is [0 1], only spike times over the
% first second are included in the PSTH. If left blank this will be
% decided by the min and the max of the data.
%
% SHOW_PLOT is a boolean variable (0/1 or true/false) that indicates whether
% to plot the count data as a bar histogram or not. Default is true; set to
% false / 0 if you only need the counts but don't want to plot the figure.
% 
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

%
if exist('bin_size','var') && length(bin_size) > 1
    % If the bin_size input instead contains BIN EDGES, set flag q_edges to
    % true and rename variables appropriately
    binedges        = bin_size;
    
    if nargin < 3
        show_plot   = true;
    else
        show_plot 	= psth_win; % Otherwise, 'psth_win' should now contain the 'show_plot' variable
    end
    
    bin_size        = median(diff(binedges));
    psth_win        = [binedges(1) binedges(end)];
    q_edges         = true;
    
else
    q_edges         = false;
    
    if nargin < 4
        show_plot   = true; % Default is to show the PSTH plot
    end
end

% Vectorise spike_times into an N spikes x 1 vector
spike_times     = spike_times(:);

% Let input variables determine how to run histcounts:
if q_edges
    [counts binedges] 	= histcounts(spike_times,binedges);
elseif ~exist('psth_win','var') && ~exist('bin_size','var')
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

% If requested, make plot
if show_plot
    % Plot the histogram using Matlab's bar function; offset binedges by .5
    % binsize so that the edge of the first bar starts at 0, set colour to
    % black, and set bar width to 1 so there are no gaps between bars
    
    psth_handle     = bar(binedges(1:end-1)+0.5*bin_size,counts,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'BarWidth',1); % plot the spike counts vs time bins
    % Tighten the x-axis around the bar plot
    xlim(psth_win)
    
    % Basic labels for axes
    xlabel('Time')
    ylabel('Spike count')
else
    psth_handle     = [];
end

