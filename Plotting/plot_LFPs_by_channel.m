function plot_handle = plot_LFPs_by_channel(LFP_traces,LFP_timestamps,spacing, alpha_level)
% function PLOT_HANDLE = plot_LFPs_by_channel(LFP_TRACES,LFP_TIMESTAMPS,SPACING,ALPHA_LEVEL)
% LFP plotting function that shows the LFP response by channel, with an
% offset between channels determined by 'spacing'
% 
% LFP_TRACES: An N_channels by N_samples matrix of LFP data.
% 
% LFP_TIMESTAMPS: The time stamps for the values in LFP_traces
%
% SPACING: Specifies how far traces are separated, default value is 0.1.
% 
% ALPHA_LEVEL: specifies the alpha level (transparency) of the traces,
% default is 0.3.


% Set default 'time stamps' to x = 1:n_samples if not specified by user
if nargin < 2
    LFP_timestamps = 1:length(LFP_traces);
end

% Set default spacing to 1
if nargin < 3
    spacing = 0.1;
end

if nargin < 4
    alpha_level = .3;
end

LFP_traces          = LFP_traces'; % LFP_traces needs to be in a different orientation for the functions to follow

% Determine the 'median range' of the traces which is used for determining the spacing between traces
median_LFP_min      = median(min(LFP_traces)); % for each trace, determine the max, and then the median of all trace maxima
median_LFP_max      = median(max(LFP_traces)); % for each trace, determine the min, and then the median of all trace minima
LFP_range           = median_LFP_max - median_LFP_min; 

% Normalise traces by the range so they should now have a median extent of ~1
LFP_traces          = LFP_traces / LFP_range;

% Make a matrix of offsets of the same size as LFP_traces, but that decrements by 1 for every next trace:
LFP_offsets         = meshgrid((size(LFP_traces,2)+1)-(1:size(LFP_traces,2)),1:size(LFP_traces,1));

% Scale the traces for spacing, then add the offsets
LFP_offset_traces   = LFP_traces/spacing + LFP_offsets;

% Plot the traces
plot_handle = plot(LFP_timestamps,LFP_offset_traces,'LineWidth',2,'Color',[0 0 0 alpha_level]);

% Some plot aesthetics
axis tight
axis ij
ylimz   = ylim;
ylim([ylimz(1)-2*spacing ylimz(2)+2*spacing])

% Plot labeling
ylabel('Trace number')
xlabel('Time')
set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2)

% To do - create scale bar?


