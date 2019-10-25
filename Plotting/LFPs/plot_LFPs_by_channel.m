function plot_handle = plot_LFPs_by_channel(LFP_traces,split_dim,LFP_timestamps,spacing,alpha_level)
% function PLOT_HANDLE = plot_LFPs_by_channel(LFP_TRACES,SPLIT_DIM,LFP_TIMESTAMPS,SPACING,ALPHA_LEVEL)
% LFP plotting function that shows the LFP response by channel, with an
% offset between channels determined by 'spacing'
% 
% LFP_TRACES: An N_channels (by N_trials) by N_samples matrix of LFP data.
% If LFP_traces is a 3D matrix, it is assumed to be M x N x N_Samples
% SPLIT_DIM will determine whether to plot M traces (SPLIT_DIM == 1) or 
% N traces (SPLIT_DIM == 2), and average over the other dimension.
% 
% If LFP_TRACES is a 2D matrix, it is assumed to be M x N_Samples.
% 
% LFP_TIMESTAMPS: The time stamps for the values in LFP_traces. Default is
% 1:length(LFP_TRACES).
%
% SPACING: Specifies how far traces are separated, default value is 0.1. 
% Spacing is a function of the range of values in the traces, so a given 
% spacing should give a somewhat consistent visualisation even for LFP 
% signals with a larger or smaller signal.
% 
% ALPHA_LEVEL: specifies the alpha level (transparency) of the traces,
% default is 0.3.


% Set default 'time stamps' to x = 1:n_samples if not specified by user
if nargin < 2
    LFP_timestamps = 1:length(LFP_traces);
end

if nargin < 3
    split_dim   = 1;
end

% Set default spacing to 1
if nargin < 4
    spacing     = 0.1;
end

if nargin < 5
    alpha_level = .3;
end

LFP_traces          = squeeze(LFP_traces);

avg_dim             = 3 - split_dim; % what dimension to average over; this should be 2 if split_dim is 1 and vice versa

if ndims(LFP_traces) > 2
    LFP_traces  = mean(LFP_traces,avg_dim);
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


