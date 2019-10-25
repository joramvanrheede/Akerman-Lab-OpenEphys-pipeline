function plot_handle = plot_LFP_traces(LFP_traces,split_dim,LFP_timestamps,spacing,alpha_level)
% function PLOT_HANDLE = plot_LFP_traces(LFP_TRACES,SPLIT_DIM,LFP_TIMESTAMPS,SPACING,ALPHA_LEVEL)
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
% SPACING: Specifies how far traces are separated, default value is 0; this
% is a special case which plots the overlaid individual traces and then a 
% mean trace on top. 
% Any values of spacing greater than 0 will plot the 
% 
% Spacing is a function of the range of values in the traces, so a given 
% spacing should give a somewhat consistent visualisation even for LFP 
% signals with a larger or smaller signal.
% 
% ALPHA_LEVEL: specifies the alpha level (transparency) of the traces,
% default is 0.2.


% Set default dimension for splitting traces (and take average over the other)
if nargin < 2
    split_dim   = 1;
end

% Set default 'time stamps' to x = 1:n_samples if not specified by user
if nargin < 3
    LFP_timestamps = 1:length(LFP_traces);
end

% Set default spacing to 0 to plot the traces with mean overlaid
if nargin < 4
    spacing     = 0;
end

% Set default alpha level
if nargin < 5
    if spacing == 0
        alpha_level = .2; % Default for overlaid traces
    else
        alpha_level = .4; % Default for separated traces
    end
end

LFP_traces          = squeeze(LFP_traces);

avg_dim             = 3 - split_dim; % what dimension to average over; this should be 2 if split_dim is 1 and vice versa

if ndims(LFP_traces) > 2
    LFP_traces  = squeeze(mean(LFP_traces,avg_dim));
end

LFP_traces          = LFP_traces'; % LFP_traces needs to be in a different orientation for the functions to follow

if spacing ~= 0
    % Determine the 'median range' of the traces which is used for determining the spacing between traces
    median_LFP_min      = median(min(LFP_traces)); % for each trace, determine the max, and then the median of all trace maxima
    median_LFP_max      = median(max(LFP_traces)); % for each trace, determine the min, and then the median of all trace minima
    LFP_range           = median_LFP_max - median_LFP_min;
    
    % Normalise traces by the range so they should now have a median extent of ~1
    LFP_traces          = LFP_traces / LFP_range;
    
    % Invert traces because we will invert the y-axis
    LFP_traces          = -LFP_traces;
    
    % Make a matrix of offsets of the same size as LFP_traces, but that decrements by 1 for every next trace:
    LFP_offsets         = meshgrid((1:size(LFP_traces,2)),1:size(LFP_traces,1));
    
    % Scale the traces for spacing, then add the offsets
    LFP_offset_traces   = LFP_traces/spacing + LFP_offsets;
    
    % Plot the traces
    plot_handle         = plot(LFP_timestamps,LFP_offset_traces,'LineWidth',2,'Color',[0 0 0 alpha_level]);
    
    % Some plot aesthetics
    axis tight
    axis ij
    ylimz   = ylim;
    ylim([ylimz(1)-2*spacing ylimz(2)+2*spacing])
    
    % Y axis label
    ylabel('Trace number')
else
    % spacing == 0; special case - plot overlaid transparent individual traces with mean trace on top
    
    % Plot individual traces with transparency
    plot(LFP_timestamps,LFP_traces,'LineWidth',2,'Color',[0 0 0 alpha_level]); 
    
    hold on
    
    % Overlay mean trace
    plot(LFP_timestamps,mean(LFP_traces,2),'k-','LineWidth',2); % plot mean trace
    
    axis tight
end

% Plot labeling
xlabel('Time')
set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2,'Box','off','TickDir','out','box','off')


