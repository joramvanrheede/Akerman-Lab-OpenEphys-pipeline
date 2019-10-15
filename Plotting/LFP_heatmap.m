function plot_handle = plot_LFP_traces(LFP_traces,LFP_timestamps)
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
    LFP_timestamps = (1:length(LFP_traces));
end

LFP_traces          = mean(LFP_traces,2);
LFP_traces          = squeeze(LFP_traces);

image_handle        = imagesc(LFP_traces);

colormap(LFP_colormap)

% Plot labeling
ylabel('Trace number')
xlabel('Time')

% To do - create scale bar?


