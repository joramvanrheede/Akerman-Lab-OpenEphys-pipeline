function image_handle = LFP_heatmap(LFP_traces)
% function image_handle = LFP_heatmap(LFP_traces)
% 
% Plot matrix of LFP data (N_CHANNELS x N_TRIALS x N_TIME_POINTS) as a 
% heatmap with channels on the Y axis and time on the X axis.
% 
% 


LFP_traces          = mean(LFP_traces,2);
LFP_traces          = squeeze(LFP_traces);

% Plot
image_handle        = imagesc(LFP_traces);

colormap(LFP_colormap)

% Set symmetrical colour limits so 0 is the middle colour of the colourmap
max_val             = robust_max(LFP_traces,0.5,'all');
min_val             = robust_min(LFP_traces,0.5,'all');

range_val           = max(max_val, abs(min_val));

set(gca,'CLim',[-range_val range_val])

colorbar

% Plot labeling
ylabel('Trace number')
xlabel('Time')


