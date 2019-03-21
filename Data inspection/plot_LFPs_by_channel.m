function fighandle = plot_LFPs_by_channel(LFP_traces,LFP_timestamps,spacing)
% function fighandle = plot_LFPs_by_channel(LFP_traces,spacing)
% LFP plotting function that shows the LFP response by channel, with an
% offset between channels determined by 'spacing'
% 
% LFP traces is an N_channels by N_samples matrix of LFP data
% spacing is a value ranging from 0 (traces overlap completely) to 1 ()



LFP_traces      = LFP_traces';

LFP_traces      = notch_filt(LFP_traces,1000,50);

median_LFP_min  = median(min(LFP_traces));
median_LFP_max  = median(max(LFP_traces));

LFP_range       = median_LFP_max - median_LFP_min; 

LFP_traces      = LFP_traces / LFP_range;

LFP_offsets  	= meshgrid((size(LFP_traces,2)+1)-(1:size(LFP_traces,2)),1:size(LFP_traces,1));

LFP_offset_traces   = LFP_traces/spacing + LFP_offsets;

% (1:size(LFP_offset_traces,1))/1000,LFP_offset_traces

plot(LFP_timestamps,LFP_offset_traces,'LineWidth',2)
axis tight
ylimz   = ylim;
ylim([ylimz(1)-spacing ylimz(2)+spacing])
ylabel('Channel number')
xlabel('Time(s)')
set(gca,'FontName','Helvetica','FontSize',24,'LineWidth',2)
