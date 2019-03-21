function aligned_traces = align_spikes(datafile, thresh_val, plotrange, numtraces)

% function plot_aligned_spikes(datafile, thresh_val, plotrange, numtraces)
% 
% Plots spikes detected in signal aligned by negative local peak
% 
% INPUTS:
% datafile:         Full path to openephys '.continuous' data file
% thresh_val:       Threshold in standard deviations of the signal
% plotrange:       	Sample range for plotting spikes (samples from -range to +range are plotted)
% numtraces:        How many traces to plot
% 
% OUTPUT:
% aligned_traces:   An M x N matrix, where M is the number of traces, and N is 2*plotrange+1
%

plotrange       = [-plotrange:plotrange];
samplerange     = [(-2*plotrange):2*plotrange];

[data, timestamps, info]    = load_open_ephys_data(datafile);
[filt_b,filt_a]             = butter(2, [500 6000]/(30000/2));

data    = filter(filt_b,filt_a,data);

thresh          = std(data) * thresh_val;

spikes          = find(diff(-data > thresh) == 1);

[a,b]           = meshgrid(range,spikes);

spikesampleinds = a+b;

spikesampleinds(spikesampleinds < 1) = 1;
spikesampleinds(spikesampleinds > length(data)) = length(data);

tracedata       = data(spikesampleinds);

spikemin        = min(tracedata,[],2);

spikeminmat     = repmat(spikemin,1,length(range));

minloc          = tracedata == spikeminmat;

[ymininds,xmininds]     = find(minloc);

[sortyinds sortinds]    = sort(ymininds);
sortxinds               = xmininds(sortinds);

% q_xinds                 = sortxinds > 7 & sortxinds < 22;

% tracedata               = tracedata(q_xinds,:);
% xinds                   = sortxinds(q_xinds);

[a,b]                   = meshgrid(plotrange,sortxinds);

xindmat                 = a+b;

xindmat(xindmat < 1) = 1;
xindmat(xindmat > length(range)) = length(range);

yindmat = repmat(sortyinds,1,size(xindmat,2));

lin_trace_inds          = sub2ind(size(tracedata),yindmat(:),xindmat(:));

aligned_traces          = reshape(tracedata(lin_trace_inds),size(xindmat));



if numtraces > size(aligned_traces,1)
    numtraces = size(aligned_traces,1);
end

plot(aligned_traces(1:numtraces,:)','LineWidth',3,'Color',[0 0 0 .05])
hold on
plot(mean(aligned_traces),'LineWidth',2,'Color',[1 1 1])
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16,'FontWeight','bold')
ylabel('Value')
xlabel('Sample number')
title('Spikes aligned by (negative) peak')
hold off









