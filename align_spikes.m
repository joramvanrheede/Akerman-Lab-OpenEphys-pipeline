function aligned_traces = align_spikes(data, thresh_val, range, qplot, numtraces)

% function plot_aligned_spikes(datafile, thresh_val, plotrange, numtraces)
% 
% Plots spikes detected in signal aligned by negative local peak
% 
% INPUTS:
% datafile:         Full path to openephys '.continuous' data file
% thresh_val:       Threshold in standard deviations of the signal
% range:            Sample range for aligned traces (samples from -range to +range are included)
% qplot:            Plot?
% numtraces:        How many traces to plot (numeric value, or string 'all')
% 
% OUTPUT:
% aligned_traces:   An M x N matrix, where M is the number of traces, and N is 2*plotrange+1
%

% Create range indices
samplerange                 = [(-2*range):(2*range)];   % sample waveform data over larger range than ploting range
peakrange                 	= [0:10];                   % which samples to consider for determining negative peak for aligning, counting from the sample that crossed threshold

% Create band pass filter object, passing 500Hz-6kHz, and filter data
[filt_b,filt_a]             = butter(2, [500 6000]/(30000/2));
data                        = filter(filt_b,filt_a,data);

% Determine threshold
thresh                      = std(data) * thresh_val;

% Find threshold crossings (negative)
spikes                      = find(diff(-data > thresh) == 1);

% Create matrix of indices to sample a range of data points around the
% threshold crossings
[a,b]                       = meshgrid(samplerange,spikes);
spikesampleinds             = a+b;

% Set out-of-range values to min and max range
spikesampleinds(spikesampleinds < 1) = 1;
spikesampleinds(spikesampleinds > length(data)) = length(data);

% Retrieve relevant data from trace
tracedata                   = data(spikesampleinds);

% Find which nr is the middle of the range
thresh_cross_ind            = ceil(length(samplerange)/2);  

% Get subselection of data that is relevant for finding the peak
testdata                    = tracedata(:,thresh_cross_ind + peakrange);

% Find minimum (= negative peak)
spikemin                    = min(testdata,[],2);

% Create boolean variable representing the occurrence of the peak
spikeminmat                 = repmat(spikemin,1,length(peakrange));
minloc                      = testdata == spikeminmat;

% Find the indices of peaks (annoyingly, sorted by x as per 'find' function)
[ymininds,xmininds]         = find(minloc);

% Re-sort the indices such that they are in trace order
[sortyinds sortinds]        = sort(ymininds);
sortxinds                   = xmininds(sortinds);

% Make grid of indices to retrieve traces
[a,b]                       = meshgrid(samplerange,sortxinds);
xindmat                     = a+b+thresh_cross_ind; % also add back original offset

% Take care of out-of-range values
xindmat(xindmat < 1)        = 1;
xindmat(xindmat > length(samplerange)) = length(samplerange);

% Generate matrix of y-indices
yindmat                     = repmat(sortyinds,1,size(xindmat,2));

% Convert sub-indices to linear indices
lin_trace_inds              = sub2ind(size(tracedata),yindmat(:),xindmat(:));

% Retrieve the aligned traces from tracedata
aligned_traces              = reshape(tracedata(lin_trace_inds),size(xindmat));

% Now only get the part of the trace that was requested, according to range
range_select                = [2 * range + (-range:range)];
aligned_traces              = aligned_traces(:,range_select);

% If plot is requested, make plots
if qplot
    
    % If fewer spikes have been detected than number of traces requested
    % for plotting, plot all available traces
    % Also plot all traces if 'all' have been requested
    if strcmp(numtraces,'all')
        numtraces = size(aligned_traces,1);
    elseif numtraces > size(aligned_traces,1)
        numtraces = size(aligned_traces,1);
    end
    
    % Plot required traces in black with alpha level .05
    plot(-range:range,aligned_traces(1:numtraces,:)','LineWidth',3,'Color',[0 0 0 .05])
    hold on
    
    % Overlay plot of mean trace in white
    plot(-range:range,mean(aligned_traces),'LineWidth',2,'Color',[1 1 1])
    
    % Aesthetic stuff
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16,'FontWeight','bold')
    ylabel('Value')
    xlabel('Sample number')
    title('Spikes aligned by (negative) peak')
    axis tight
    hold off
end







