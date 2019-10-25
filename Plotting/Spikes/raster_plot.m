function line_handles = raster_plot(spike_data,split_dim,spacing)
% function LINE_HANDLES = RASTER_PLOT(SPIKE_DATA,SPLIT_DIM,SPACING)
% 
% SPIKE_DATA: M*N*O matrix of spike times. In current Akerman lab
% 'ephys_data' format, the dimensions represent:
% M = number of channels
% N = number of trials
% O = number of spikes
% 
% 
% NaNs are ignored (it is assumed that NaNs represent empty slots in the 
% SPIKE_DATA matrix.
% 
%
% SPLIT_DIM: Determines which dimension of SPIKE_DATA to display as individual 
% rows, and which to pool together. E.g. for SPLIT_DIM == 1, output is a
% row for each channel and pools spikes over N trials. SPLIT_DIM == 2 gives
% a row for each individual trial and pools spikes over M channels.
% (the only options are 1 or 2)
% 
% SPACING: Determines the space between rows. Default for a sensible
% visualisation is 0.2.
% 
% LINE_HANDLES: Returns the line handles for the line objects in the raster
% plot.
% 
% Joram van Rheede - Jan 2019


% Default to displaying the data by channel, pooling over N trials
if nargin < 2
    split_dim = 1;
end

% Default vertical spacing between raster plot lines
if nargin < 3
    spacing = 0.2;
end

% Determine how many rows the raster plot will have
n_rows  = size(spike_data, split_dim);

% Plot empty point to make sure we overwrite existing plots if 'hold' is not 'on'
plot([0],[0],'w-')

% For each row...
for a = 1:n_rows
    
    % Get relevant spikes, splitting by split_dim and pooling over other dims
    if split_dim == 1
        trial_spikes      	= spike_data(a,:,:);
    elseif split_dim == 2
        trial_spikes      	= spike_data(:,a,:);
    end
    
    % turn into an N spikes * 1 vector
    xvals               = trial_spikes(:);
    
    % Generate a vector of X values that goes like this:
    % [X1 X1 NaN X2 X2 NaN ... Xn Xn NaN]
    % Each set of X values is used to draw a vertical line (where the 
    % X-position does not vary, and the NaN is used to generate an empty
    % line segment so that the line for one spike is not connected to the 
    % next spike. This allows for using a single line object for each row
    % of spikes and is much faster and memory-efficient than generating
    % a line object for each individual spike.
    x_raster            = NaN(length(xvals)*3,1);
    x_raster(1:3:end)   = xvals;
    x_raster(2:3:end)   = xvals;
    
    % Generate a vector of Y values that goes like this:
    % [Y1a Y1b NaN Y2a Y2b NaN ... YNa YNb NaN]
    % Each set (YNa&b) specifies Y values for the top and bottom of the
    % line within the relevant row; interspersing with NaNs generates empty
    % line segment for efficiency reasons outlined above.
    y_raster            = NaN(length(xvals)*3,1);
    y_raster(1:3:end)   = a-spacing/2;
    y_raster(2:3:end)   = a-(1-spacing/2);
    
    % Plot all spikes in this row as a single line object with empty segments
    line_handles(a)     = line(x_raster(:),y_raster(:),'Color',[0 0 0]);
    
end

% Because your spikes are tight AF:
axis tight

% 'line' function by default seems to leave 'hold' set to 'on'
hold off

% invert y axis so that figure can be read from top to bottom (i.e. top row
% is the first / most superficial channel, or the first trial / sweep)
axis ij


