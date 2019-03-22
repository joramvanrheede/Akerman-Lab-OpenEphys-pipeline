function beeswarmplot(indata,groups,labels,colours)
% function beeswarmplot(INDATA,GROUPS,LABELS,COLOURS)
% 
% INDATA: vector of input data pounts
% GROUPS: vector of same length as INDATA, indicating group membership
% LABELS: Optional group labels (default is 'Group 1', 'Group 2', etc)
% COLOURS: Optional colour specification for each group, as an N GROUPS * 3
% RGB specification

% Some constants that govern the size of the spread of the dots
dotwidth            = 1/3;
meanwidth           = 1/3;

% Find number of unique groups
groupids            = unique(groups);

% Set default labels
if ~exist('labels','var') || isempty(labels)
    for a = 1:length(groupids)
        labels{a} = ['Group ' num2str(groupids(a))];
    end
end

% Set default colours
if ~exist('colours','var')
    colours = lines;
end

% Open main figure for plotting
fig_handle = figure;

xshift              = 0;
for b = 1:length(groupids)

    % shift beeswarm along X axis 
    xshift          = xshift + 1;
    
    % select data for this group
    thiscol         = indata(groups == groupids(b));
    
    % determine a sensible number of bins for the data
    nbins           = ceil(range(thiscol)/serr(thiscol))/2;
    
    % count number of data points that fall within these bins to determine
    % spread of points at this location
    [counts, bins]  = histc(thiscol,linspace(min(thiscol),max(thiscol),nbins));
    
    % start generating x values for points
    xvals           = NaN(size(thiscol));
    for a = 1:length(counts)
        
        nevents             = counts(a); % How many data points are we dealing with in this y space
        
        locations           = linspace(-1,1,nevents + 2); % generate
        locations           = locations(2:end-1);
        
        locations           = locations * dotwidth;
        
        datainds            = find(bins == a);
        
        for c = 1:length(datainds)
            xvals(datainds)     = locations;
        end
    end
    
    % Shift x values according to the group
    xvals       = xvals + xshift;
    
    % Use plot function to plot the points
    plothandle  = plot(xvals,thiscol,'ko','MarkerSize',8,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor', colours(b,:));
    
    hold on
    
    % Plot median line
    line([(xshift - meanwidth) (xshift + meanwidth) ], [nanmedian(thiscol) nanmedian(thiscol)],'LineWidth',3,'Color',[0 0 0])
    
    % Set sensible x extent for axes
    xlim([0 length(groupids)+1])
    
end

xlim([0 length(groupids)+1])
set(gca,'XTick',1:length(groupids),'XTicklabel',labels)

hold off
