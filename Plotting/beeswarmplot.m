function plothandle = beeswarmplot(indata,groups,labels,yaxislabel)
% function plothandle = beeswarmplot(indata,groups,labels,yaxislabel)
% 

% make histogram of data
% 1: determine right bin size of histogram (max 5-8 data points per bin)
% 
% make histogram, 
% determine where in histogram each data point falls
% for each bin in the histogram, determine optimal spreading
% e.g. 2 data points: 0/1, 1/1, 3 data points: 0/1, .5/1, 1/1
% 

% locat
% testn = 50;
% 
% randdata        = rand(testn,100);
% indata          = mean(randdata,2)
% indata          = rand(testn,1) .* rand(testn,1) + rand(testn,1) %./ rand(100,1);
% labels              = {'group1' 'group2'}
%

figure
dotwidth            = 1/3;
meanwidth           = 1/3;

xincrement          = 1;
xshift              = 0;
groupids            = unique(groups);
for b = 1:length(groupids)
    b
    % shift beeswarm along X axis 
    xshift          = xshift + xincrement;
    
    % find data column
    thiscol         = indata(groups == groupids(b));
    
    % determine a sensible number of bins for the data
    nbins           = ceil(range(thiscol)/serr(thiscol))/2;
    
    % 
    [counts, bins]  = histc(thiscol,linspace(min(thiscol),max(thiscol),nbins));
    
    xvals           = NaN(size(thiscol));
    
    for a = 1:length(counts)
        
        nevents             = counts(a);
        
        locations           = linspace(-1,1,nevents + 2);
        locations           = locations(2:end-1);
        
        locations           = locations * dotwidth;
        
        datainds            = find(bins == a);
        
        for b = 1:length(datainds)
            xvals(datainds)     = locations;
        end
    end
    
    xvals       = xvals + xshift;

    plothandle  = plot(xvals,thiscol,'ko','MarkerSize',8,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor', [.5 .5 .5]);
    
    hold on
    
    line([(xshift - meanwidth) (xshift + meanwidth) ], [nanmedian(thiscol) nanmedian(thiscol)],'LineWidth',3,'Color',[0 0 0])
    
    xshift+xincrement
    
    xlim([0 length(groupids)+1])
    
    scf
    
end

xlim([0 length(groupids)+1])
set(gca,'XTick',1)
set(gca,'XTicklabel',{'' labels{:} ''})
set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','bold')
set(gca,'LineWidth',2)
ylabel(yaxislabel)
set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','bold')
xlabel('Group')
set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','bold')





hold off
