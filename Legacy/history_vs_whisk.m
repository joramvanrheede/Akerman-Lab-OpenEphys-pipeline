
save_folder = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output/Recent activity vs whisk plots'

mode        = 'P'; % ('P' = principal, A = Adjacent)

% look at recent activity history and correlate with whisker response
close all
figure
set(gcf,'Units','Normalized','Color',[1 1 1])
set(gcf,'Position',[.1 .3 .8 .5])

for b = 1:length(sdata)
    
    experiment      = sdata(b).expt;
    
    control_conds  	= [15 16]; % which condition, which whisker?
    
    prewhiskbins    = analysisparams.prewhiskbins;
    
    win_nr          = 1:(length(prewhiskbins)-1);
    chan_nr         = 1:16;
    test            = [];
    for c = 1:length(control_conds)
        test(c).thesewhiskrates  	= experiment.whiskrates(control_conds(c),:);
        test(c).thesewhiskrates 	= cell2mat(test(c).thesewhiskrates');
        
        test(c).meanrate            = mean(test(c).thesewhiskrates(:));
        
    end
    
    meanrates   = [test.meanrate];
    P_ind       = find(meanrates == max(meanrates(:)));
    whiskrates  = test(P_ind).thesewhiskrates;
    
    prewhiskcounts  = experiment.prewhiskcounts(control_conds(P_ind),:);
    prewhiskmat     = [];
    for a = 1:length(prewhiskcounts)
        prewhiskmat(a,:,:) = prewhiskcounts{a};
    end
    
    correlations    = [];
    pvals           = [];
    for a = 1:length(win_nr)
        prewhiskvals   	= squeeze(prewhiskmat(chan_nr,:,a));
        whiskratevals   = whiskrates(:,:);
        
        
        [correlations(a),pvals(a)] = corr(prewhiskvals(:), whiskratevals(:));
    end
    
    subplot(1,2,1)
    plot(prewhiskbins(1:end-1),correlations,'LineWidth',2)
    hold on
    % hold on
    subplot(1,2,2)
    plot(prewhiskbins(1:end-1),log10(1./pvals),'LineWidth',2)
    hold on
end

subplot(1,2,1)
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
title('Correlation between activity at time x and whisk response')
ylabel('Correlation coefficient (r)')
xlabel('Time bin prior to whisker')

subplot(1,2,2)
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
title('Correlation between activity at time x and whisk response')
ylabel('Log 1 over p-value (p)')
xlabel('Time bin prior to whisker')

if ~isdir(save_folder)
    mkdir(save_folder)
end

print(gcf,[save_folder filesep 'History vs whisk resp plot'],'-dpng','-r300')