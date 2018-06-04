% single experiment plots for frequency experiment

close all

experiment          = sdata(2).expt(1);

summarise_channels  = [1:16]; % include these channels
split_conditions    = [1 6]; % split by these conditions, summarise over others

split_plots         = [1 2]; % [4 6] works


%% 

condition_mat       = experiment.condition_mat;

split_cond_mat      = condition_mat(:,split_conditions);
[split_cond_rows, indxa, cond_inds] = unique(split_cond_mat,'rows');

split_plot_mat      = condition_mat(:,split_plots);
[split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');

[a,b,cond_plot_inds] = unique(cond_plot_inds);

figure
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .1 1 .9])
for a = 1:size(split_cond_rows,1)
    
    these_conds         = split_cond_rows(a,:);
    sum_inds            = cond_inds == a;
    
    this_whisk_psth  	= mean(experiment.whisk_win_rel(sum_inds,:,:),1);
    
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),5),'LineWidth',2)
    xlim([-.5 1])
    ylim([0 8])
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
    
end

set(gcf,'Color',[1 1 1])


% determine peak IROF

%
contrast_conds      = [1 6];
contrast_cond_mat   = condition_mat(:,contrast_conds);
[contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
contrast_rates      = [];
contrast_times      = [];
for a  = unique(cond_inds)'
    
    contrast_rates  = [contrast_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
    contrast_times  = [contrast_times; mean(mean(experiment.whisk_peak_time(cond_inds == a,summarise_channels)))];
end

if contrast_rates(end-1) > contrast_rates(end)
    P_whisk_stim = 1;
    A_whisk_stim = 2;
else
    P_whisk_stim = 2;
    A_whisk_stim = 1;
end

P_rates             = contrast_rates(P_whisk_stim:2:end);
A_rates             = contrast_rates(A_whisk_stim:2:end);

PA_ratios           = P_rates ./ A_rates;
LED_delays       	= unique(condition_mat(:,1) - condition_mat(:,2));

figure
title('Peak instantaneous rate of fire')

% plots
plot(LED_delays(1:end-1),P_rates(1:end-1),'k.-','LineWidth',4,'MarkerSize',30)
hold on
plot(LED_delays(1:end-1),A_rates(1:end-1),'r.-','LineWidth',4,'MarkerSize',30)

xlim([(LED_delays(1) - .05) (LED_delays(end-1) + .05)])

line(xlim,[P_rates(end) P_rates(end)],'LineWidth',4,'Color',[.5 .5 .5],'LineStyle','--');
line(xlim,[A_rates(end) A_rates(end)],'LineWidth',4,'Color',[1 .5 .5],'LineStyle','--');


legend({'Principal whisker' 'Adjacent whisker', 'Principal baseline', 'Adjacent baseline'},'Location','SouthEast')

ylimits = ylim;
ylim([0 ylimits(2)]);
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16,'FontWeight','bold')
set(gcf,'Color',[1 1 1])

ylabel('Whisker response after LED stimulus at time X')
xlabel('Time')
title('Effect of timed LED pulse on whisker response')

%%
figure
plot(LED_delays(1:end-1),PA_ratios(1:end-1),'k.-','LineWidth',4,'MarkerSize',20)
hold on
xlim([(LED_delays(1) - .05) (LED_delays(end-1) + .05)])
line(xlim,[PA_ratios(end) PA_ratios(end)],'LineWidth',4,'Color',[.5 .5 .5],'LineStyle','--');
ylimits = ylim;
ylim([0 ylimits(2)]);
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16,'FontWeight','bold')
set(gcf,'Color',[1 1 1])
legend({'P/A ratio', 'P/A ratio baseline'},'Location','SouthEast')
title('Effect of timed LED pulse on P/A response ratio')

