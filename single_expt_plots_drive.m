% single experiment plots for frequency experiment

close all

experiment          = sdata(5).expt;

summarise_channels  = [1:16]; % include these channels
split_conditions    = [1 5 6]; % split by these conditions, summarise over others

split_plots         = [5 6]; % [4 6] works

save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
save_figs        	= false;

%% 

condition_mat       = experiment.condition_mat;

split_cond_mat      = condition_mat(:,split_conditions);
[split_cond_rows, indxa, cond_inds] = unique(split_cond_mat,'rows');

split_plot_mat      = condition_mat(:,split_plots);
[split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');

[a,b,cond_plot_inds] = unique(cond_plot_inds);

figure
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .6])
for a = 1:size(split_cond_rows,1)
    
    these_conds         = split_cond_rows(a,:);
    sum_inds            = cond_inds == a;
    
    this_whisk_psth  	= mean(experiment.whisk_win_rel(sum_inds,:,:),1);
    
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),5),'LineWidth',2)
    xlim([-1 1.5])
    ylim([0 10])
    ylabel('Spike rate (relative to spontaneous)')
    xlabel('Time (s)')
    hold on
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',30)
    
end
set(gcf,'Color',[1 1 1])
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
legend({'LED on' 'LED off'})
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')
legend({'LED on' 'LED off'})

if save_figs
    experiment_plot_folder      = [save_folder filesep experiment.filename(1:end-4)];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
    print(gcf,[experiment_plot_folder filesep ' PSTH plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end

%% 


%% Look at peak instantaneous firing rate between LED and no LED
contrast_conds      = [1 6];
contrast_cond_mat   = condition_mat(:,contrast_conds);
[contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
contrast_rates      = [];
contrast_times      = [];
for a  = unique(cond_inds)'
    
    contrast_rates  = [contrast_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
    contrast_times  = [contrast_times; mean(mean(experiment.whisk_peak_time(cond_inds == a,summarise_channels)))];
end

% Find out which whisker is principal and which is adjacent
if contrast_rates(3) > contrast_rates(4)
    P_whisk_stim = 1;
    A_whisk_stim = 2;
else
    P_whisk_stim = 2;
    A_whisk_stim = 1;
end

P_rate_LED_on       = contrast_rates(P_whisk_stim);
P_rate_LED_off      = contrast_rates(P_whisk_stim+2);
A_rate_LED_on       = contrast_rates(A_whisk_stim);
A_rate_LED_off      = contrast_rates(A_whisk_stim+2);

PA_ratio_LED_on     = P_rate_LED_on / A_rate_LED_on;
PA_ratio_LED_off    = P_rate_LED_off / A_rate_LED_off;

P_LED_onoff_ratio   = P_rate_LED_on / P_rate_LED_off;
A_LED_onoff_ratio 	= A_rate_LED_on / A_rate_LED_off;

%% 
figure
bar_handle  = bar([P_rate_LED_on P_rate_LED_off; A_rate_LED_on A_rate_LED_off],'LineWidth',2);
legend(bar_handle,'LED on','LED off')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
ylim([0 1000])
title('Response size, P vs. A whisker, LED ON vs. OFF')
ylabel('Stimulus-evoked instantaneous firing rate')
xlabel('1 = Principal, 2 = Adjacent')

if save_figs
    print(gcf,[experiment_plot_folder filesep ' Bar graph P vs A - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end

