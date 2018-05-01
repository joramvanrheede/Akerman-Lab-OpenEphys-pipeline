% single experiment plots for frequency experiment

close all

experiment          = sdata(11).expt(1);

summarise_channels  = [1:16]; % include these channels
split_conditions    = [1 4 6]; % split by these conditions, summarise over others
split_plots         = [4 6]; % [4 6] works

% for rasterplots:
split_raster_plots  = [4 1];
split_figures       = [6];
trialrange          = [1 10]; % [min max] - don't exceed max nr of trials; else errors result.


save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
save_figs        	= false;

respwinsize         = .04;

%% Make PSTH plots split by condition

condition_mat       = experiment.condition_mat;
frequencies         = condition_mat(:,4);

split_cond_mat      = condition_mat(:,split_conditions);
[split_cond_rows, indxa, cond_inds] = unique(split_cond_mat,'rows');

split_plot_mat      = condition_mat(:,split_plots);
[split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');

[a,b,cond_plot_inds] = unique(cond_plot_inds);

figure(1)
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .6])
figure(2)
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .6])
figure(3)
set(gcf,'Units','normalized')
set(gcf,'Position',[.3 0 .4 1])
for a = 1:size(split_cond_rows,1)
    
    these_conds         = split_cond_rows(a,:);
    sum_inds            = cond_inds == a;
    
    this_whisk_psth  	= mean(experiment.whisk_win_rel(sum_inds,:,:),1);
    
    figure(1)
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),5),'LineWidth',2)
    hold on
    this_frequency      = mean(frequencies(sum_inds));
    
    stim_times          = [0:1/this_frequency:5.5];
    
    plot(stim_times,zeros(size(stim_times))+0.2,'k^','MarkerSize',5,'LineWidth',2)
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    
    xlim([-1 4.5])
    ylim([0 10])
    
    
    profile_plots(a).stim_times     = stim_times;
    for b = 1:length(stim_times)
        qstimwin            = experiment.whiskwinedges > stim_times(b) & experiment.whiskwinedges < (stim_times(b) + respwinsize);
        profile_segment     = squeeze(mean(mean(experiment.whisk_profile(sum_inds,summarise_channels,qstimwin),1),2));
        if isempty(profile_segment)
            profile_segment = NaN;
        end
        
        profile_plots(a).resp_peaks(b)  = max(profile_segment);
        
    end
    
    figure(2)
    ylim([0 .3])
    xlim([0 4])
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(profile_plots(a).stim_times,profile_plots(a).resp_peaks,'LineWidth',2,'MarkerSize',25);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
    
    figure(3)
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    ylim([0 .3])
    xlim([0 6])
    profile_plots(a).resp_peaks = profile_plots(a).resp_peaks(~isnan(profile_plots(a).resp_peaks));
    plot(1:5,[profile_plots(a).resp_peaks([1:4]) profile_plots(a).resp_peaks(end)],'LineWidth',2,'MarkerSize',25);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
end

figure(1)
set(gcf,'Color',[1 1 1])
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')

figure(2)
set(gcf,'Color',[1 1 1])
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')

figure(1)
if save_figs
    experiment_plot_folder      = [save_folder filesep experiment.filename(1:end-4)];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
    print(gcf,[experiment_plot_folder filesep ' PSTH plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end

%% rasterplots

channels_to_rasterplots(experiment.filename,summarise_channels,split_raster_plots,split_figures,trialrange,x_ax_lims)


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

% plot

figure
bar_handle  = bar([P_rate_LED_on P_rate_LED_off; A_rate_LED_on A_rate_LED_off],'LineWidth',2);
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','bold')
ylim([0 1000])
title('Response size, P vs. A whisker, LED ON vs. OFF')
ylabel('Peak stimulus-evoked firing rate (Hz)')
xlabel('Whisker (1 = principal, 2 = adjacent)')
legend({'LED ON' 'LED OFF'})

if save_figs
    print(gcf,[experiment_plot_folder filesep ' Bar graph P vs A - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end

