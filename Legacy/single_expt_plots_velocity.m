% single experiment plots for velocity experiment

close all

experiment          = sdata(11).expt(1);

summarise_channels  = [1:16]; % include these channels
split_conditions    = [1 5 6]; % split by these conditions, summarise over others

split_plots         = [5 6]; % [4 6] works

% for rasterplots:
split_raster_plots  = [1 5];    % split plots by these conditions
split_figures       = [6];      % split


P_whisk_stim        = 1;
A_whisk_stim        = 2;


qchannelim          = true;
qraster             = true;
save_figs           = true;

%%

if save_figs
    fileseps = regexp(experiment.filename,'/.');
    experiment_plot_folder      = [save_folder filesep experiment.filename((fileseps(end)+1):end-4)];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
end

%% Make PSTH plots split by condition

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
    
    this_whisk_psth  	= mean(experiment.whisk_win_rates(sum_inds,:,:),1);
    
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),5),'LineWidth',2)
    xlim([-.1 1])
    ylim([0 750])
    hold on
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    
end
set(gcf,'Color',[1 1 1])
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')

if save_figs
    print(gcf,[experiment_plot_folder filesep 'PSTH plots ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end

%% Look at peak instantaneous firing rate between LED and no LED conditions
contrast_conds      = [1 6];
contrast_cond_mat   = condition_mat(:,contrast_conds);
[contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
contrast_rates      = [];
contrast_times      = [];
for a  = unique(cond_inds)'
    
    contrast_rates  = [contrast_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
    contrast_times  = [contrast_times; mean(mean(experiment.whisk_peak_time(cond_inds == a,summarise_channels)))];
end

P_rate_LED_on       = contrast_rates(P_whisk_stim);
P_rate_LED_off      = contrast_rates(P_whisk_stim+2);
A_rate_LED_on       = contrast_rates(A_whisk_stim);
A_rate_LED_off      = contrast_rates(A_whisk_stim+2);

PA_ratio_LED_on     = P_rate_LED_on / A_rate_LED_on;
PA_ratio_LED_off    = P_rate_LED_off / A_rate_LED_off;

P_LED_onoff_ratio   = P_rate_LED_on / P_rate_LED_off;
A_LED_onoff_ratio 	= A_rate_LED_on / A_rate_LED_off;

% plotting
figure
bar_handle  = bar([P_rate_LED_on P_rate_LED_off; A_rate_LED_on A_rate_LED_off],'LineWidth',2);
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','bold')
set(gcf,'Color',[1 1 1])
title('Response size, P vs. A whisker, LED ON vs. OFF')
ylabel('Peak stimulus-evoked firing rate (Hz)')
xlabel('Whisker (1 = principal, 2 = adjacent)')
legend({'LED ON' 'LED OFF'})

%% rasterplots
if qraster
    % use channels_to_rasterplots figure to generate
    fighandles = channels_to_rasterplots(experiment.filename,summarise_channels,split_raster_plots,split_figures,1:8,x_ax_lims,condition_name, condition_units);
    
    % save figures
    if save_figs
        for a = 1:length(fighandles)
            figure(fighandles(a))
            print(gcf,[experiment_plot_folder filesep 'Rasterplot stimulator ' num2str(a) ' chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
        end
    end
end

%% visual representation of responses across channels
if qchannelim
    whisk_resps     = experiment.whisk_peak_rate;
    stimulator      = experiment.condition_mat(:,6);
    LEDtime         = experiment.condition_mat(:,1) < 2;
    
    spont_rates     = experiment.spont_rate;
    
    P_whisk_resps   = whisk_resps(P_whisk_stim:2:end,:);
    A_whisk_resps   = whisk_resps(A_whisk_stim:2:end,:);
    
    P_whisk_resps_LED   = whisk_resps(stimulator == P_whisk_stim & LEDtime,:);
    P_whisk_resps_nLED  = whisk_resps(stimulator == P_whisk_stim & ~LEDtime,:);
    
    A_whisk_resps_LED   = whisk_resps(stimulator == A_whisk_stim & LEDtime,:);
    A_whisk_resps_nLED  = whisk_resps(stimulator == A_whisk_stim & ~LEDtime,:);
    
    LED_resps           = experiment.LED_rate;
    
    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.2 .2 .6 .6])

    subplot(1,8,1)
    imagesc(P_whisk_resps_LED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('PW + LED')
    ylabel('Channel number')
    % colorbar
    axis off
    
    subplot(1,8,2)
    imagesc(P_whisk_resps_nLED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('PW - LED')
    ylabel('Channel number')
    % colorbar
    axis off
    
    subplot(1,8,3)
    imagesc(A_whisk_resps_LED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('AW + LED')
    %colorbar
    axis off
    
    subplot(1,8,4)
    imagesc(A_whisk_resps_nLED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('AW - LED')
    % colorbar
    axis off
    
    subplot(1,8,5)
    imagesc(log(P_whisk_resps_LED' ./ A_whisk_resps_LED'))
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('PAR + LED')
    % colorbar
    axis off
    
    subplot(1,8,6)
    imagesc(log(P_whisk_resps_nLED' ./ A_whisk_resps_nLED'))
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('PAR - LED')
    %colorbar
    axis off
    
    subplot(1,8,7)
    imagesc(LED_resps')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('LED')
    % colorbar
    axis off
    
    subplot(1,8,8)
    imagesc(spont_rates(:))
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Spont')
    ylabel('Channel number')
    % colorbar
    axis off
    
    PA_ratio = mean(A_whisk_resps(:)) / mean(P_whisk_resps(:))
    
    [r_P_whisk_LED, p_P_whisk_LED] = corr(mean(LED_resps)',mean(P_whisk_resps)','type','Pearson')
    
    [r_A_whisk_LED, p_A_whisk_LED] = corr(mean(LED_resps)',mean(A_whisk_resps)','type','Pearson')
    
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'Channel heatmap' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
    end
    
end
