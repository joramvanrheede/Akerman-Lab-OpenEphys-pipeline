% single experiment plots for frequency experiment

close all

experiment          = sdata(1).expt;

summarise_channels  = [1:16]; % include these channels
split_conditions    = [1 6 7]; % split by these conditions, summarise over others

split_plots         = [6 7]; % [ works

save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
save_figs        	= false;

resp_measure        = 'peak'; % 'peak' or 'area'

qchannelim          = true;
qLEDplot            = true;
qraster             = true;
qpsth               = true;
qrespsizeplot       = true;

% for rasterplots:
split_raster_plots  = [1 7];    % split plots by these conditions
split_figures       = [6];      % split 
% for raster plots:
trialrange          = [1 10]; % [min max] - don't exceed max nr of trials; else errors result.
x_ax_lims           = [0 4]; % limits for x-axes
condition_name    	= 'LED timing';
condition_units   	= 's';

%%
if save_figs
    fileseps = regexp(experiment.filename,'/.');
    experiment_plot_folder      = [save_folder filesep experiment.filename((fileseps(end)+1):end-4)];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
end

%%

condition_mat       = experiment.condition_mat;

split_cond_mat      = condition_mat(:,split_conditions);
[split_cond_rows, indxa, cond_inds] = unique(split_cond_mat,'rows');

split_plot_mat      = condition_mat(:,split_plots);
[split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');

[a,b,cond_plot_inds] = unique(cond_plot_inds);

if qpsth
    figure
    set(gcf,'Units','normalized')
    set(gcf,'Position',[0 .4 1 .6])
    for a = 1:size(split_cond_rows,1)
        
        these_conds         = split_cond_rows(a,:);
        sum_inds            = cond_inds == a;
        
        this_whisk_psth  	= mean(experiment.whisk_win_rates(sum_inds,:,:),1);
        
        subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
        plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),5),'LineWidth',2)
        
        xlim([-1.5 1.5])
        ylabel('Spike rate')
        xlabel('Time (s)')
        hold on
        set(gca,'LineWidth',2,'FontName','Garamond','FontSize',30)
        
    end
    
    % Set all y axes to the same range (based on the largest range)
    plotaxes    = get(gcf,'Children');
    maxy        = cellfun(@max,get(plotaxes,'Ylim'));
    set(plotaxes,'Ylim',[0 max(maxy)]);
    
    set(gcf,'Color',[1 1 1])
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
    title(gca,'Stimulator 1')
    legend({'LED on' 'LED off'})
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
    title(gca,'Stimulator 2')
    legend({'LED on' 'LED off'})
    
    if save_figs
        print(gcf,[experiment_plot_folder filesep ' PSTH plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
    end
end

%% LED

if qLEDplot
    % experiment.LED_win_counts
    
    LED_delays      = experiment.condition_mat(:,1);
    
    control_LED     = LED_delays == max(LED_delays);
    
    LED_rate_hist   = experiment.LED_win_rates(control_LED,summarise_channels,:);
    
    LED_rate_hist   = squeeze(mean(mean(LED_rate_hist,1),2));
    
    figure
    plot(experiment.LEDwinedges(1:end-1),smooth(LED_rate_hist,5),'LineWidth',2)
    set(gca,'FontName','Garamond','FontWeight','Bold','FontSize',16)
    set(gcf,'Color',[1 1 1])
    title('LED response')
    xlim([experiment.LEDwinedges(1) experiment.LEDwinedges(end)])
    xlabel('Time (s) relative to LED stimulus')
    ylabel('Firing rate')
    if save_figs
        print(gcf,[experiment_plot_folder filesep ' LED plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
    end
end


%% rasterplots
if qraster
    % use channels_to_rasterplots figure to generate
    fighandles = channels_to_rasterplots(experiment.filename,summarise_channels,split_raster_plots,split_figures,1:10,x_ax_lims,condition_name, condition_units);
    
    % save figures
    if save_figs
        for a = 1:length(fighandles)
            figure(fighandles(a))
            print(gcf,[experiment_plot_folder filesep 'Rasterplot stimulator ' num2str(a) ' chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
        end
    end
end

%% Look at peak instantaneous firing rate between LED and no LED
contrast_conds      = [1 6];
contrast_cond_mat   = condition_mat(:,contrast_conds);
[contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
contrast_rates      = [];
contrast_times      = [];
for a  = unique(cond_inds)'
    switch resp_measure
        case 'peak'
            contrast_rates  = [contrast_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
        case 'area'
            contrast_rates  = [contrast_rates; mean(mean(experiment.whisk_rate(cond_inds == a,summarise_channels)))];
    end
    
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

%% Response size plot
if qrespsizeplot
    figure
    bar_handle  = bar([P_rate_LED_on P_rate_LED_off; A_rate_LED_on A_rate_LED_off],'LineWidth',2);
    legend(bar_handle,'LED on','LED off')
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    title('Response size, P vs. A whisker, LED ON vs. OFF')
    ylabel('Stimulus-evoked instantaneous firing rate')
    xlabel('1 = Principal, 2 = Adjacent')
    
    if save_figs
        print(gcf,[experiment_plot_folder filesep ' Bar graph P vs A - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
    end
end


% visual representation of responses across channels
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
