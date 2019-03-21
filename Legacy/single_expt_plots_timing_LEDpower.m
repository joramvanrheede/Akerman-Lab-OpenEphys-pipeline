% single experiment plots for frequency experiment

close all

experiment          = sdata.expt(2);

resp_measure        = 'peak'; % 'peak' or 'area'
qraster             = true;
qpsth               = true;
qdelayvsresp        = true;
qLEDplot            = true;
qchannelim          = false;

save_figs           = false;
figure_dpi          = 300;


summarise_channels  = [13:16]; % include these channels
split_conditions    = [1 8]; % split by these conditions, summarise over others

split_plots         = [1 6 8]; % [4 6] works

% for rasterplots:
split_raster_plots  = [1];	% split plots by these conditions
split_figures       = [6]; 	% split 

% for raster plots:
trialrange          = [1 10]; % [min max] - don't exceed max nr of trials; else errors result.
x_ax_lims           = [0 4]; % limits for x-axes
condition_name    	= 'LED timing';
condition_units   	= 's';

save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';


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
    set(gcf,'Position',[0 .1 1 .9])
    for a = 1:size(split_cond_rows,1)
        
        these_conds         = split_cond_rows(a,:);
        sum_inds            = cond_inds == a;
        sum_inds            = find(sum_inds);
        for b = 1:length(sum_inds)
            this_whisk_psth  	= experiment.whisk_win_rates(sum_inds(b),:,:);
            
            subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
            plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),3),'LineWidth',2)
            hold on
        end
        
        xlim([-.5 1])
        set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
        hold on
        
    end
    
    set(gcf,'Color',[1 1 1])
    
    % Set all y axes to the same range (based on the largest range)
    plotaxes    = get(gcf,'Children');
    maxy        = cellfun(@max,get(plotaxes,'Ylim'));
    set(plotaxes,'Ylim',[0 max(maxy)]);
    
    % save figure
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'PSTHs chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    end
    
end





%
resp_conds      = [1 6];
resp_cond_mat   = condition_mat(:,resp_conds);
[resp_cond_rows, indxa, cond_inds] = unique(resp_cond_mat,'rows');
resp_rates      = [];
resp_times      = [];
for a  = unique(cond_inds)'
    
    switch resp_measure
        case 'area'
            resp_rates  = [resp_rates; mean(mean(experiment.whisk_rate(cond_inds == a,summarise_channels)))];
        case 'peak'
            resp_rates  = [resp_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
    end
    
    resp_times  = [resp_times; mean(mean(experiment.whisk_peak_time(cond_inds == a,summarise_channels)))];
end

if resp_rates(end-1) > resp_rates(end)
    P_whisk_stim = 1;
    A_whisk_stim = 2;
else
    P_whisk_stim = 2;
    A_whisk_stim = 1;
end

P_rates             = resp_rates(P_whisk_stim:2:end);
A_rates             = resp_rates(A_whisk_stim:2:end);

PA_ratios           = P_rates ./ A_rates;
LED_delays       	= unique(condition_mat(:,1) - condition_mat(:,2));

if qdelayvsresp
    figure
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.2 .4 .6 .4])
    
    subplot(1,2,1)
    % plots
    plot(LED_delays(1:end-1),P_rates(1:end-1),'k.-','LineWidth',4,'MarkerSize',30)
    hold on
    plot(LED_delays(1:end-1),A_rates(1:end-1),'r.-','LineWidth',4,'MarkerSize',30)
    
    xlim([(LED_delays(1) - .05) (LED_delays(end-1) + .05)])
    
    line(xlim,[P_rates(end) P_rates(end)],'LineWidth',4,'Color',[.5 .5 .5],'LineStyle','--');
    line(xlim,[A_rates(end) A_rates(end)],'LineWidth',4,'Color',[1 .5 .5],'LineStyle','--');
    
    title('Peak instantaneous rate of fire')
    
    legend({'Principal whisker' 'Adjacent whisker', 'Principal baseline', 'Adjacent baseline'},'Location','SouthEast')
    
    ylimits = ylim;
    ylim([0 ylimits(2)]);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16,'FontWeight','bold')
    set(gcf,'Color',[1 1 1])
    
    ylabel('Whisker response after LED stimulus at time X')
    xlabel('Time')
    title('Effect of timed LED pulse on whisker response')
    
    %%
    subplot(1,2,2)
    
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
    
    
    % save figure
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'Delay vs Resp chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    end
    
end

if qLEDplot
    % experiment.LED_win_counts
    
    LED_delays      = experiment.condition_mat(:,1);
    
    control_LED     = LED_delays == max(LED_delays);
    
    LED_rate_hist   = experiment.LED_win_rates(control_LED,summarise_channels,:);
    
    LED_rate_hist   = squeeze(mean(mean(LED_rate_hist,1),2));
    
    figure
    plot(experiment.LEDwinedges(1:end-1),smooth(LED_rate_hist,3),'LineWidth',2)
    set(gca,'FontName','Garamond','FontWeight','Bold','FontSize',16)
    set(gcf,'Color',[1 1 1])
    title('LED response')
    xlim([experiment.LEDwinedges(1) experiment.LEDwinedges(end)])
    xlabel('Time (s) relative to LED stimulus')
    ylabel('Firing rate')
    
    % save figure
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'LED plot ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    end
end

%% rasterplots
if qraster
    % use channels_to_rasterplots figure to generate
    fighandles = channels_to_rasterplots(experiment.filename,summarise_channels,split_raster_plots,split_figures,1:10,x_ax_lims,condition_name, condition_units);
    
    % save figure
    if save_figs
        for a = 1:length(fighandles)
            figure(fighandles(a))
            print(gcf,[experiment_plot_folder filesep 'Rasterplot stimulator ' num2str(a) ' chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
        end
    end
end


% visual representation of responses across channels
if qchannelim
    whisk_resps     = experiment.whisk_peak_rate;
    P_whisk_resps   = whisk_resps(P_whisk_stim:2:end,:);
    A_whisk_resps   = whisk_resps(A_whisk_stim:2:end,:);
    
    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.2 .2 .6 .6])
    
    subplot(1,3,1)
    imagesc(P_whisk_resps')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Principal whisker')
    ylabel('Channer number')
    colorbar
    % axis off
    
    subplot(1,3,2)
    imagesc(A_whisk_resps')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Adjacent whisker')
    colorbar
    % axis off
    
    subplot(1,3,3)
    imagesc(P_whisk_resps' ./ A_whisk_resps')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('P/A ratio')
    colorbar
    % axis off
    
    % save figure
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'Channel heatmaps'],'-dpng',['-r' num2str(figure_dpi)])
    end
end
