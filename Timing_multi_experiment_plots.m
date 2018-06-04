% Multiple experiment plots for drive experiment

%% Add LED response size
% LED ON response peak
% LED ON sustained response (in control condition!)
% LED OFF response
% some responsiveness criterion implementation

%% Add examination of tracking repeat stimuli
% Plot response with error bars to first 4-5 stimuli
% Plot response to first stimulus vs response to stimulus 3 seconds later
% (by frequency)
% Make some metrics: first/fourth stim ratio, first vs 3 secs ratio
%


save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
save_expt_name      = 'Timing 20180417';
save_figs        	= false;

summarise_channels  = [13:16]; % include these channels
LEDresp_threshold   = 2.5;

q_check_LED_resp    = true;

% Zoomed in values for visualising instantaneous LED response
z_LED_tmin          = -.01;
z_LED_tmax          = .2;

%% Set for this type of experiment
split_conditions    = [1 5 6]; % split by these conditions, summarise over others
split_plots         = [1]; % [4 6] works

%%
close all

P_rate_LED_on       = [];
P_rate_LED_off      = [];
A_rate_LED_on       = [];
A_rate_LED_off      = [];
PA_ratio_LED_on     = [];
PA_ratio_LED_off    = [];
P_LED_onoff_ratio   = [];
A_LED_onoff_ratio   = [];
P_pktime_LED_on     = [];
P_pktime_LED_off	= [];
A_pktime_LED_on     = [];
A_pktime_LED_off	= [];
LED_rel_mean        = [];
LED_sust_rel_mean   = [];
LED_OFF_rel_mean    = [];
LED_rel_traces      = [];
LED_win_edges       = [];
counter             = 0;
LEDresp           	= [];
for i = 1:length(sdata)
    experiment          = sdata(i).expt;
    
    for j = 1:length(experiment)
        counter             = counter+1;
        
        condition_mat       = experiment(j).condition_mat;
        
        split_cond_mat      = condition_mat(:,split_conditions);
        [split_cond_rows, indxa, cond_inds] = unique(split_cond_mat,'rows');
        
        split_plot_mat      = condition_mat(:,split_plots);
        [split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');
        
        [a,b,cond_plot_inds] = unique(cond_plot_inds);
        
        %% Look at peak instantaneous firing rate between LED and no LED
        contrast_conds      = [1 6];
        contrast_cond_mat   = condition_mat(:,contrast_conds);
        [contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
        contrast_rates      = [];
        contrast_times      = [];
        contrast_delays     = [];
        contrast_stim       = [];
        for a  = unique(cond_inds)'
            
            these_peak_rates    = experiment(j).whisk_peak_rate(cond_inds == a,summarise_channels);
            these_spont_rates   = experiment(j).spont_rate(summarise_channels);
            these_rel_peaks     = these_peak_rates ./ these_spont_rates;
            
            contrast_rates      = [contrast_rates; mean(these_rel_peaks)];
            contrast_times      = [contrast_times; mean(mean(experiment(j).whisk_peak_time(cond_inds == a,summarise_channels)))];
            contrast_delays     = [contrast_delays; condition_mat(cond_inds == a,1) - condition_mat(cond_inds == a,2)];
            contrast_stim       = [contrast_stim; condition_mat(cond_inds == a,6)];
        end
        
        if contrast_rates(end - 1) > contrast_rates(end)
            P_whisk_stim = 1;
            A_whisk_stim = 2;
        else
            P_whisk_stim = 2;
            A_whisk_stim = 1;
        end
        
        P_rates(:,counter)      = contrast_rates([P_whisk_stim:2:end]);
        A_rates(:,counter)      = contrast_rates([A_whisk_stim:2:end]);
        
        P_times(:,counter)      = contrast_times([P_whisk_stim:2:end]);
        A_times(:,counter)      = contrast_times([A_whisk_stim:2:end]);
        
        P_delays(:,counter)     = contrast_delays([P_whisk_stim:2:end]);
        A_delays(:,counter)     = contrast_delays([A_whisk_stim:2:end]);
        
        
        PA_ratios(:,counter)    = contrast_rates([P_whisk_stim:2:end]) ./ contrast_rates([A_whisk_stim:2:end]);
        
        
        %% LED resp values
        
        LEDresp(counter)         	= mean(mean(experiment(j).LED_rel(:,[end-1 end])));
        
        %% LED PSTHs
        LED_rel_PSTHs               = experiment(j).LED_win_rel([end-1 end],summarise_channels,:);
        
        mean_LED_rel_PSTHs          = squeeze(mean(LED_rel_PSTHs,2));
        
        LED_rel_traces(:,counter) 	= mean(mean_LED_rel_PSTHs);
        LED_win_edges(:,counter)  	= experiment(j).LEDwinedges(1:end-1);
        
    end
end

qLEDresp            = LEDresp > LEDresp_threshold;


if q_check_LED_resp
	P_rates     = P_rates(:,qLEDresp);
	A_rates     = A_rates(:,qLEDresp);
        
	P_times     = P_times(:,qLEDresp);
	A_times     = A_times(:,qLEDresp);

	P_delays    = P_delays(:,qLEDresp);
	A_delays    = A_delays(:,qLEDresp);
        
	PA_ratios   = PA_ratios(:,qLEDresp);
    
end


%% Response size plotting

figure
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .4])
set(gcf,'Color',[1 1 1])
subplot(1,2,1)
plot(P_delays(1:end-1,:),P_rates(1:end-1,:),'k.-','LineWidth',2,'MarkerSize',25)
hold on
xlimits = xlim;
%line(xlimits,[P_rates(end), P_rates(end)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
%line(xlimits,[A_rates(end), A_rates(end)],'LineWidth',2,'LineStyle','--','Color',[1 0 0])
title('Principal whisker response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

subplot(1,2,2)
plot(A_delays(1:end-1,:),A_rates(1:end-1,:),'r.-','LineWidth',2,'MarkerSize',25)
hold on
title('Adjacent whisker response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])


if save_figs
    experiment_plot_folder      = [save_folder filesep save_expt_name];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
    print(gcf,[experiment_plot_folder filesep 'Principal vs Adjacent resp paired plots - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end


return

figure
set(gcf,'Color',[1 1 1])
pairedlineplot(P_LED_onoff_ratio,A_LED_onoff_ratio,{'Principal' 'Adjacent'},'Whisker','LED effect on whisker')
title('Effect of LED on P vs A whisker')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])


%% Peak time plotting

figure
set(gcf,'Units','normalized')
set(gcf,'Position',[.167 .4 .67 .4])
set(gcf,'Color',[1 1 1])

subplot(1,2,2)
pairedlineplot(A_pktime_LED_off,A_pktime_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Adjacent Whisker peak time')
title('Adjacent whisker peak time')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

subplot(1,2,1)
pairedlineplot(P_pktime_LED_off,P_pktime_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal Whisker peak time')
title('Principal whisker peak time')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylim([0 ylimits(2)*1.2])


%% LED plotting
figure
set(gcf,'Units','normalized')
set(gcf,'Position',[.167 .4 .67 .4])
set(gcf,'Color',[1 1 1])

subplot(1,2,1)
plot_handle = plot(LED_win_edges(:,~qLEDresp),LED_rel_traces(:,~qLEDresp),'r-','LineWidth',2);
hold on
plot_handle = plot(LED_win_edges(:,qLEDresp),LED_rel_traces(:,qLEDresp),'k-','LineWidth',2);
title('Average LED response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
xlabel('Time (s)')
ylabel('Spike rate relative to spontaneous')
axis tight
ylimits = ylim;
ylim([0 ylimits(2)*1.2])



qzoom           = LED_win_edges(:,1) >= z_LED_tmin & LED_win_edges(:,1) <= z_LED_tmax;
z_LED_tvals     = LED_win_edges(qzoom,1);
z_LED_traces    = LED_rel_traces(qzoom,:);

subplot(1,2,2)
plot(z_LED_tvals,z_LED_traces(:,~qLEDresp),'r-','LineWidth',2)
hold on
plot(z_LED_tvals,z_LED_traces(:,qLEDresp),'k-','LineWidth',2)
title('Average LED response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
xlabel('Time (s)')
ylabel('Spike rate relative to spontaneous')
axis tight
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

% subplot(1,2,2)
% pairedlineplot(LED_rel_mean,LED_sust_rel_mean,{'LED ON' 'LED SUSTAINED'},'Opto dynamics','Spike rate relative to mean')
% title('Adjacent whisker peak time')
% set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
% ylimits = ylim;
% ylim([0 ylimits(2)*1.2])


