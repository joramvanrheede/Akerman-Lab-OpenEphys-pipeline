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


summarise_channels      = [1:16]; % include these channels

resp_measure            = 'area'; % area or peak

save_figs        		= false;
    save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
    save_expt_name      = 'Timing 20180627';

qpktimeplots            = false;

qtimingplots            = true;
    q_norm              = true;

q_check_resp            = true;
    LEDresp_threshold   = 100;
    PAratio_theshold    = 1.3;

% Zoomed in values for visualising instantaneous LED response
q_LED_plots             = true;
    z_LED_tmin          = -.01;
    z_LED_tmax          = .05;

%% Set for this type of experiment
split_conditions    = [1 5 6]; % split by these conditions, summarise over others
split_plots         = [1]; % [4 6] works

%%
close all

P_rates             = [];
P_times             = [];
P_delays            = [];
A_rates             = [];
A_times             = [];
A_delays            = [];
PA_ratios           = [];
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
for i = [2 5 6] % 1:length(sdata)
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
            %%
            switch resp_measure
                case 'peak'
                    these_peak_rates    = experiment(j).whisk_peak_rate(cond_inds == a,summarise_channels);
                case 'area'
                    these_peak_rates    = experiment(j).whisk_rate(cond_inds == a,summarise_channels);
            end
            these_spont_rates   = experiment(j).spont_rate(summarise_channels);
            these_rel_peaks     = these_peak_rates - these_spont_rates;
            
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
        
        LEDresp(counter)         	= mean(mean(experiment(j).LED_rate(:,[end-1 end])));
        
        %% LED PSTHs
        LED_rel_PSTHs               = experiment(j).LED_win_rel([end-1 end],summarise_channels,:);
        
        mean_LED_rel_PSTHs          = squeeze(mean(LED_rel_PSTHs,2));
        
        LED_rel_traces(:,counter) 	= mean(mean_LED_rel_PSTHs);
        LED_win_edges(:,counter)  	= experiment(j).LEDwinedges(1:end-1);
        
    end
end


qLEDresp            = LEDresp > LEDresp_threshold;
qPAratio            = PA_ratios(end,:) > PAratio_theshold;

if q_check_resp
	P_rates     = P_rates(:,qLEDresp & qPAratio);
	A_rates     = A_rates(:,qLEDresp & qPAratio);
        
	P_times     = P_times(:,qLEDresp & qPAratio);
	A_times     = A_times(:,qLEDresp & qPAratio);

	P_delays    = P_delays(:,qLEDresp & qPAratio);
	A_delays    = A_delays(:,qLEDresp & qPAratio);
        
	PA_ratios   = PA_ratios(:,qLEDresp & qPAratio);
    
end


%% Response size plotting
if qtimingplots
    
    if q_norm
        % PW response divided by control condition
        norm_P_rates    = P_rates(1:end-1,:)./repmat(P_rates(end,:),size(P_rates,1)-1,1);
        norm_P_delays   = P_delays(1:end-1,1);
        
        % AW response divided by control condition
        norm_A_rates    = A_rates(1:end-1,:)./repmat(A_rates(end,:),size(A_rates,1)-1,1);
        norm_A_delays   = A_delays(1:end-1,1);
    else
        % PW response divided by control condition
        norm_P_rates    = P_rates(1:end-1,:);
        norm_P_delays   = P_delays(1:end-1,1);
        
        % AW response divided by control condition
        norm_A_rates    = A_rates(1:end-1,:);
        norm_A_delays   = A_delays(1:end-1,1);
    end
    
    disp('LED effects for the different delays on P whisker:')
    [hP,pP] = ttest2(P_rates(1:end-1,:)',repmat(P_rates(end,:),size(P_rates,1)-1,1)')
    
    disp('LED effects for the different delays on A whisker:')
    [hA,pA] = ttest2(A_rates(1:end-1,:)',repmat(A_rates(end,:),size(A_rates,1)-1,1)')
    
    % Comparison of LED effect on P vs A over the different delays
    disp('Comparison of LED effect on P vs A:')
    [hPA,pPA] = ttest2(norm_A_rates', norm_P_rates')
    
    % Figure setup
    figure
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.2 .4 .6 .4])
    set(gcf,'Color',[1 1 1])
    
    % Plot individual lines for PW
    subplot(1,2,1)
    plot(norm_P_delays,norm_P_rates,'Color',[0 0 0 .2],'LineWidth',4,'MarkerSize',25)
    hold on
    
    % Errorbar plot
    errorbar(norm_P_delays,mean(norm_P_rates,2),serr(norm_P_rates'),'LineWidth',2,'Color',[0 0 0])
    
    % More aesthetics
    xlimits = xlim;
    line(xlimits,[1 1],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
    title('Principal whisker response')
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
    ylimits = ylim;
    ylim([0 ylimits(2)*1.2])
    
    % Plot individual lines for AW
    subplot(1,2,2)
    plot(A_delays(1:end-1,:),norm_A_rates,'Color',[1 0 0 .2],'LineWidth',4,'MarkerSize',25)
    hold on
    
    % Errorbar plot
    errorbar(norm_A_delays,mean(norm_A_rates,2),serr(norm_A_rates'),'LineWidth',2,'Color',[1 0 0])
    
    
    xlimits = xlim;
    line(xlimits,[1 1],'LineWidth',2,'LineStyle','--','Color',[1 0 0])
    title('Adjacent whisker response')
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
    ylimits = ylim;
    ylim([0 ylimits(2)*1.2])
    
    if save_figs
        experiment_plot_folder      = [save_folder filesep save_expt_name];
        if ~isdir(experiment_plot_folder)
            mkdir(experiment_plot_folder)
        end
        print(gcf,[experiment_plot_folder filesep 'Timing plots PW and AW - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
    end
end



%% Peak time plotting
if qpktimeplots
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
end


%% LED plotting

if q_LED_plots
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
    if sum(~qLEDresp) ~= 0
        plot(z_LED_tvals,z_LED_traces(:,~qLEDresp),'r-','LineWidth',2)
    end
    hold on
    if sum(qLEDresp) ~= 0
        plot(z_LED_tvals,z_LED_traces(:,qLEDresp),'k-','LineWidth',2)
    end
    title('Average LED response')
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    xlabel('Time (s)')
    ylabel('Spike rate relative to spontaneous')
    axis tight
    ylimits = ylim;
    ylim([0 ylimits(2)*1.2])
    
end
% subplot(1,2,2)
% pairedlineplot(LED_rel_mean,LED_sust_rel_mean,{'LED ON' 'LED SUSTAINED'},'Opto dynamics','Spike rate relative to mean')
% title('Adjacent whisker peak time')
% set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
% ylimits = ylim;
% ylim([0 ylimits(2)*1.2])


