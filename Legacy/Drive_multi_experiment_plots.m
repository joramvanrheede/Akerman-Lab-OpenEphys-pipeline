% Multiple experiment plots for drive experiment


save_figs        		= false;
    save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
    save_expt_name      = 'Drive 20180626';

qpairedlineplots    = true;
resp_measure        = 'area'; % 'area' / 'peak'

qLEDplots           = true;
    z_LED_tmin      = -0.05;
    z_LED_tmax      = 0.15;
    
summarise_channels  = [1:16]; % include these channels

LEDresp_threshold   = 50;
PAratio_threshold   = 1.1;

%% Set for this type of experiment
split_conditions    = [1 5 6]; % split by these conditions, summarise over others
split_plots         = [5 6]; % [4 6] works

%%

if save_figs
    experiment_plot_folder      = [save_folder filesep save_expt_name];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
end


%%
close all

LEDresp             = [];
LED_traces          = [];
LED_win_edges       = [];
P_rate_LED_on       = []; % NaN(length(sdata),1);
P_rate_LED_off      = []; % NaN(length(sdata),1);
A_rate_LED_on       = []; % NaN(length(sdata),1);
A_rate_LED_off      = []; % NaN(length(sdata),1);
PA_ratio_LED_on     = []; % NaN(length(sdata),1);
PA_ratio_LED_off    = []; % NaN(length(sdata),1);
P_LED_onoff_ratio   = []; % NaN(length(sdata),1);
A_LED_onoff_ratio   = []; % NaN(length(sdata),1);
counter             = 0;
for i =  1:length(sdata)
    experiment          = sdata(i).expt;
    
    for j = 1:length(experiment)
        counter     = counter + 1;
        %%
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
        end
        
        if contrast_rates(3) > contrast_rates(4)
            P_whisk_stim = 1;
            A_whisk_stim = 2;
        else
            P_whisk_stim = 2;
            A_whisk_stim = 1;
        end
        
        P_rate_LED_on(counter)     	= contrast_rates(P_whisk_stim);
        P_rate_LED_off(counter)  	= contrast_rates(P_whisk_stim+2);
        A_rate_LED_on(counter)    	= contrast_rates(A_whisk_stim);
        A_rate_LED_off(counter)   	= contrast_rates(A_whisk_stim+2);
        
        PA_ratio_LED_on(counter)  	= P_rate_LED_on(counter) / A_rate_LED_on(counter);
        PA_ratio_LED_off(counter) 	= P_rate_LED_off(counter) / A_rate_LED_off(counter);
        
        P_LED_onoff_ratio(counter)	= P_rate_LED_on(counter) / P_rate_LED_off(counter);
        A_LED_onoff_ratio(counter)	= A_rate_LED_on(counter) / A_rate_LED_off(counter);
        
        %% LED resp values
        
        LED_traces                  = [LED_traces squeeze(mean(mean(experiment(j).LED_win_rates(3:4,:,:))))];
        LED_win_edges               = [LED_win_edges experiment(j).LEDwinedges(:)];
        
        LEDresp(counter)         	= mean(mean(experiment(j).LED_rate(:,[end-1 end])));
    end
end


qLEDresp            = LEDresp > LEDresp_threshold;
qPAratio            = PA_ratio_LED_off > PAratio_threshold;

P_rate_LED_off      = P_rate_LED_off(qLEDresp & qPAratio);
P_rate_LED_on       = P_rate_LED_on(qLEDresp & qPAratio);
A_rate_LED_off      = A_rate_LED_off(qLEDresp & qPAratio);
A_rate_LED_on       = A_rate_LED_on(qLEDresp & qPAratio);
PA_ratio_LED_off    = PA_ratio_LED_off(qLEDresp & qPAratio);
PA_ratio_LED_on     = PA_ratio_LED_on(qLEDresp & qPAratio);
P_LED_onoff_ratio   = P_LED_onoff_ratio(qLEDresp & qPAratio);
A_LED_onoff_ratio   = A_LED_onoff_ratio(qLEDresp & qPAratio);


if qpairedlineplots
    figure
    set(gcf,'Units','normalized')
    set(gcf,'Position',[0 .6 1 .4])
    subplot(1,3,1)
    pairedlineplot(P_rate_LED_off,P_rate_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal whisker firing rate')
    title('Principal whisker firing rate')
    ylimz = ylim;
    ylim([0 ylimz(2)])
    
    subplot(1,3,2)
    pairedlineplot(A_rate_LED_off,A_rate_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Adjacent whisker firing rate')
    title('Adjacent whisker firing rate')
    ylim([0 ylimz(2)])
    
    subplot(1,3,3)
    pairedlineplot(PA_ratio_LED_off,PA_ratio_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal/Adjacent ratio')
    title('Principal/Adjacent ratio')
    
    % save figure
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'Pairedlineplots chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    end
    
    figure
    pairedlineplot(P_LED_onoff_ratio,A_LED_onoff_ratio,{'Principal' 'Adjacent'},'Whisker','LED ON / OFF ratio')
    title('LED ON/OFF ratio')
    
end



if qLEDplots
    figure
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.167 .4 .67 .4])
    set(gcf,'Color',[1 1 1])
    
    subplot(1,2,1)
    if sum(~qLEDresp) ~= 0
        plot_handle = plot(LED_win_edges(1:end-1,~qLEDresp),LED_traces(:,~qLEDresp),'r-','LineWidth',2);
    end
    hold on
    if sum(qLEDresp) ~= 0
        plot_handle = plot(LED_win_edges(1:end-1,qLEDresp),LED_traces(:,qLEDresp),'k-','LineWidth',2);
    end
    title('Average LED response')
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    xlabel('Time (s)')
    ylabel('Spike rate')
    axis tight
    ylimits = ylim;
    ylim([0 ylimits(2)*1.2])
    
    
    qzoom           = LED_win_edges(:,1) >= z_LED_tmin & LED_win_edges(:,1) <= z_LED_tmax;
    qzoom           = qzoom(1:end-1);
    z_LED_tvals     = LED_win_edges(qzoom,1);
    z_LED_traces    = LED_traces(qzoom,:);
    
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
    ylabel('Spike rate')
    axis tight
    ylimits = ylim;
    ylim([0 ylimits(2)*1.2])
    
    % save figure
    if save_figs
        print(gcf,[experiment_plot_folder filesep 'LED plots chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    end
    
end
