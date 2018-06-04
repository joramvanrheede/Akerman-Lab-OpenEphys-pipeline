% Multiple experiment plots for drive experiment



save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
save_figs        	= false;

summarise_channels  = [1:16]; % include these channels

%% Set for this type of experiment
split_conditions    = [1 5 6]; % split by these conditions, summarise over others
split_plots         = [5 6]; % [4 6] works

%%
close all

P_rate_LED_on       = NaN(length(sdata),1);
P_rate_LED_off      = NaN(length(sdata),1);
A_rate_LED_on       = NaN(length(sdata),1);
A_rate_LED_off      = NaN(length(sdata),1);
PA_ratio_LED_on     = NaN(length(sdata),1);
PA_ratio_LED_off    = NaN(length(sdata),1);
P_LED_onoff_ratio   = NaN(length(sdata),1);
A_LED_onoff_ratio   = NaN(length(sdata),1);
for i = 1:length(sdata)
    experiment          = sdata(i).expt;
    
    %%
    condition_mat       = experiment.condition_mat;
    
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
        contrast_rates      = [contrast_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
        contrast_times      = [contrast_times; mean(mean(experiment.whisk_peak_time(cond_inds == a,summarise_channels)))];
    end
    
    if contrast_rates(3) > contrast_rates(4)
        P_whisk_stim = 1;
        A_whisk_stim = 2;
    else
        P_whisk_stim = 2;
        A_whisk_stim = 1;
    end
    
    P_rate_LED_on(i)        = contrast_rates(P_whisk_stim);
    P_rate_LED_off(i)       = contrast_rates(P_whisk_stim+2);
    A_rate_LED_on(i)        = contrast_rates(A_whisk_stim);
    A_rate_LED_off(i)       = contrast_rates(A_whisk_stim+2);
    
    PA_ratio_LED_on(i)      = P_rate_LED_on(i) / A_rate_LED_on(i);
    PA_ratio_LED_off(i)     = P_rate_LED_off(i) / A_rate_LED_off(i);
    
    P_LED_onoff_ratio(i)    = P_rate_LED_on(i) / P_rate_LED_off(i);
    A_LED_onoff_ratio(i)    = A_rate_LED_on(i) / A_rate_LED_off(i);
    
end

P_rate_LED_on
figure
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .6 1 .4])
subplot(1,3,1)
pairedlineplot(P_rate_LED_off,P_rate_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal whisker firing rate')
title('Principal whisker firing rate')
subplot(1,3,2)
pairedlineplot(A_rate_LED_off,A_rate_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Adjacent whisker firing rate')
title('Adjacent whisker firing rate')
subplot(1,3,3)
pairedlineplot(PA_ratio_LED_off,PA_ratio_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal/Adjacent ratio')
title('Principal/Adjacent ratio')

figure
pairedlineplot(P_LED_onoff_ratio,A_LED_onoff_ratio,{'Principal' 'Adjacent'},'Whisker','LED ON / OFF ratio')
title('LED ON/OFF ratio')


