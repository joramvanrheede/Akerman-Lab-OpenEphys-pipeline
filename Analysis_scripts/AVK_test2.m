close all;

psth_max_thresh = 5;
psth_mean_thresh = 1;
count =0;
psth_grid(1:7,1:20,1:50) = NaN;
psth_delta_grid(1:7,1:20,1:50) = NaN;
for a = 1 %: numel(Analysed_group)
    for b = 1 : numel(Analysed(a).exp_preps);
for k =  2 : numel(Analysed(a).exp_preps(b).Units);
    unit = Analysed(a).exp_preps(b).Units(k);
    opto_r = sum(unit.conditions(end).opto_behaviour.psth_100,1);
    whisk_r = sum(unit.conditions(end).whisk_behaviour.psth_100,1);
    opto_activity = (max(opto_r)  > psth_max_thresh) & (mean(opto_r)  > psth_mean_thresh);
    whisk_activity = (max(whisk_r) > psth_max_thresh) &(mean(whisk_r) > psth_mean_thresh);
    isactive(k) = opto_activity | whisk_activity;
   
    if isactive(k) & unit.unit_depth <300
    count = count+1;
    for j = 1:7
    psth_grid(j,:,count) =  sum(unit.conditions(j).whisk_behaviour.psth_100,1);
    cond_temp =sum(unit.conditions(j).whisk_behaviour.psth_100,1);
    control_temp =sum(unit.conditions(end).whisk_behaviour.psth_100,1);
    norm_temp = cond_temp-control_temp;
    psth_delta_grid(j,:,count) = norm_temp;
    end;
    close all;
    % quick_plot_unit(unit,'grid');
    % pause();
    end;
    
end;
end;
end;    
    
figure();
for k = 1:7
    subplot(1,7,k);
    hold on;
    disp(k);
    n = ones(20,1)*count;
    cond_psth_delta = squeeze(psth_delta_grid(k,:,:));
    mean_psth_delta(:,k) = nanmean(cond_psth_delta,2);
    std_psth_delta(:,k) = nanstd(cond_psth_delta,[],2);
    
    bar(mean_psth_delta(:,k))
    
    
end;

   


