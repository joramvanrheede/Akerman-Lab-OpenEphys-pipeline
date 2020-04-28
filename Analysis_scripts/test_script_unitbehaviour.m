
clear opto_effect opto_effect_pom L2_3 unit med neg;
hold on;
%close all;
clear opto_effect opto_effect_pom L2_3 unit med neg; 

all_times(1:2,1:7,1:200) = NaN;
sig_times(1:2,1:7,1:200) = NaN;
      
for b = 1: 2
    figure('Name',['Group : ' num2str(b)]);
for a = 1: 6
    delta_cum_psth(1:numel(Analysed_group(b).all_units),1:20) = NaN;
for k = 1 : numel(Analysed_group(b).all_units)
    unit = Analysed_group(b).all_units(k);
if (unit.whisk_responsive & unit.unit_depth <400)
    trial_spike_times = unit.conditions(a).whisk_behaviour.second_spikes_by_trial ;
    control_spike_times = unit.conditions(end).whisk_behaviour.second_spikes_by_trial;
     delta_spike_times = nanmean((trial_spike_times)-(control_spike_times));
     all_times(b,a,k) = delta_spike_times;
     
     h =0 
     %[h,p] = ranksum(trial_spike_times,control_spike_times);
     hold on;
     if h ==1;
         plot(a,delta_spike_times,'g*');
          sig_times(b,a,k) = delta_spike_times;
        else
         plot(a,delta_spike_times,'ro');
        
     end;
end;
end;
end;
end;

pom_jitter = squeeze(all_times(1,:,:));
cortical_jitter = squeeze(all_times(2,:,:));

for k = 1 : 6
    [p(k),h(k)] = ranksum(pom_jitter(k,:),cortical_jitter(k,:));
end;


