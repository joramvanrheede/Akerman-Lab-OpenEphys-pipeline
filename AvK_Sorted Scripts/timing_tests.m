clc; clear all; close all;
reload = 1;
animal = '2019_11_28_2'
sorted_units = 'F:\Sorted_POM\Timing\2019_11_27_1\2019_11_27-3-Timing.mat'
multi_unit = 'F:\Multi_unit Coalated\1_POM\Timing\2019_11_27_1\2019_11_27-3-Timing.mat'

experiment_type = ['Timing'] % either 'Drive' or 'Timing';
Outputdirect = 'F:\Coalated_analysis\POM'
experiment = [animal '_timing_1'];
OutputFn = [Outputdirect experiment_type];
Layer_bounds_width = [178 241 314 273]; % L2/3 L4 L5 L6

%% script starts
if ~exist([OutputFn '\Figs'] , 'dir')
       mkdir([OutputFn '\Figs'])
end;

if (reload == 1)
    load(sorted_units);
    sorted_data = ephys_data;
    clear ephys_data;
    load(multi_unit);
    MUE_data =  ephys_data;
    clear ephys_data;
end;

LFP_out = LFP_Depth(experiment_type,MUE_data,OutputFn,true);
mean_sinks = 0 - mean(LFP_out.Sinks);
[~,a] = findpeaks(0-LFP_out.Sinks,'MinPeakHeight',mean_sinks);
[~,b] = min(LFP_out.LFP_Sinks);
[~,c] = min(abs(a-b));
L4_sink = a(c)+1; % Defines the sink in Layer 4 as the Peak CSD sink closes to the LFP sink. 

Channel_depths = L4_sink

clear a b c;

Layers.L4_sink_depth = 25*L4_sink;
Layers.L4 = [L4_sink-4:L4_sink+4];
Layers.L4_depth = [int16(Layers.L4_sink_depth-Layer_bounds_width(3)/2);int16(Layers.L4_sink_depth+Layer_bounds_width(3)/2)];
Layers.L5_depth = [Layers.L4_depth(2);Layers.L4_depth(2)+Layer_bounds_width(4)];
Layers.L2_3_depth = [Layers.L4_depth(1)-Layer_bounds_width(2);Layers.L4_depth(1)+50]; % gave 50um boundry to l2/3 for overlap with upper L4?
Layers.L5 = [Layers.L4(end)+1:Layers.L4(end)+12];
Layers.L2_3 = [Layers.L4(1)-7:Layers.L4(1)+1];% gave 50um boundry to l2/3 for overlap with upper L4?
Layers.L2_3 = Layers.L2_3(Layers.L2_3>0); % only includes up to highest channel
Layers.L5 = Layers.L5(Layers.L5<33); % only includes up to lowest channel

%% Whisk and Light Responsiveness

L2_3_units = sorted_data.unit_depths < Layers.L2_3_depth(2);


resp_win            = [0.005 0.055];
opto_resp_win       = [0.006 0.055];
control_win(1)      = 0-resp_win(2);
control_win(2)      = 0-resp_win(1);
opto_control_win(1) = 0 - opto_resp_win(2);
opto_control_win(2) = 0 - opto_resp_win(1);
this_t_opto         = sorted_data.conditions(end).LED_onset; 
this_t_whisk        = sorted_data.conditions(end).whisk_onset;


spikes          = sorted_data.conditions(end).spikes(:,:,:);
spikes          = spikes - this_t_whisk;
opto_spikes     = sorted_data.conditions(end).spikes(:,:,:) - this_t_opto;
delta_t         = this_t_opto - this_t_whisk;
n_trials        =  sorted_data.conditions(end).n_trials;  

L5_MUA_Spikes   = MUE_data.conditions(end).spikes(Layers.L5,:,:);
L5_MUA_Spikes   = L5_MUA_Spikes - this_t_opto;

Whisk_Resp = Stim_Responsive(spikes,resp_win,control_win,n_trials,delta_t,true,'Whisker Stim Response',Outputdirect); % tests whether spike rate or probablity is signifactly larger following stimulus
Opto_Resp = Stim_Responsive(opto_spikes,opto_resp_win,opto_control_win,n_trials,delta_t,true,'Opto Stim Response',Outputdirect);

L5_opto_recruitment = spike_rate_in_win(L5_MUA_Spikes,opto_resp_win);
L5_opto_control = spike_rate_in_win(L5_MUA_Spikes,opto_control_win);
L5_delta_recruitment = L5_opto_recruitment - L5_opto_control;

%%
window = [0.004 0.1]

unit_of_interest = squeeze(sorted_data.conditions);

for k = 1 : numel(unit_of_interest)
    unit_of_interest(k).spikes = squeeze(unit_of_interest(k).spikes(1,:,:));
    unit_of_interest(k).spikes = unit_of_interest(k).spikes - this_t_whisk; 
end;

clear k; 
[unit_of_interest(end).first_spike,unit_of_interest(end).second_spike] = find_spike_times(unit_of_interest(end).spikes,window);

unit_of_interest(end).first_spike_prob = sum(~isnan(unit_of_interest(end).first_spike))/numel(unit_of_interest(end).first_spike);
unit_of_interest(end).second_spike_prob = sum(~isnan(unit_of_interest(end).second_spike))/numel(unit_of_interest(end).second_spike);
unit_of_interest(end).second_spike_prob = unit_of_interest(end).first_spike_prob*unit_of_interest(end).second_spike_prob;

for k = 1 : numel(unit_of_interest)-1;
[unit_of_interest(k).first_spike,unit_of_interest(k).second_spike] = find_spike_times(unit_of_interest(k).spikes,window);
unit_of_interest(k).first_spike_prob = sum(~isnan(unit_of_interest(k).first_spike))/numel(unit_of_interest(k).first_spike);
unit_of_interest(k).second_spike_prob = sum(~isnan(unit_of_interest(k).second_spike))/numel(unit_of_interest(k).second_spike);
unit_of_interest(k).second_spike_prob = unit_of_interest(k).first_spike_prob*unit_of_interest(k).second_spike_prob;

spike1_test = unit_of_interest(k).first_spike';
spike1_control = unit_of_interest(end).first_spike';

[n,i] = max([numel(spike1_test) numel(spike1_control)]);
if i == 1
    pad = nan*(numel(spike1_control)+1:n)
    spike1_control(numel(spike1_control)+1:n) = pad;
else
    pad = nan*(numel(spike1_test)+1:n)
    spike1_test(numel(spike1_test)+1:n) = pad;   
end;

spike2_test = unit_of_interest(k).second_spike';
spike2_control = unit_of_interest(end).second_spike';

[n,i] = max([numel(spike2_test) numel(spike2_control)]);
if i == 1
    pad = NaN*(numel(spike2_control)+1:n)
    spike2_control(numel(spike2_control)+1:n) = pad;
else
    pad = NaN*(numel(spike2_test)+1:n)
    spike2_test(numel(spike2_test)+1:n) = pad;   
end;


[unit_of_interest(k).first_spike_diff_p,unit_of_interest(k).first_spike_diff_h] = ranksum(unit_of_interest(k).first_spike,unit_of_interest(end).first_spike);
[unit_of_interest(k).second_spike_diff_p,unit_of_interest(k).second_spike_diff_h] = ranksum(unit_of_interest(k).second_spike,unit_of_interest(end).second_spike);
[unit_of_interest(k).first_spike_var_p,unit_of_interest(k).first_spike_var_stats] = vartestn([spike1_test spike1_control]);
[unit_of_interest(k).second_spike_var_p,unit_of_interest(k).second_spike_var_stats] = vartestn([spike2_test spike2_control]);

clear n i pad spike1_control spike1_test;

 [unit_of_interest(k).spike1_prob_diff_p,unit_of_interest(k).spike1_prob_diff_h] = Chi2_test(unit_of_interest(k).first_spike_prob,numel(unit_of_interest(k).first_spike),unit_of_interest(end).first_spike_prob,numel(unit_of_interest(end).first_spike),0.05);
 [unit_of_interest(k).spike2_prob_diff_p,unit_of_interest(k).spike2_prob_diff_h] = Chi2_test(unit_of_interest(k).second_spike_prob,numel(unit_of_interest(k).first_spike),unit_of_interest(end).second_spike_prob,numel(unit_of_interest(end).first_spike),0.05);


end;
clear k;

close all;

 ha = figure('Name','Spike Timing','NumberTitle','off');
    set(gcf,'Units','normalized','Position',[.2 .1 .4 .8])% left bottom width height
   
subplot(4,1,1);
for k = 1: numel(unit_of_interest);
hold on;
    for j = 1: numel(unit_of_interest(k).first_spike)
        plot(k,unit_of_interest(k).first_spike(j),'k+','MarkerSize',1);
    end;
plot(k,nanmean(unit_of_interest(k).first_spike),'ko')
    if ~isempty(unit_of_interest(k).first_spike_diff_h);
        plot(k,unit_of_interest(k).first_spike_diff_h*window(2),'k*')
        if(unit_of_interest(k).first_spike_var_p <0.05)
            plot(k,window(2)-0.1*window(2),'r*');
        end;
    end;
    
end;
xlim([0 numel(unit_of_interest)+1]);
ylabel('Time (s)');
xlabel('Delay condition');
title('Time to first spike');

subplot(4,1,2);
for k = 1: numel(unit_of_interest);
hold on;
    for j = 1: numel(unit_of_interest(k).second_spike)
        plot(k,unit_of_interest(k).second_spike(j),'k+','MarkerSize',1);
    end;

    plot(k,nanmean(unit_of_interest(k).second_spike),'ko')
    if ~isempty(unit_of_interest(k).second_spike_diff_h);
        plot(k,unit_of_interest(k).second_spike_diff_h*0.5,'k*');
        if(unit_of_interest(k).second_spike_var_p <0.05)
            plot(k,0.4,'r*');
        end;
   
    end;
    
end;
xlim([0 numel(unit_of_interest)+1]);
ylabel('Time (s)');
xlabel('Delay condition');
title('Time to second spike');

subplot(4,1,3);
for k = 1: numel(unit_of_interest);
   hold on;
   plot(k,unit_of_interest(k).first_spike_prob,'ko');      
end;
xlim([0 numel(unit_of_interest)+1]);
xlabel('Delay condition');
ylabel('Time (s)');
title('Probability of first spike');



subplot(4,1,4)
for k = 1: numel(unit_of_interest);
hold on;
        plot(k,unit_of_interest(k).second_spike_prob,'ko');
end;
xlim([0 numel(unit_of_interest)+1]);
xlabel('Delay condition');
ylabel('Time (s)');
title('Time to first spike');
title('Probability of second spike');
clear k j;

%{

figure();
hold on;
plot(unit_of_interest(end).first_spike,'b*');
plot(unit_of_interest(end).second_spike,'r*');



for k = 4:7;
figure();
subplot(3,1,1);
raster_plot(sorted_data.conditions(k).spikes(2,:,:),2,0.2);
title(['Trial : ' num2str(k) ]);
xlim([0 4]);

subplot(3,1,2);
raster_plot(MUE_data.conditions(k).spikes(Layers.L5,:,:),2,0.2);
xlim([0 4]);

subplot(3,1,3)
LFP = MUE_data.conditions(k).LFP_trace(L2_3_units,:,:);
LFP = squeeze(nanmean(LFP,1));
LFP = nanmean(LFP,1);
plot(LFP);
xlim([0 4000]);

end;

%}
