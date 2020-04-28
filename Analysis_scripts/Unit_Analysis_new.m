%% Unit_analysis_new
close all

%%
reload  = 1;
%%
if reload == 1
    clear all
load_folder = 'F:\Synced\Synced_individual'
group_folders     	= {'POM';'Cortical'}; 
reload =1;
end;


%
bin_size = 0.005;
spont_win = [0 0.05]
resp_win = [0.0001 0.03]; 
% Max and min depth for units to include
max_depth           = -50;  % Depth relative to L4 sink (negative numbers = towards pia)
min_depth           = -300; % Depth relative to L4 sink (negative numbers = towards pia)

L5_max_depth        =  400;   % Depth relative to L4 sink (negative numbers = towards pia)
L5_min_depth        =  100;  % Depth relative to L4 sink (negative numbers = towards pia)

L5_chans            = 5:12; % Layer 5 channels relative to the L4 sink (sink channel + [L5_chans])

%% Script starts : 
if reload == 1;
    for k = 1: numel(group_folders)
    load_fn = [load_folder '\' group_folders{k} '_selected_exps.mat']
    disp(['Loading...' load_fn]);
    load(load_fn);
    groups(k) = selected_exps;
    disp('Done');
  
    
    end;
    %% corrects unit depth to be same as that for multiunit channels 
    %% NB only do this once!
    
    for a = 1:numel(groups)
        for b = 1 : numel(groups(a).ephys_data)
            groups(a).ephys_data(b).unit_depths = 800-(groups(a).ephys_data(b).unit_depths)+25; 
        end;     
    end;

else 
    disp('No reload');    
end;

clear selected_exps;
a = 1;
b = 1;


%% Sanity check on unit correlations
% unit depth should positively correlate with depth of channel that
% activity is best related to. 
[unit_depths_out,chan_depths_out,corr_check] = unit_chan_corr_check(groups);
%%
for a = 1: numel(groups);
    Analysed_group(a).all_units = []
    for b = 1: numel(groups(a).ephys_data)
        if ~isempty(groups(a).ephys_data(b).conditions)
            Units = Prep_Unit_Analysis(groups(a).ephys_data(b),resp_win);
            Analysed(a).exp_preps(b).Units = Units;
            Analysed_group(a).all_units = [Analysed_group(a).all_units Units];
        end;
    end;
end;


%{
figure('Name','Spike_rate');

clear sig_delta;
    sig_delta(1:2,1:128,1:7) = NaN;
    sig_combined(1:2,1:128) = NaN;
     sig_delta_supergran(1:2,1:128,1:7) = NaN;
for a = 1:2;
    unit_count =numel(Analysed_group(a).all_units);
    cond_count = numel(Analysed_group(a).all_units(1).conditions);
 for b = 1: unit_count;
     for k = 1 : cond_count;
     subplot(2,1,a);
     title(group_folders{a});
         hold on;
         sig_index = Analysed_group(a).all_units(b).conditions(k).trial_rate_h;
         unit_depth = Analysed_group(a).all_units(b).unit_depth;
     hold on;
     refline(0,1);
     plot(k,Analysed_group(a).all_units(b).conditions(k).delta_whisk_rate,'bo');
     if sig_index == true;
     sig_delta(a,b,k) = Analysed_group(a).all_units(b).conditions(k).delta_whisk_rate;
       if unit_depth >400
           plot(k,Analysed_group(a).all_units(b).conditions(k).delta_whisk_rate,'r*');
       end;
         if unit_depth <400
            plot(k,Analysed_group(a).all_units(b).conditions(k).delta_whisk_rate,'g*');
            sig_delta_supergran(a,b,k) = Analysed_group(a).all_units(b).conditions(k).delta_whisk_rate;
         end;
     end;
     end;
     subplot(2,1,a);
     hold on;
     if Analysed_group(a).all_units(b).combined_rate_h == 1 & Analysed_group(a).all_units(b).unit_depth <400;
         plot(4.5,nanmean(Analysed_group(a).all_units(b).combined_response),'g*');
         sig_combined(a,b) = nanmean(Analysed_group(a).all_units(b).combined_response);
     end;
     
 end;
end;

figure('Name','Spike prob');
 sig_delta_prob(1:2,1:128,1:7) = NaN;
   
for a = 1:2;
    unit_count =numel(Analysed_group(a).all_units);
    cond_count = numel(Analysed_group(a).all_units(1).conditions);
 for b = 1: unit_count;
     for k = 1 : cond_count;
     subplot(2,1,a);
     title([group_folders{a} 'Probability']);
         hold on;
         sig_index = Analysed_group(a).all_units(b).conditions(k).prob_p;
         unit_depth = Analysed_group(a).all_units(b).unit_depth;
     hold on;
     refline(0,1);
     plot(k,Analysed_group(a).all_units(b).conditions(k).prob_delta,'bo');
     if sig_index == true;
     sig_delta_prob(a,b,k) = Analysed_group(a).all_units(b).conditions(k).prob_delta;
       if unit_depth >400
           plot(k,Analysed_group(a).all_units(b).conditions(k).prob_delta,'r*');
       end;
         if unit_depth <400
            plot(k,Analysed_group(a).all_units(b).conditions(k).prob_delta,'g*');
         end;
     end;
     end; 
 end;
end;

figure('Name','Spike latency');
   
for a = 1:2;
    unit_count =numel(Analysed_group(a).all_units);
    cond_count = numel(Analysed_group(a).all_units(1).conditions);
 for b = 1: unit_count;
     for k = 1 : cond_count;
     subplot(2,1,a);
     title([group_folders{a} 'Spike_latency']);
         hold on;
   refline(0,0);
         plot(k,Analysed_group(a).all_units(b).conditions(k).Spike_1_delta,'bo');
     end; 
 end;
end;


for k = 1 : 6
stat_h(k) = ranksum(sig_delta_supergran(1,:,k),sig_delta(2,:,k));
stat_prob_h(k) = ranksum(sig_delta_prob(1,:,k),sig_delta_prob(2,:,k));
end;

[h,p] = ranksum(sig_combined(1,:),sig_combined(2,:));
%}