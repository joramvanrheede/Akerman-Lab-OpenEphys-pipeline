%% Unit_analysis_new
close all

%%
reload  = 1;
%%
if reload == 1
    clear all
load_folder = 'F:\Synced\Synced_individual'
group_folders     	= {'POM';'Cortical'}; 
target_opto = '60'
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
    load_fn = [load_folder '\' group_folders{k} '_selected_exps_' target_opto '.mat']
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




