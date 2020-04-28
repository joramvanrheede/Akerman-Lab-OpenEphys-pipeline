function [Unit_behaviour] = all_unit_behavior(groups,spont_win,resp_win,target_dt)


all_unit_isi = [];
all_unit_rate = [];
all_unit_depth = [];
all_unit_opto_rate =[];
all_unit_opto_isi = [];
all_unit_whisk_rate = [];
all_unit_whisk_isi = [];


for b = 1 : numel(groups);
for a = 1 : numel(groups(b).ephys_data)
    if ~isempty(groups(b).ephys_data(a).conditions)
        all_unit_depth = [all_unit_depth;groups(b).ephys_data(a).unit_depths];
        clear mean_isi mean_spike_rate;
        for j = 1 : numel(groups(b).ephys_data(a).conditions) % for all conditions 
           spikes = groups(b).ephys_data(a).conditions(j).unit_spikes;
                
                    clear spont_trial_rate opto_trial_rate whisk_trial_rate;
                    spont_trial_rate(1:size(spikes,1),1:size(spikes,2)) = NaN; % creates an empty NaN array size no Units x no_trials
                    opto_trial_rate(1:size(spikes,1),1:size(spikes,2)) = NaN;
                    whisk_trial_rate(1:size(spikes,1),1:size(spikes,2)) = NaN;
                    
                for k = 1:size(spikes,1); % for all units
                    
                    unit_spikes = squeeze(spikes(k,:,:)); % gets spikes for single unit(k) in given condition (j).
                    [mean_isi(k,j),mean_spike_rate(k,j),s_r] =  Unit_behaviour_analysis(unit_spikes,spont_win); %analysises unit response behavour
                     spont_trial_rate(k,1:numel(s_r)) = squeeze(s_r); clear s_r;
                                         
                     if j == target_dt;%numel(groups(b).ephys_data(a).conditions) % if control condition
                         Unit_Opto_spikes = unit_spikes - groups(b).ephys_data(a).conditions(j).LED_onset; % gets spikes in opto window
                         Unit_Whisk_spikes = unit_spikes - groups(b).ephys_data(a).conditions(j).whisk_onset; % gets spikes in whisk window
                 
                       [opto_isi(k),opto_spike_rate(k),o_r] =  Unit_behaviour_analysis(Unit_Opto_spikes,resp_win);
                      [whisk_isi(k),whisk_spike_rate(k),w_r] = Unit_behaviour_analysis(Unit_Whisk_spikes,resp_win);    
                      opto_trial_rate(k,1:numel(o_r)) = squeeze(o_r); clear o_r;
                      whisk_trial_rate(k,1:numel(w_r)) = squeeze(w_r); clear w_r;
                      
                      Unit_behaviour.groups(b).ephys_data(a).Spont = spont_trial_rate; % spont_trial_rate n(no units)x m(no trials) matrix for given condition 
                      Unit_behaviour.groups(b).ephys_data(a).Opto = opto_trial_rate;
                      Unit_behaviour.groups(b).ephys_data(a).Whisk = whisk_trial_rate;
                     end; % ends if control condition
           
                 end; % ends looping through each unit
              
        end; % ends looping through all conditions
         
         mean_isi = nanmean(mean_isi,2);
         mean_spike_rate = nanmean(mean_spike_rate,2);
         all_unit_isi = [all_unit_isi;mean_isi];
        all_unit_rate = [all_unit_rate;mean_spike_rate];
        
        all_unit_opto_rate =[all_unit_opto_rate opto_spike_rate];
        all_unit_opto_isi = [all_unit_opto_isi opto_isi];
        all_unit_whisk_rate = [all_unit_whisk_rate whisk_spike_rate];
        all_unit_whisk_isi = [all_unit_whisk_isi whisk_isi];
    end;    
    clear mean_isi mean_spike_rate opto_isi opto_spike_rate whisk_isi whisk_spike_rate;
end; % ends looping through all experimental preps
end; % ends looping through each group

clear mean_unit_rate std_unit_rate rate_count mean_unit_isi std_unit_isi max_unit_rate;

Unit_behaviour.all_unit_depth = 800-all_unit_depth; %% NB takes unit depth and inverts it from kilosort output.

Unit_behaviour.all_unit_rate = all_unit_rate;
Unit_behaviour.all_unit_isi = all_unit_isi;

Unit_behaviour.all_unit_opto_rate = all_unit_opto_rate;
Unit_behaviour.all_unit_opto_isi = all_unit_opto_isi;

Unit_behaviour.all_unit_whisk_rate = all_unit_whisk_rate;
Unit_behaviour.all_unit_whisk_isi = all_unit_whisk_isi;


opto_delta_rate = all_unit_opto_rate - all_unit_rate';
whisk_delta_rate = all_unit_whisk_rate - all_unit_rate';

n_bins = 16
bin_edges = 0:50:50*n_bins;
[~,~,depth_bins] = histcounts(all_unit_depth,bin_edges);

%output_depths
Unit_behaviour.depth_bins = depth_bins; 
Unit_behaviour.depth_edges = bin_edges;

for k = 1 : n_bins  %analysis by depth
    iq = depth_bins == k; % for each depth bin ( n = around 16)
    
    %% mean rates
    Unit_behaviour.mean_unit_rate(k) = nanmean(all_unit_rate(iq)); % mean unit rate for a bin
    Unit_behaviour.mean_opto_rate(k) = nanmean(all_unit_opto_rate(iq));
    Unit_behaviour.mean_opto_delta(k) = nanmean(opto_delta_rate(iq));
    Unit_behaviour.mean_whisk_rate(k) = nanmean(all_unit_whisk_rate(iq));
    Unit_behaviour.mean_whisk_delta(k) = nanmean(whisk_delta_rate(iq));
    
    % max rates
    z = nanmax(all_unit_rate(iq));
     if ~isempty(z)
         Unit_behaviour.max_unit_rate(k) = z; 
     else
         Unit_behaviour.max_unit_rate(k) = 0;
     end;
     clear z;
     z = nanmax(all_unit_opto_rate(iq));
     if ~isempty(z)
         Unit_behaviour.max_opto_rate(k) = z; 
     else
         Unit_behaviour.max_opto_rate(k) = 0;
     end;
    
     clear z;
     z = nanmax(all_unit_whisk_rate(iq));
     if ~isempty(z)
         Unit_behaviour.max_whisk_rate(k) = z; 
     else
         Unit_behaviour.max_whisk_rate(k) = 0;
     end;
     
     
     %std on rates
    Unit_behaviour.std_unit_rate(k) = nanstd(all_unit_rate(iq));
    Unit_behaviour.std_unit_opto_rate(k) = nanstd(all_unit_opto_rate(iq));
    Unit_behaviour.std_opto_delta(k) = nanstd(opto_delta_rate(iq));
    Unit_behaviour.std_unit_whisk_rate(k) = nanstd(all_unit_whisk_rate(iq));
    Unit_behaviour.std_whisk_delta(k) = nanstd(whisk_delta_rate(iq));
    
    %n for each depth bin
    Unit_behaviour.rate_count(k)       = numel(all_unit_rate(iq));
    Unit_behaviour.opto_rate_count(k)  = numel(all_unit_opto_rate(iq));
    Unit_behaviour.opto_delta_count(k)      = numel(opto_delta_rate(iq));
    Unit_behaviour.whisk_rate_count(k)  = numel(all_unit_whisk_rate(iq));
    Unit_behaviour.whisk_delta_count(k)      = numel(whisk_delta_rate(iq));
    
    % isi
    Unit_behaviour.mean_unit_isi(k) = nanmean(all_unit_isi(iq));
    Unit_behaviour.std_unit_isi(k) = nanstd(all_unit_isi(iq));
    
end
