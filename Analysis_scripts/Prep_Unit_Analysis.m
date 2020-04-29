function units = Prep_Unit_Analysis(ephys_data,resp_win);


% whisk response threshold
whisk_resp_p_thresh = 0.2;
opto_resp_p_thresh = 0.05;
whisk_resp_rate_thresh = 10;
opto_resp_rate_thresh = 10;

L5_chans = ephys_data.LFP_min_chan+4:ephys_data.LFP_min_chan+16; % mid layer 4+100um : mid layer 4 + 400um assumes layer width of +/- 300um for layer 5;
L5_chans(L5_chans>32) = [];
L5_chans = int16(L5_chans);

for j = 1 : size(ephys_data.conditions(end).unit_spikes,1) % for all units
    units(j).unit_depth = ephys_data.unit_depths(j); 
    units(j).L4_depth = ephys_data.LFP_min_chan*25;

    for k = 1 : numel(ephys_data.conditions) % for each condition
        whisk_onset = ephys_data.conditions(k).whisk_onset; % gets whisk onset for condition
        opto_onset = ephys_data.conditions(k).LED_onset; % gets LED onset for condition
        
        spikes = squeeze(ephys_data.conditions(k).unit_spikes(j,:,:)); % gets spikes for this unit for this condition
        opto_spikes = spikes - opto_onset; % gets opto spikes for this condition, this unit
        whisk_spikes = spikes - whisk_onset; % gets whisk spikes, for this coniditon this unit. 
        
        multi_spikes                     = ephys_data.conditions(k).multiunit_spikes;
        lc                               =  int16(units(j).unit_depth/25);
        clear local_chan
        local_chan                       =  lc-2 : lc+2 ; % 125 um range;
        local_chan(local_chan >32)       = [];
        local_chan(local_chan < 1)       = [];
        
        multispikes_local                = multi_spikes(local_chan,:,:);
        multispikes_l5                   = multi_spikes(L5_chans,:,:);
        
        opto_multispikes_local           = multispikes_local - opto_onset;
        opto_multispikes_l5              = multispikes_l5    - opto_onset;
        
        whisk_multispikes_local          = multispikes_local - whisk_onset;
        whisk_multispikes_l5             = multispikes_l5 - whisk_onset;
        
    
        [units(j).conditions(k).spont_behaviour] = Unit_behaviour_analysis(spikes,resp_win,multispikes_local,multispikes_l5); % gets unit behaviour for this condition in spont window
        [units(j).conditions(k).opto_behaviour] = Unit_behaviour_analysis(opto_spikes,resp_win,opto_multispikes_local,opto_multispikes_l5); % LED window
        [units(j).conditions(k).whisk_behaviour] = Unit_behaviour_analysis(whisk_spikes,resp_win,whisk_multispikes_local,whisk_multispikes_l5); %whisker window
       
        units(j).conditions(k).whisk_onset = whisk_onset;
        units(j).conditions(k).opto_onset = opto_onset;
        units(j).conditions(k).delta_t = opto_onset-whisk_onset;
        units(j).conditions(k).norm_trial_whisk_rate = units(j).conditions(k).whisk_behaviour.trial_spike_rate - units(j).conditions(k).spont_behaviour.trial_spike_rate; 
        units(j).conditions(k).norm_whisk_rate = units(j).conditions(k).whisk_behaviour.Mean_spike_rate;%-units(j).conditions(k).spont_behaviour.Mean_spike_rate;%
        units(j).conditions(k).norm_whisk_prob = units(j).conditions(k).whisk_behaviour.spike_prob;
        units(j).conditions(k).norm_whisk_1st_spike = nanmean(units(j).conditions(k).whisk_behaviour.first_spikes_by_trial);
        units(j).conditions(k).whisk_1st_spike_jitter = (nanstd(units(j).conditions(k).whisk_behaviour.first_spikes_by_trial))^2;
    end;
    for k = 1 : numel(ephys_data.conditions) % for each condition
    units(j).conditions(k).delta_whisk_rate = units(j).conditions(k).norm_whisk_rate/units(j).conditions(end).norm_whisk_rate;
   % units(j).conditions(k).trial_delta_whisk_rate = units(j).conditions(k).whisk_behaviour.trial_spike_rate-units(j).conditions(end).whisk_behaviour.trial_spike_rate;
    %[units(j).conditions(k).trial_rate_p,units(j).conditions(k).trial_rate_h] = ranksum(units(j).conditions(k).whisk_behaviour.trial_spike_rate,units(j).conditions(end).whisk_behaviour.trial_spike_rate);
    
    units(j).conditions(k).prob_delta = units(j).conditions(k).norm_whisk_prob - units(j).conditions(end).norm_whisk_prob;
   % test for change in probability of spiking vs control
    n_trials = size(units(j).conditions(k).whisk_behaviour.raster,1);
    trial_prob = int16(units(j).conditions(k).whisk_behaviour.spike_prob);
    control_prob = int16(units(j).conditions(end).whisk_behaviour.spike_prob);
    x = table([trial_prob*n_trials;n_trials-trial_prob*n_trials],[control_prob*n_trials;n_trials-control_prob*n_trials]);
    units(j).conditions(k).prob_p = fishertest(x);
    
    units(j).conditions(k).Spike_1_delta = units(j).conditions(k).norm_whisk_1st_spike - units(j).conditions(end).norm_whisk_1st_spike;
    end;
    
    %whisker responsiveness test for increase in probability of spiking
    n_trials = size(units(j).conditions(end).whisk_behaviour.raster,1);
    whisk_prob = int16(units(j).conditions(end).whisk_behaviour.spike_prob*n_trials);
    spont_prob = int16(units(j).conditions(end).spont_behaviour.spike_prob*n_trials);
    clear x;
    x = table([whisk_prob;n_trials-whisk_prob],[spont_prob;n_trials-spont_prob]);
    [whisk_resp_h,whisk_resp_p] = fishertest(x,'Tail','right');
    
    %Opto_responsiveness test for increase in probability of spiking
    n_trials = size(units(j).conditions(end).opto_behaviour.raster,1);
    opto_prob = int16(units(j).conditions(end).opto_behaviour.spike_prob*n_trials);
    clear x;
    x = table([opto_prob;n_trials-opto_prob],[spont_prob;n_trials-spont_prob]);
    [opto_resp_h,opto_resp_p] = fishertest(x,'Tail','right');
    
    
    whisk_resp_rate = units(j).conditions(end).whisk_behaviour.Mean_spike_rate;
    opto_resp_rate = units(j).conditions(end).opto_behaviour.Mean_spike_rate;
    
    
    units(j).whisk_responsive = (whisk_resp_p <= whisk_resp_p_thresh) & (whisk_resp_rate >= whisk_resp_rate_thresh);
    units(j).Opto_responsive = (opto_resp_p < opto_resp_p_thresh) & (opto_resp_rate >= opto_resp_rate_thresh);
    units(j).whisk_resp_p = whisk_resp_p;
    units(j).whisk_resp_rate = whisk_resp_rate;
    units(j).combined_response = [units(j).conditions(2).norm_trial_whisk_rate units(j).conditions(3).norm_trial_whisk_rate units(j).conditions(4).norm_trial_whisk_rate units(j).conditions(5).norm_trial_whisk_rate] ;
    [units(j).combined_rate_p,units(j).combined_rate_h] = ranksum(units(j).combined_response,units(j).conditions(end).whisk_behaviour.trial_spike_rate);
end;
end
