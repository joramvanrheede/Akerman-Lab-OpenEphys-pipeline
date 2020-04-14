function [Mean_Light_Delta] = Layer_response(cond_data,channels);
    lowest_power = 1;
        opto_powers     = [cond_data(:).LED_power]';
    opto_powers(isnan(opto_powers)) = lowest_power;
    artifact_win    = [-0.002 0.002];
     
    % if any powers were NaN and replaced with a minimum nr, powers and conditions struct
    % will not be in ascending order; fix:
    [opto_powers]           = unique(opto_powers);
    cond_data_powers      	= [cond_data.LED_power];
    cond_data_powers(isnan(cond_data_powers)) = lowest_power;

    
    no_Powers = numel(opto_powers);
    opto_onset              = nanmedian([cond_data(:).LED_onset]);
    opto_duration           = nanmedian([cond_data(:).LED_duration]);         
    resp_win        = [artifact_win(2)+opto_duration artifact_win(2)+opto_duration+0.02];

    Light_evoked = NaN*ones(no_Powers,200);
    Spont = NaN*ones(no_Powers,200);
    Light_Delta =NaN*ones(no_Powers,200);

    for k = 1 : no_Powers
    this_power                  = opto_powers(k);
    Powers(k) = this_power;
    
    q_power                     = cond_data_powers == this_power;
    
    
    power_inds                  = find(q_power);
    spikes                      = [];
    for c = 1:length(power_inds)
        
        these_spikes                = cond_data(power_inds(c)).spikes(channels,:,:);

        spikes(1:size(these_spikes,1),size(spikes,2)+(1:size(these_spikes,2)),1:size(these_spikes,3)) = these_spikes;
        
    end
    spikes(spikes == 0)     = NaN; % Default var gets padded with zeros; change padding zeros to NaN
    
    spont_spikes            = spikes;
    spikes                  = spikes - opto_onset;
    
    n_trials                = size(spikes,2);
       
        %% removes artifact
        q_artifact            	= spikes > artifact_win(1) & spikes < artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > artifact_win(1) & spont_spikes < artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
        q_artifact            	= spikes > opto_duration + artifact_win(1) & spikes < opto_duration +artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > opto_duration + artifact_win(1) & spont_spikes < opto_duration +artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
        
            % Binned spike rate
            spike_rate(1:n_trials)                               = spike_rate_by_trial(spikes(:,:,:), resp_win);
    
            % Spontaneous spike rate
            spont_rate(1:n_trials)                               = spike_rate_by_trial(spont_spikes(:,:,:), resp_win);
     
      Light_evoked(k,1:numel(spike_rate)) = spike_rate;
       Spont(k,1:numel(spont_rate)) = spont_rate;
      Light_Delta(k,1:numel(spont_rate)) = spike_rate-spont_rate;
    
    Mean_Light_Delta(k) = nanmean(Light_Delta(k,:));
   
    

end;
end