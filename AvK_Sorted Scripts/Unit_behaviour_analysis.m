function Win_behaviour = Unit_behaviour_analysis(spikes,resp_win,multi_spikes_local,multi_spikes_l5);
% Win_behaviour  : returns Unit Behavioural measures from Spikes() for a
% given condition during the time window specified by resp_win
%% INPUTS
 % spikes : n(trials)x m (spiketimes) matrix with all spiketimes for a unit
 %          in a given condition.
 %resp_win : time window in seconds to assess; 
  
 %% OUTPUTS
 % Win_behaviour strucure: 
    %Mean_spike_rate         : the mean spike rate in Hz across all trials
    %Mean_spike_isi          : the mean spike isi in seconds across all trials
    %Mean_spike_probablity   : Probability of spiking in given window 
    %Mean_ 1st spike_latency : time to first spike in the time window. 
    %
    %Trial_spike_rate        : n (number of trials) x1 array with the spike rate
    %                           in Hz for all trials during the given response
    %                           window. 
    %% create PSTH and Raster for single unit spiking activity for 100ms from onset of window
    %artifact_win                        = 0.005;
   [Win_behaviour.spike_prob]           = spike_prob_by_channel(spikes, resp_win);
   spikes_100                           = spikes;
   spikes_in_win_100                    = spikes >= resp_win(1) & spikes <= resp_win(1)+0.1;
   spikes_100(~spikes_in_win_100)       =NaN;
   [Win_behaviour.raster_100]           = spikes_100; 
   [Win_behaviour.psth_100]             = trial_bins_singleu(Win_behaviour.raster_100,0.005,[resp_win(1) resp_win(1)+0.1]);
   Win_behaviour.cum_psth = cumsum(sum(Win_behaviour.psth_100,1)); %cumulative sum of spikes in 100ms
   
   %% create PSTH and Raster for both local multiunit and l5 activity 
   multi_spikes_local_100                             = multi_spikes_local;
   multi_spikes_in_win_100                            = multi_spikes_local >= resp_win(1) & multi_spikes_local <= resp_win(1)+0.1;
   multi_spikes_local_100(~multi_spikes_in_win_100)    = NaN;
   [Win_behaviour.multi_local_raster_100]             = multi_spikes_local_100; 
   [~, Win_behaviour.multi_local_psth_100, ~] = psth_by_trial(multi_spikes_local_100,0.005,[resp_win(1) resp_win(1)+0.1],false);
   clear multi_spikes_in_win_100;
   
   multi_spikes_l5_100                                = multi_spikes_l5;
   multi_spikes_in_win_100                            = multi_spikes_l5 >= resp_win(1) & multi_spikes_l5 <= resp_win(1)+0.1;
   multi_spikes_l5_100(~multi_spikes_in_win_100)      = NaN;
   [Win_behaviour.multi_l5_raster_100]                = multi_spikes_l5_100;
   [~, Win_behaviour.multi_l5_psth_100, ~]            = psth_by_trial(multi_spikes_l5_100,0.005,[resp_win(1) resp_win(1)+0.1],false);
   
   
   
   spikes_in_win                       = spikes >= resp_win(1) & spikes <= resp_win(2); % boolean of which spike times fall within window
   spike_times_in_win                  = spikes; % Copy original spike times
   spike_times_in_win(~spikes_in_win)  = NaN; % Set spike times out of window to 0

                                              % Find first spike in window by trial by channel
   Win_behaviour.raster = spike_times_in_win;
   [first_spikes_by_trial,i]    = nanmin(spike_times_in_win,[],2);
    spike_times_in_win(:,i) = NaN;
   [second_spikes_by_trial]    = nanmin(spike_times_in_win,[],2);
   
    
   Win_behaviour.first_spikes_by_trial = first_spikes_by_trial;
   Win_behaviour.second_spikes_by_trial = second_spikes_by_trial;
   Win_behaviour.first_isi = second_spikes_by_trial - first_spikes_by_trial;
   
   Win_behaviour.Mean_Spike_latency = nanmean(first_spikes_by_trial);
   
     for k = 1 : size(spikes,1)
            trial  = squeeze(spikes(k,:)); % gets all unit spi
            trial = trial(resp_win(1) < trial & resp_win(2) > trial);
            trial_isi(k) = nanmean(diff(trial));
            trial_spike_rate(k) = numel(trial)/(resp_win(2)-resp_win(1));
    end;
  Win_behaviour.Mean_isi = nanmean(trial_isi);
  Win_behaviour.Mean_spike_rate = nanmean(trial_spike_rate);
  Win_behaviour.trial_spike_rate = trial_spike_rate;
  Win_behaviour.rate_norm = kstest(trial_spike_rate);
end
