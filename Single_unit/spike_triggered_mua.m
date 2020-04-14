function [triggered_psths, norm_psths, n_valid_spikes] = spike_triggered_mua(ref_spikes, multi_unit_spikes, resp_win, psth_bins)

% function [TRIGGERED_PSTHS, NORM_PSTHS, N_VALID_SPIKES] = SPIKE_TRIGGERED_MUA(REF_SPIKES, MULTI_UNIT_SPIKES, RESP_WIN, PSTH_BINS)
% 
% Makes a spike-triggered average post-stimulus time histogram for each unit
% or channel in REF_SPIKES.
% For each trial in REF_SPIKES, the FIRST spike found within RESP_WIN is used
% to make the TRIGGERED_PSTHS.
% 
% INPUTS:
% 
% REQUIRED:
% 
% REF_SPIKES: An U x M x N matrix of spike times where U = n_units / channels;
% M = n_trials and N = n_spikes (padded with NaNs for empty positions). 
% 
% MULTI_UNIT_SPIKES: An U x M x N matrix of spike times where U = n_units / 
% channels; M = n_trials and N = n_spikes (padded with NaNs for empty positions).
% Multiunit activity is taken from ALL channels / units U. The number of trials
% should correspond to the number of trials for REF_SPIKES.
% 
% OPTIONAL:
% 
% RESP_WIN: A 2-element vector indicating the start and end time of the response
% window. Defaults to [0.006 0.030].
% 
% PSTH_BINS: Bin edges for the PSTH. Defaults to [-0.1:0.001:0.1] for 1ms bins
% ranging from -100ms to +100ms
% 
% OUTPUTS:
% 
% TRIGGERED_PSTHS: an M x N matrix of post stimulus time histogram bins where
% M = n_units / n_channels as determined by the number of units / channels
% in REF_SPIKES, and N = the number of psth bins. There will be NaNs if no 
% REF_SPIKES were found within RESP_WIN for any trials such that a spike 
% triggered average of multiunit activity could not be generated.
% 
% NORM_PSTHS: TRIGGERED_PSTHS normalised such that the total of each row adds 
% up to 100 (or NaN if no REF_SPIKES were found within RESP_WIN for any trials 
% and a spike triggered average of multiunit activity could not be generated).
% 
% N_VALID_SPIKES: the number of valid spikes found for each unit, ranging 
% from 0 to n_trials, the number of trials in REF_SPIKES & MULTIUNIT_SPIKES. 
% 
% 

% Set default psth_bins from -100ms to +100ms in 1ms increments
if nargin < 4
    psth_bins   = [-0.1:0.001:0.1];
end

% Set default response window
if nargin < 3 || isempty(resp_win)
    resp_win    = [0.006 0.030];
end

% Pre-allocate psth counts (NaNs will remain in trials where there was no reference spike)
psth_counts         = NaN(size(ref_spikes,1),size(ref_spikes,2),length(psth_bins)-1);

% Loop over all trials
for a = 1:size(ref_spikes,2)
    
    % Get the multiunit spikes for this trial, across all channels
    trial_multi_spikes      = multi_unit_spikes(:,a,:);
    trial_multi_spikes      = trial_multi_spikes(:);
    
    % Loop over all units / channels in ref_spikes
    for b = 1:size(ref_spikes,1)
        % Get ref_spikes for this trial and this unit
        unit_trial_spikes       = squeeze(ref_spikes(b,a,:));
        
        % Select spikes in response window
        q_resp_win              = unit_trial_spikes > resp_win(1) & unit_trial_spikes < resp_win(2);
        resp_spikes             = unit_trial_spikes(q_resp_win);
        
        % Select first spike time in response window
        first_spike             = min(resp_spikes);
        
        % If there is indeed a (first) spike...
        if ~isempty(first_spike)
            % Subtract the first spike time from the multiunit spikes to get
            % spike times relative to this reference spike
            spike_triggered_spikes      = trial_multi_spikes - first_spike;
            
            % Make a spike-triggered PSTH for this trial
            [~, psth_counts(b,a,:)]     = psth(spike_triggered_spikes,psth_bins,0);
        end
    end
end

% Sum over trials to generate the n_units x n_psth_bins output psths
triggered_psths = squeeze(nansum(psth_counts,2));

% Normalise psths so that they sum to 100 (for easy cross-comparison and e.g. 
% determining percentile location of the reference spikes)
norm_psths      = squeeze(nanmean(psth_counts,2));
norm_psths      = norm_psths ./ sum(norm_psths,2) * 100;

% Use the absence of NaNs to establish where there were and were not any valid ref spikes
n_valid_spikes = sum(~isnan(psth_counts(:,:,1)),2);

