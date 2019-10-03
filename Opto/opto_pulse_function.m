function pulse_data = opto_pulse_function(data_dir, resp_win, psth_bins, artifact_win)
% pulse_data = opto_pulse_function(data_dir, channels, resp_win, psth_bins, artifact_win)
% Laser_Pulse analysis


if exist('artifact_win','var')
    kill_artifact   = 1;
else
    kill_artifact   = 0;
end

if exist('resp_thresh','var')
    do_thresh       = true;
end

%%

folder_files    = dir([data_dir filesep '*.mat']); % read .mat files from folder

all_count_data  = [];
for a = 1:length(folder_files)
    
    this_file       = folder_files(a).name;
    full_filenm     = fullfile(data_dir, this_file);
    
    disp(['Loading ' this_file]);
    load(full_filenm);  % this should load variable 'ephys_data' and'
    
    laser_power(a)  = ephys_data.parameters.LED_power;  % What laser/LED power was used?
    spikes          = ephys_data.conditions(1).spikes;  % get spikes
    opto_onset      = ephys_data.conditions(1).LED_onset; % when did opto stimulus come on?
    
    spont_spikes    = spikes;
    spikes          = spikes - opto_onset; % make spike times relative to opto_onset
    
    spont_win       = [0 (opto_onset - 0.005)]; % get spontaneous spiking from entire time win before stimulus onset
    
    if kill_artifact
        q_artifact      = spikes >= artifact_win(1) & spikes <= artifact_win(2); %
        spikes(q_artifact) = NaN;
    end
    
    [im_h, laser_psth(a,:)]             = psth(spikes,psth_bins,false); % get psth in chose time window
    
    spike_resp_spont(:,a)            	= spike_rate_by_channel(spont_spikes,spont_win);
    peak_ROF_spont(:,a)                 = peak_ROF_by_channel(spont_spikes,spont_win);
    for b = 1:size(resp_win,2)
        [raw_spike_resp(:,a,b), spike_std(:,a,b)]  = spike_rate_by_channel(spikes,resp_win(b,:));
        corr_spike_resp(:,a,b)                      = raw_spike_resp(:,a,b) - spike_resp_spont(:,a);
        
        [raw_peak_ROF(:,a,b), peak_time(:,a,b)]     = peak_ROF_by_channel(spikes,resp_win(b,:));
        corr_peak_ROF(:,a,b)                        = raw_peak_ROF(:,a,b) - spike_resp_spont(:,a); % correct peak by baseline spike rate (not by baseline peak)
    end
    
    % Get spike counts over time by channel (don't show spike density plot)
    [imh, all_count_data(:,:,a)]     	= spike_density_plot(spikes,1,psth_bins, false);
end

uniq_powers     = unique(laser_power);
n_powers        = length(uniq_powers);

%% Loop to make PSTH and get response for each laser power
for b = 1:n_powers
    this_power                          = uniq_powers(b);
    q_power                             = laser_power == this_power;
    
    corr_spike_resp_by_power(:,b,:)     = nanmean(corr_spike_resp(:,q_power,:),2);
    raw_spike_resp_by_power(:,b,:)      = nanmean(raw_spike_resp(:,q_power,:),2);
    spont_spike_resp_by_power(:,b)      = nanmean(spike_resp_spont(:,q_power),2);
    
    corr_peak_ROF_by_power(:,b,:)     	= nanmean(corr_peak_ROF(:,q_power,:),2);
    raw_peak_ROF_by_power(:,b,:)     	= nanmean(raw_peak_ROF(:,q_power,:),2);
    peak_time_by_power(:,b,:)           = nanmean(peak_time(:,q_power,:),2);
    spont_peak_ROF_by_power(:,b)        = nanmean(peak_ROF_spont(:,q_power),2);
    
    
    spike_std_by_power(:,b,:)           = nanmean(spike_std(:,q_power,:),2);
    
    these_psths                         = laser_psth(q_power,:);     % Select all PSTHs for this laser power
    mean_psth_by_power(b,:)             = mean(these_psths,1);    % Make mean PSTH trace
    
    these_counts                        = all_count_data(:,:,q_power);
    mean_counts_by_power(:,:,b)         = nanmean(these_counts,3);
end

pulse_data.opto_power           = uniq_powers;
pulse_data.resp_win             = resp_win;
pulse_data.delta_spike_rate     = corr_spike_resp_by_power;
pulse_data.raw_spike_rate     	= raw_spike_resp_by_power;
pulse_data.spont_spike_rate     = spont_spike_resp_by_power;
pulse_data.spike_std            = spike_std_by_power;
pulse_data.peak_spike_rate      = raw_peak_ROF_by_power;
pulse_data.corr_peak_spike_rate = corr_peak_ROF_by_power;
pulse_data.peak_time            = peak_time_by_power;
pulse_data.psths                = mean_psth_by_power;
pulse_data.density_counts       = mean_counts_by_power;



