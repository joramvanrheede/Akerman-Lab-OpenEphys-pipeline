
% Load kilosorted units
% Load synced structs
% /Volumes/Akermanlab/Joram/Spike_sorting/2019_03_09-_1_3_4_5_6


Kilosort_dir    = ['/Volumes/Akermanlab/Joram/Spike_sorting/2019_09_16-_12'];
raw_data_dir    = []; % 

do_units_only   = true;

%% Constants
samplerate      = 30000;

%% Load Kilosorted data 
sp              = loadKSdir(Kilosort_dir);

spike_times     = sp.st;
cluster_ids     = sp.clu;
cluster_groups  = sp.cgs;
cluster_numbers = sp.cids;




[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

% Get depth of cluster
cluster_depths = unique([cluster_ids, spikeDepths],'rows');

is_unit         = cluster_groups == 2;
is_mua          = cluster_groups == 1;

unit_clusters   = cluster_numbers(is_unit);
unit_depths     = cluster_depths(is_unit,2);

[sort_depths, depth_order]  = sort(unit_depths);

is_unit_spike   = ismember(cluster_ids,unit_clusters);

if do_units_only
    spike_times     = spike_times(is_unit_spike);
    cluster_ids     = cluster_ids(is_unit_spike);
end

uniq_clusters   = unique(cluster_ids);
n_clusters      = length(uniq_clusters);

%% Load files for stimulus sync
load(fullfile(Kilosort_dir,'sync_data.mat')); % loads 'sync_data' variable
load(fullfile(Kilosort_dir,'concatenated_data_medianTrace.mat')) % loads medianTrace variable

% time offset for each struct
rec_durations   = [sync_data.rec_duration];



%% 
median_no_diff    	= diff(medianTrace) == 0; % true for where there are identical values in a row

diff_one            = diff(median_no_diff == 1);
start_gap           = find(diff_one == 1);
end_gap             = find(diff_one == -1);

if start_gap(end) > end_gap(end)
    start_gap(end) = [];
end

rec_end_times       = NaN(size(rec_durations));
cum_duration        = 0;
for a = 1:length(rec_durations)-1

    this_duration               = rec_durations(a);
    
    total_duration_temp         = this_duration + cum_duration;
    
    total_duration_temp_ind     = total_duration_temp * samplerate; % find index by multiplying by sample rate
    min_ind                     = total_duration_temp_ind - samplerate; % -10sec
    max_ind                     = total_duration_temp_ind + samplerate; % +10sec
    
    these_start_gaps            = start_gap(start_gap > min_ind & start_gap < max_ind);
    these_end_gaps              = end_gap(end_gap > min_ind & start_gap < max_ind);
    
    if these_start_gaps(end) > these_end_gaps(end)
        these_start_gaps(end)   = [];
    end
    
    if these_end_gaps(1) < these_start_gaps(1)
        these_end_gaps(1)       = [];
    end
    
    flat_line_lengths           = these_end_gaps - these_start_gaps;
    [max_gap, max_ind]          = max(flat_line_lengths);
    rec_end_ind                 = these_end_gaps(max_ind);
    rec_end_times(a)            = rec_end_ind / samplerate;
    
    cum_duration                = rec_end_times(a);
end

rec_end_times(end)  = numel(medianTrace) / samplerate;
rec_start_times     = [0 rec_end_times(1:end-1)];

%% 

for a = 1:length(sync_data)
    disp(['Loading recording nr ' num2str(a) ' of ' num2str(length(sync_data))])
    
    Kilosort_rec_start_time     = rec_start_times(a);
    original_rec_start_time     = sync_data(a).rec_start_time;
    for b = 1:length(sync_data(a).conditions)
        for c = 1:length(sync_data(a).conditions(b).trial_starts)
            for d = 1:n_clusters
                this_cluster    = uniq_clusters(depth_order(d));
                q_cluster       = cluster_ids == this_cluster;
                
                % Make relative to Kilosort concatenated data time, not openephys recording time stamp
                trial_start     = sync_data(a).conditions(b).trial_starts(c) - original_rec_start_time + Kilosort_rec_start_time;
                trial_end       = sync_data(a).conditions(b).trial_ends(c) - original_rec_start_time + Kilosort_rec_start_time;
                
                q_trial         = spike_times >= trial_start & spike_times < trial_end;
                
                % get relevant spike times and make them relative to trial onset
                unit_trial_spikes   = spike_times(q_cluster & q_trial) - trial_start; 
                
                
                sync_data(a).conditions(b).spikes(d,c,1:length(unit_trial_spikes))  = unit_trial_spikes;
                
            end
        end
        
        % set empty values to NaN instead of 0
        sync_data(a).conditions(b).spikes(sync_data(a).conditions(b).spikes == 0) = NaN;
    end
end

sync_data(a).unit_depths = sort_depths;
