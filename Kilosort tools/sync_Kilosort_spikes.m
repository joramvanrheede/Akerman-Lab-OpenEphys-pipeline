function [sorted_ephys_data] = sync_Kilosort_spikes(Kilosort_dir)
% function [sorted_ephys_data] = sync_Kilosort_spikes(Kilosort_dir)
% Takes output from Kilosort and pre-made sync_data structure and puts them together
% sorted_ephys_data should be equivalent to the 'ephys_data' structure for multichannel
% multiunit data, only in 'spikes', the first dimension is sorted 'good' units rather
% than channels.
% 

%% Constant - unlikely to change unless we change acquisition system
samplerate      = 30000;

%% Load Kilosorted data using CortexLab/spikes functions:

sp              = loadKSdir(Kilosort_dir); % CortexLab/spikes function to read in Kilosorted, phy-curated data quickly

% Unpack some of the loadKSdir data
spike_times     = sp.st;    % Spike times for all spikes, N spikes * 1 vector
cluster_ids     = sp.clu;   % cluster ID (number) for each spike time, N spikes * 1 vector
cluster_numbers = sp.cids;  % IDs of all clusters (1 * N clusters vector of ascending integers, starting at 0, with potential gaps for merged clusters)
cluster_groups  = sp.cgs;   % Group corresponding to each cluster_number, 1 * N clusters vector --> 1 = MUA, 2 = 'Good' unit

% CortexLab/spikes function to get more info about the clusters using the loadKSdir data
% We need spikeDepths, a N spikes * 1 vector, to determine the depth of the clusters
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

%% Select 

% Get depth of cluster by identifying the unique cluster / depth pairs
cluster_depths = unique([cluster_ids, spikeDepths],'rows');

% boolean for selecting units
is_unit         = cluster_groups == 2;

% get cluster numbers and depths for confirmed 'Good' units only
unit_clusters   = cluster_numbers(is_unit);
unit_depths     = cluster_depths(is_unit,2);

% Sort by depth, from superficial to deep (this will be the order in which units are added to the sorted_ephys_data struct)
[sort_depths, depth_order]  = sort(unit_depths);

% Make boolean to select only spikes that came from confirmed 'Good' units
is_unit_spike   = ismember(cluster_ids,unit_clusters);

% Select only spikes and corresponding cluster IDs from good units
spike_times     = spike_times(is_unit_spike);
cluster_ids     = cluster_ids(is_unit_spike);

% See which clusters are left
uniq_clusters   = unique(cluster_ids);
n_clusters      = length(uniq_clusters);

%% Load files for stimulus sync
load(fullfile(Kilosort_dir,'sync_data.mat')); % loads 'sync_data' variable
load(fullfile(Kilosort_dir,'concatenated_data_medianTrace.mat')) % loads medianTrace variable

% Copy to new struct that will be populated with spike data as the output of this function
sorted_ephys_data   = sync_data;

% time offset for each struct
rec_durations       = [sorted_ephys_data.rec_duration];


%% Hack - use the median trace to identify gaps between recordings to correct for accumulating spike time delay
% This requires common average reference to have been run.
% This will be problematic if common average reference is not carried out -
% Are there any situations in which this would be the case?
% 
% Essentially, the problem is this: When concatenating recordings, the end 
% of each recording will be padded with (zero) values to fill up a batch size.
% This leads to an increasing temporal offset between sync_data and the kilosorted
% spike times with an increasing spike time 'delay' in later recordings.
%
% Current hack: identify instances where the common average reference trace
% is completely flat for a long period of time. Determine the length of this
% gap and subtract it from the spike times for each subsequent recording.
%
% Not a robust solution - This is bound to occasionally break down when the
% recording ends (nearly) exactly on the end of a batch.

median_no_diff    	= diff(medianTrace) == 0; % true for where there are identical values in a row

% find instances where the gaps start and end
diff_one            = diff(median_no_diff == 1);
start_gap           = find(diff_one == 1);
end_gap             = find(diff_one == -1);

% Each start has to have a subsequent end
if start_gap(end) > end_gap(end)
    start_gap(end) = [];
end

% Build a cumulative duration variable tha
rec_end_times       = NaN(size(rec_durations));
cum_duration        = 0;
for a = 1:length(rec_durations)-1

    this_duration               = rec_durations(a);
    
    total_duration_temp         = this_duration + cum_duration; % Add to cumulative duration
    
    total_duration_temp_ind     = total_duration_temp * samplerate; % find index by multiplying by sample rate
    min_ind                     = total_duration_temp_ind - samplerate; % -1sec
    max_ind                     = total_duration_temp_ind + samplerate; % +1sec
    
    these_start_gaps            = start_gap(start_gap > min_ind & start_gap < max_ind);
    these_end_gaps              = end_gap(end_gap > min_ind & start_gap < max_ind);
    
    if these_start_gaps(end) > these_end_gaps(end)
        these_start_gaps(end)   = [];
    end
    
    if these_end_gaps(1) < these_start_gaps(1)
        these_end_gaps(1)       = [];
    end
    
    flat_line_lengths           = these_end_gaps - these_start_gaps; % how long is CAR trace flatlining?
    [max_gap, max_ind]          = max(flat_line_lengths); % Find max gap
    rec_end_ind                 = these_end_gaps(max_ind); % 
    rec_end_times(a)            = rec_end_ind / samplerate;
    
    cum_duration                = rec_end_times(a); 
end

rec_end_times(end)  = numel(medianTrace) / samplerate;
rec_start_times     = [0 rec_end_times(1:end-1)];

%% Assign Kilosorted spike times to trials and conditions in sorted_ephys_data

for a = 1:length(sorted_ephys_data)
    disp(['Loading recording nr ' num2str(a) ' of ' num2str(length(sorted_ephys_data))])
    
    Kilosort_rec_start_time     = rec_start_times(a);
    original_rec_start_time     = sorted_ephys_data(a).rec_start_time;
    for b = 1:length(sorted_ephys_data(a).conditions)
        for c = 1:length(sorted_ephys_data(a).conditions(b).trial_starts)
            for d = 1:n_clusters
                this_cluster    = uniq_clusters(depth_order(d));
                q_cluster       = cluster_ids == this_cluster;
                
                % Make relative to Kilosort concatenated data time, not openephys recording time stamp
                trial_start     = sorted_ephys_data(a).conditions(b).trial_starts(c) - original_rec_start_time + Kilosort_rec_start_time;
                trial_end       = sorted_ephys_data(a).conditions(b).trial_ends(c) - original_rec_start_time + Kilosort_rec_start_time;
                
                q_trial         = spike_times >= trial_start & spike_times < trial_end;
                
                % get relevant spike times and make them relative to trial onset
                unit_trial_spikes   = spike_times(q_cluster & q_trial) - trial_start; 
                
                % Assign spikes for cluster to this condition
                sorted_ephys_data(a).conditions(b).spikes(d,c,1:length(unit_trial_spikes))  = unit_trial_spikes;
                
            end
        end
        
        % set empty values to NaN instead of 0
        sorted_ephys_data(a).conditions(b).spikes(sorted_ephys_data(a).conditions(b).spikes == 0) = NaN;
    end
    
    % Add unit depth information
    sorted_ephys_data(a).unit_depths = sort_depths;
end


