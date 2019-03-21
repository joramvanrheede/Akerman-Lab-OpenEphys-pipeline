function [cluster_spikes, cluster_ids, cluster_groups] = load_Kilosorted_spikes(spike_ID_npy_file, spike_time_npy_file, cluster_csv_file)
% function [cluster_spikes, cluster_ids, cluster_groups] = load_Kilosorted_spikes(spike_ID_npy_file, spike_time_npy_file, cluster_csv_file)
%
%
%
fid = fopen(cluster_csv_file);
headers         = textscan(fid,'%s%s',1);
cluster_data    = textscan(fid,'%d%s'); % read in cluster as numeric and grouping as string
fclose(fid);

cluster_ids     = cluster_data{1};
cluster_groups  = cluster_data{2};

spike_IDs       = readNPY(spike_ID_npy_file);
spike_times     = readNPY(spike_time_npy_file);


