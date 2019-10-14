function sync_curated_data(Kilosort_dir, data_output_dir)
% function sync_curated_data(KILOSORT_DIR, DATA_OUTPUT_DIR)
% 
% Matches up curated sorted spikes in KILOSORT_DIR with sync_data file to 
% generate a data structure organised by experiment, by experimental condition.
% 
% Saves the data structure by experiment in a hierarchically organised folder
% structure (DATA_OUTPUT_DIR/'experiment name'/'date_implantation_nr').
% 

% loads data structure sync_data.mat and matches up with Kilosort spike times
sorted_ephys_data   = sync_Kilosort_spikes(Kilosort_dir);

% distributes the data by experiment type and by date in DATA_OUTPUT_DIR
distribute_sorted_data(sorted_ephys_data, data_output_dir)