% get_opto_resp

data_folder         = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data/2019_04_23/';
expt_folder         = 'RBSN_2019-04-23_14-07-22_8';


events_chans        = [1 1 3 2]; % order of events channels: [trial, whisk, opto, whisk_stim_nr]

spike_thresh        = 4;

%% 
channel_map         = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];


%% 
full_expt_folder    = fullfile(data_folder, expt_folder);


ephys_data          = get_ephys_data(full_expt_folder,events_chans,channel_map, spike_thresh);

