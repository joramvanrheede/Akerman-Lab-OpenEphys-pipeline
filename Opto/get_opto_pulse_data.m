function [pulse_data] = get_opto_pulse_data(data_folder, resp_win, psth_bins, artifact_win)
% [pulse_data] = get_opto_pulse_data(data_folder, channels, resp_win, psth_bins, artifact_win)
% retrieve optogenetic pulse spiking data from target data folder
% 
% 

expt_dirs                   = dir(data_folder);
expt_folders                = {expt_dirs.name};
qremove                     = ismember(expt_folders,{'.','..','.DS_Store'});
expt_folders(qremove)       = [];

for a = 1:length(expt_folders)

    pulse_data(a)           = opto_pulse_function(fullfile(data_folder,expt_folders{a}), resp_win, psth_bins, artifact_win);

end
