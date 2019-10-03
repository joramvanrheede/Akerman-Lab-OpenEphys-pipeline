function [timing_data] = get_opto_timing_data(data_folder, resp_win, psth_bins)
% function [timing_data] = get_opto_timing_data(data_folder, resp_win, psth_bins)

expt_dirs                   = dir(data_folder);
expt_folders                = {expt_dirs.name};
qremove                     = ismember(expt_folders,{'.','..','.DS_Store'});
expt_folders(qremove)       = [];

timing_data                 = struct;
for a = 1:length(expt_folders)
    this_folder     = expt_folders{a};
    
    fullfolder      = fullfile(data_folder, this_folder);
    this_exp_files  = dir([fullfolder '/*.mat']);
    
    for b = 1:length(this_exp_files)
        
        this_expt_name                  = this_exp_files(b).name;
        disp(['Loading ' this_expt_name '...'])
        load(fullfile(fullfolder, this_expt_name));
        
        % 
        timing_data(a).experiment(b)    = timing_function(ephys_data, resp_win, psth_bins);

    end
end