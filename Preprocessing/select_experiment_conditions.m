function select_experiment_conditions(experiment_filepath, save_dir, target_whisk_nr, target_opto_power)
% function select_experiment_conditions(experiment_filepath, save_dir, target_whisk_nr, target_opto_power)
% 
% Select only conditions with a TARGET_WHISK_NR and a TARGET_OPTO_POWER. 
% New experiment file with only these conditions will be saved in target SAVE_DIR.
% 

% Load experiment
load(experiment_filepath)

% Get parameters for the different conditions
cond_vals                   = ephys_data.condition_values;

% Column 6 contains whisker stim nr
whisk_stim_nrs              = cond_vals(:,6);
q_whisk_stim                = whisk_stim_nrs == target_whisk_nr;

% Column 8 contains opto power
opto_powers                 = cond_vals(:,8);
q_opto_power                = opto_powers == target_opto_power;

% Select only relevant conditions
ephys_data.conditions       = ephys_data.conditions(q_whisk_stim & q_opto_power);
ephys_data.condition_values = ephys_data.condition_values(q_whisk_stim & q_opto_power,:);

% Deconstruct full file name
[dirnm, filenm, extn]     	= fileparts(experiment_file);

% Save in target location
save([save_dir filesep filenm extn 'w' num2str(target_whisk_nr) '_pw' num2str(target_opto_power)], 'ephys_data')

