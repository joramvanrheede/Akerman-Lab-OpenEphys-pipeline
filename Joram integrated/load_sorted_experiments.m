function groups = load_sorted_experiments(group_folders)

tic

% Loop over experiment groups (e.g. POM / M1 / S1)
for a = 1:length(group_folders)
    
    prep_folders    = dir(group_folders{a}); %
    prep_names      = {prep_folders([prep_folders.isdir]).name}; % Get directories only
    prep_names(ismember(prep_names, {'.' '..'})) = []; % remove . and ..
    
    % Loop over experimental preps (usually organised by date)
    for b = 1:length(prep_names)
        prep_dir    = fullfile(group_folders{a},prep_names{b});
        
        experiments = dir(fullfile(prep_dir,'*.mat'));
        
        % Loop over experiments within prep
        for c = 1:length(experiments)
            
            expt_file_path  = fullfile(prep_dir,experiments(c).name);
            disp(['Loading ' expt_file_path '...'])
            toc
            groups(a).prep(b).expt_data(c)    = load(expt_file_path);
            
            %
        end
    end
end
