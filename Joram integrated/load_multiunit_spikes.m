function groups = load_multiunit_spikes(group_folders)

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
            disp(['Loading ' expt_file_path])
            
            groups(a).prep(b).expt_data(c)    = load(expt_file_path);
            for d = 1:length(groups(a).prep(b).expt_data(c).ephys_data.conditions)
                
                if d == length(groups(a).prep(b).expt_data(c).ephys_data.conditions)
                    cond_data = groups(a).prep(b).expt_data(c).ephys_data.conditions(d);
                    
                    [CSD,w,sinks,sources,analysis_win] = Current_Source_Density(cond_data.LFP_trace(:,:,900:1200),[1 200]);
                    
                    LFP_means       = squeeze(mean(cond_data.LFP_trace,2));
                    LFP_min_profile = min(LFP_means(:,900:1200),[],2);
                    LFP_min_profile = smooth(LFP_min_profile,5);
                    LFP_min_chan  	= find(LFP_min_profile == min(LFP_min_profile));
                    
                    sink_profile    = smooth(sinks,5);
                    max_sink_chan  	= find(sink_profile == min(sink_profile))+1;
                    
                    groups(a).prep(b).expt_data(c).ephys_data.sink_profile      = sink_profile;
                    groups(a).prep(b).expt_data(c).ephys_data.max_sink_chan     = max_sink_chan;
                    
                    groups(a).prep(b).expt_data(c).ephys_data.LFP_min_profile 	= LFP_min_profile;
                    groups(a).prep(b).expt_data(c).ephys_data.LFP_min_chan      = LFP_min_chan;
                end
                
                rmfield(groups(a).prep(b).expt_data(c).ephys_data.conditions, 'LFP_trace');
                
            end
            toc
        end
    end
end
