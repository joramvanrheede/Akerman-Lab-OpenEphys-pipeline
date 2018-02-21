% Process_channels_generic
% 
% Creates a new data structure, organised by date by experiment, that
% contains PSTH data for experiments
%
% It also contains all spike data in a 4d matrix (.all_spikes), organised as 
% condition * channel * episode * spike events
% 
% Intended to only process data for one experiment type at a time
% Loops through a target folder for a particular experiment type (e.g. 'Frequency')
% 

% location of preprocessed files
sdata_folder    = '/Users/Joram/Dropbox/Akerman Postdoc/Data/Extracted data 2018_02_07/16_02_Frequency';

%% Script parameters

% responsiveness criteria
whisk_resp_threshold          	= 1.5; % response threshold in x increase from baseline/spontaneous (i.e. 2 = 2x spontaneous rate)
LED_resp_threshold              = 1.5; % response threshold in x increase from baseline/spontaneous (i.e. 2 = 2x spontaneous rate)

% Parameters for analyse_sdata_function
analysisparams.whiskwinedges    = [-0.5:0.001:12]; % in seconds
analysisparams.LEDwin           = [0.005 0.030]; % in seconds
analysisparams.whiskwin         = [0 0.040]; % in seconds
analysisparams.LED_sust_win     = [0.100 0.400]; % in seconds
analysisparams.LEDwinedges      = [-0.250:0.005:0.750]; % in seconds
analysisparams.samplerate       = 30000; % in Hz
analysisparams.profile_smoothing = [0.03]; % size of gaussian window (in seconds) (span covers -3xSD to + 3xSD)

%% Running code starts here

date_folders = dir(sdata_folder);   % Get list of folders corresponding to the dates of experiments
date_folders = date_folders(3:end); % Get rid of '.' and '..' folders

% initialise 'sdata' struct
sdata	= [];
for a = 1:length(date_folders)
    
    % go to date folder and get file list
    this_date_folder    = date_folders(a).name;
    sdata_files         = dir([sdata_folder filesep this_date_folder]);
    
    % get rid of '.' and '..' directories
    sdata_files         = sdata_files(3:end);
    
    % get file names
    file_names          = {sdata_files.name};
    
    % Get rid of annoying '.DS_store' files that pop up everywhere
    if any(strcmp(file_names,'.DS_Store'))
        sdata_files     = sdata_files(~strcmp(file_names,'.DS_Store'));
    end
    
    % loop over data files
    for b = 1:length(sdata_files)
        this_sdata_file = sdata_files(b).name;
        
        % loads 'channels' and 'parameters' from saved data
        load([sdata_folder filesep this_date_folder filesep this_sdata_file])        
        disp([sdata_folder filesep this_date_folder filesep this_sdata_file])
        
        % Use analyse_channels_function to add PSTH data to 'channels'
        % struct
        channels            = analyse_channels_function(channels,analysisparams);
        
        % get conditions 
        condition_mat    	= cell2mat({channels(1).conditions.timings}');
        
        % store filename and the matrix of conditions in sdata.expt struct
        sdata(a).expt(b).filename                 	= this_sdata_file;
        sdata(a).expt(b).condition_mat             	= condition_mat;
        
        %% Start analysis by channel
        
        for c = 1:length(select_channels)
            this_channel                            = select_channels(c);
            
            % store current channel in temporary variable
            temp_channel                            = channels(this_channel).conditions;
            
            % store spontaneous rate of this channel in sdata.expt struct
            sdata(a).expt(b).spont_rate(c)         	= channels(this_channel).spontspikerate;
            
            % transfer all PSTH type data to the sdata.expt struct
            sdata(a).expt(b).whiskwinedges          = analysisparams.whiskwinedges;
            sdata(a).expt(b).LEDwinedges            = analysisparams.LEDwinedges;
            
            sdata(a).expt(b).whisk_counts(:,c)    	= [temp_channel.whisk_spike_count]';
            sdata(a).expt(b).whisk_rate(:,c)     	= [temp_channel.whisk_spike_rate]';
            sdata(a).expt(b).whisk_rel(:,c)       	= [temp_channel.whisk_spike_rel]';
            
            sdata(a).expt(b).whisk_peak_rate(:,c)  	= [temp_channel.whisk_peak_rate]';
            sdata(a).expt(b).whisk_peak_time(:,c)  	= [temp_channel.whisk_peak_time]';
            sdata(a).expt(b).whisk_profile(:,c,:) 	= [temp_channel.whisk_profile]';
            
            sdata(a).expt(b).whisk_win_counts(:,c,:)= [temp_channel.whisk_win_counts]';
            sdata(a).expt(b).whisk_win_rates(:,c,:) = [temp_channel.whisk_win_rates]';
            sdata(a).expt(b).whisk_win_rel(:,c,:)   = [temp_channel.whisk_rel_rates]';
            
            sdata(a).expt(b).LED_counts(:,c)      	= [temp_channel.LED_spike_count]';
            sdata(a).expt(b).LED_rate(:,c)        	= [temp_channel.LED_spike_rate]';
            sdata(a).expt(b).LED_rel(:,c)         	= [temp_channel.LED_spike_rel]';
            
            sdata(a).expt(b).LED_sust_counts(:,c) 	= [temp_channel.LED_sust_spike_count]';
            sdata(a).expt(b).LED_sust_rates(:,c)  	= [temp_channel.LED_sust_spike_rate]';
            sdata(a).expt(b).LED_sust_rel(:,c)    	= [temp_channel.LED_sust_spike_rel]';
            
            sdata(a).expt(b).LED_OFF_counts(:,c) 	= [temp_channel.LED_OFF_spike_count]';
            sdata(a).expt(b).LED_OFF_rates(:,c)  	= [temp_channel.LED_OFF_spike_rate]';
            sdata(a).expt(b).LED_OFF_rel(:,c)    	= [temp_channel.LED_OFF_spike_rel]';
            
            sdata(a).expt(b).LED_win_counts(:,c,:) 	= [temp_channel.LED_win_counts]';
            sdata(a).expt(b).LED_win_rates(:,c,:) 	= [temp_channel.LED_win_rates]';
            sdata(a).expt(b).LED_win_rel(:,c,:)   	= [temp_channel.LED_rel_rates]';
            
            % classify as 'responsive' if showing certain level of above
            % baseline activity
            % should this be done using standard deviations?
            sdata(a).expt(b).whisk_resp(:,c)        = [temp_channel.whisk_spike_rel]' > whisk_resp_threshold;
            sdata(a).expt(b).LED_resp(:,c)          = [temp_channel.LED_spike_rel]' > LED_resp_threshold;
            
            % make big matrix of spike times to facilitate further analysis
            % down the line (rasterplots etc):
            
            % loop over conditions
            for d = 1:length(temp_channel)
                temp_cond       = temp_channel(d);
                
                % loop over episodes
                for e = 1:length(temp_cond.episodes)
                    
                    % put all spike data for this experiment in one big
                    % matrix of spike times
                    sdata(a).expt(b).all_spikes(d,c,e,1:length(temp_cond.episodes(e).spikes)) = temp_cond.episodes(e).spikes;
                end
            end
            
            % the matrix will have lots of empty points, by default these
            % become zeros. Turn any value that is exactly zero into a NaN.
            sdata(a).expt(b).all_spikes(sdata(a).expt(b).all_spikes == 0)  = NaN;
            
        end
    end
end
