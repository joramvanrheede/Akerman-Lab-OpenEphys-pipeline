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

experiment_type     = 'Frequency';
experiment_folder   = '/Users/Joram/Dropbox/Akerman Postdoc/Data/Extracted data/Dual Stim';

% location of preprocessed files
sdata_folder        = [experiment_folder filesep experiment_type]; % '/Users/Joram/Dropbox/Akerman Postdoc/Data/Extracted data 2018_02_07/Dual stim test';

% responsiveness criteria
whisk_resp_threshold                = 1.5; % response threshold in x increase from baseline/spontaneous (i.e. 2 = 2x spontaneous rate)
LED_resp_threshold                  = 1.5; % response threshold in x increase from baseline/spontaneous (i.e. 2 = 2x spontaneous rate)

%% Script parameters

% set general parameters for use in analyse_channels_function
analysisparams.samplerate           = 30000;            % data sample rate in Hz
analysisparams.profile_smoothing    = [0.02];           % size of gaussian window (in seconds) (span covers -3xSD to + 3xSD)

analysisparams.whiskwin             = [0 0.040];        % window for assessing LED response
analysisparams.whiskwinedges        = [-1:0.001:3];     % histogram bin edges for whisker PSTH

analysisparams.LEDwin               = [0.002 0.040];    % in seconds
analysisparams.LED_sust_win         = [1 2]; % in seconds
analysisparams.LEDwinedges          = [-0.250:0.001:3]; % in seconds

% adjust individual parameters dependent on experiment_type
switch experiment_type
    case 'Frequency'
        analysisparams.LED_sust_win         = [1 4];            % LED is on for longer, can take larger sample
        analysisparams.whiskwinedges        = [-1:0.001:5];     % multiple whisks - needs PSTH with more time bins
        analysisparams.LEDwinedges          = [-0.5:0.001:6];   % Longer LED stimulation needs PSTH with more time bins
    case 'Drive'
        % no adjustments needed from general settings?
    case 'Timing'
        analysisparams.LEDwin               = [0.002 0.025];    % LED pulse only lasts 25ms
        analysisparams.LED_sust_win         = [0.002 0.025];    % no real sustained response here
        analysisparams.LEDwinedges          = [-0.5:0.001:2];   % shorten PSTH
    case 'Velocity'
        analysisparams.whiskwin             = [0 1];            % slow whisks need longer time bin to capture peak response
end

%% Running code starts here

date_folders            = dir(sdata_folder);   % Get list of folders corresponding to the dates of experiments
date_folders            = date_folders([date_folders.isdir]);
qremove                 = ismember({date_folders.name},{'.','..'});
date_folders(qremove)   = [];
% initialise 'sdata' struct
sdata	= [];
for a = 1:length(date_folders)
    
    % go to date folder and get file list
    this_date_folder        = date_folders(a).name;
    sdata_files             = dir([sdata_folder filesep this_date_folder]);
    sdata_files             = sdata_files(~[sdata_files.isdir]); % only look at files, not folders
    qremove                 = ismember({sdata_files.name},{'.DS_Store'}); % get rid of annoying  '.DS_store' files
    sdata_files(qremove)    = [];
    file_names              = {sdata_files.name}; % get file names
    
    % loop over data files
    for b = 1:length(sdata_files)
        this_sdata_file = sdata_files(b).name;
        
        full_channels_file      = [sdata_folder filesep this_date_folder filesep this_sdata_file];
        
        % load 'channels' and 'parameters' from saved data
        load(full_channels_file)        
        disp(['Loaded ' full_channels_file '...'])
        
        % Use analyse_channels_function to add PSTH data to 'channels'
        % struct
        channels            = analyse_channels_function(channels,analysisparams);
        
        % get conditions 
        condition_mat    	= cell2mat({channels(1).conditions.timings}');
        
        % store filename and the matrix of conditions in sdata.expt struct
        sdata(a).expt(b).filename                 	= full_channels_file;
        sdata(a).expt(b).condition_mat             	= condition_mat;
        
        %% Start analysis by channel
        
        for c = 1:length(channels)
            
            % store current channel in temporary variable
            temp_channel                            = channels(c).conditions;
            
            % store spontaneous rate of this channel in sdata.expt struct
            sdata(a).expt(b).spont_rate(c)         	= channels(c).spontspikerate;
            
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
            
        end
        sdata(a).expt(b).whisk_profile          = sdata(a).expt(b).whisk_profile(:,:,1:30:end);
    end
end
