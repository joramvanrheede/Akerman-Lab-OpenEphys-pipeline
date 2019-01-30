%% 

metadata_file       = '/Volumes/PS2Akermanlab/Joram/Data/in vivo metadata/Metadata File.xlsx'; % Which metadata file to use?

data_folder         = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data';   %'/Volumes/Akermanlab-1/Joram/In vivo mouse data'; % Where are the data? 

save_folder         = '/Volumes/Akermanlab/Joram/Dysrhythmia';  %/Volumes/Akermanlab-1/Joram/Extracted data'; % Where to save output?

start_date          = '2019_01_24';     % format: 'yyyy_mm_dd'; Process files from this date onwards
end_date            = '2019_01_24';     % format: 'yyyy_mm_dd'; Process files up until this date

process_expts       = {'All'};          % indicate which experiment types to run, e.g.: {'Drive', 'Timing'}, or use {'All'}

get_LFP             = true;             % get LFP traces? This does increase output data size

data_output         = 'new';            % 'new': improved data structure, or 'old': 'channels' style data structure do
%% global analysis parameters

q_spike_detection   = 1;      % do spike detection, or use detected spikes from openephys?

% set this up for different probe types:
channel_order_16ch  = [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]; % 16ch linear silicon A16 probe neuronexus: [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]
channel_order_32ch  = [16 32 1 17 14 30 3 19 8 24 7 29 9 25 15 20 10 23 2 28 6 26 5 21 11 31 4 27 12 22 13 18]; % 32Ch linear A32 probe neuronexus [18 2 31 15 19 3 30 14 17 1 32 16 24 8 25 9 23 7 26 10 22 6 27 11 21 5 28 12 20 4 29 13];

%% running code starts here

% Unpack date strings to screen for which dates should be processed
start_year   	= str2num(start_date(1:4));
start_month  	= str2num(start_date(6:7));
start_day    	= str2num(start_date(9:10));
end_year     	= str2num(end_date(1:4));
end_month   	= str2num(end_date(6:7));
end_day         = str2num(end_date(9:10));

% Read excel file containing metadata
[num, text, metadata] = xlsread(metadata_file,1);
headers         = metadata(1,:); % top line contains the headers for the columns
metadata        = metadata(2:end,:); % the rest is data

% Find relevant columns in excel file; note that by using strcmpi the column
% order of the excel file is not fixed and it is easy to add extra columns,
% however the script will need updating if you change the text of a header
date_col        = find(strcmpi('Date',headers)); % Date
ID_col          = find(strcmpi('Animal ID',headers)); % Animal ID (Sanger system preferred, to trace animal / electroporation notes)
type_col        = find(strcmpi('Animal type',headers)); % EG 'IUE ChR2-YFP-ON' or 'CAMKII-GCAMP Transgenic'
age_col         = find(strcmpi('Age (w)',headers)); % Age of animal in weeks

probe_col       = find(strcmpi('Probe type',headers));
pen_col         = find(strcmpi('Penetration nr',headers)); % Penetration number (for this animal)
caud_col        = find(strcmpi('Caudal',headers)); % Caudal position (mm)
lat_col         = find(strcmpi('Lateral',headers)); % Lateral position (mm)
depth_col       = find(strcmpi('Depth',headers)); % Depth of penetration (mm)
angle_col       = find(strcmpi('Angle',headers)); % Angle of penetration (0deg = vertical)

file_col        = find(strcmpi('File nr',headers)); % file number
LEDpw_col       = find(strcmpi('LED power (mA)',headers)); % LED power
expt_col        = find(strcmpi('Experiment',headers)); % Experiment type
stim_col        = find(strcmpi('Stimulator',headers)); % What stimulator is being used
whisk_col       = find(strcmpi('Whiskers',headers)); % Which whiskers are being stimulated?

trialTTL_col    = find(strcmpi('TrialTTL',headers)); % Which TTL channel has info about trials
whiskTTL_col    = find(strcmpi('WhiskTTL',headers)); % Which TTL channel has info about whisker stimulation
LEDTTL_col      = find(strcmpi('LEDTTL',headers)); % Which TTL channel has info about LED on/off
switchTTL_col   = find(strcmpi('Stim_switch_TTL',headers)); % Which TTL channel is used for switching stimulators?

override_col    = find(strcmpi('Override conditions',headers)); % If some event TTLs are broken, you can use override col to impose a nr. of conditions. only works if condition order is not randomised.
n_cond_col      = find(strcmpi('Number of conditions',headers)); % If overriding the TTL conditions, 

spkthresh_col   = find(strcmpi('Spike threshold',headers)); % Threshold for spike detection
smoothw_col     = find(strcmpi('Smoothing',headers)); % Smoothing of the signal (in number of samples) before applying threshold for spike detection
morse_col       = find(strcmpi('Morse',headers)); % Does this experiment use a 'morse' signal to signal the start and end of a trial? 
chanmap_col     = find(strcmpi('Channelmap',headers)); % Does this experiment use the channelmap module from openephys? 
dig_events_col  = find(strcmpi('Digital events',headers)); % Are events registered using Digital inputs (and auto event detection) on the acquisition system, or are they analog inputs that need manual extraction of events?
prefix_col      = find(strcmpi('File prefix',headers)); % Prefix for data files
whisk_res_col   = find(strcmpi('Whisker timing minimal increment',headers)); % Minimal increment for automatic extraction of whisker trigger onset. For old data (2017), recommended to keep this at 250 ms
LED_res_col     = find(strcmpi('LED timing minimal increment',headers)); % Minimal increment for automatic extraction of LED timing conditions
spontwin_col    = find(strcmpi('Spontwin max',headers)); % Window at the start of each trial but before any whisk or LED stim, to use as a measure of spontaneous activity

for a = 1:size(metadata,1)
    
    % find experiment date; see if this date is within the processing range
    this_date       = metadata{a,date_col};
    
    this_year    	= str2num(this_date(1:4));
    this_month   	= str2num(this_date(6:7));
    this_day        = str2num(this_date(9:10));
    
    if this_year < start_year
        continue
    elseif this_year == start_year && this_month < start_month
        continue
    elseif this_year == start_year && this_month == start_month && this_day < start_day
        continue
    elseif this_year > end_year
        continue
    elseif this_year == end_year && this_month > end_month
        continue
    elseif this_year == end_year && this_month == end_month && this_day > end_day
        continue
    end
    
    % find this experiment type; see if it needs processing
    this_expt       = metadata{a,expt_col};
    
    if ~any(strcmp(this_expt,process_expts)) && ~any(strcmpi('All',process_expts))
        continue
    end
    
    % Extract other relevant data
    this_ID         = metadata{a,ID_col};
    this_type       = metadata{a,type_col};
    this_probe      = metadata{a,probe_col};
    this_age        = metadata{a,age_col};
    this_pen        = metadata{a,pen_col};
    this_caud       = metadata{a,caud_col};
    this_lat        = metadata{a,lat_col};
    this_depth      = metadata{a,depth_col};
    this_angle      = metadata{a,angle_col};
    this_file       = metadata{a,file_col};
    this_LEDpw      = metadata{a,LEDpw_col};
    
    this_stim       = metadata{a,stim_col};
    this_whisk      = metadata{a,whisk_col};
    
    switch this_probe
        case '16ch'
            parameters.get_channels         = channel_order_16ch;
        case '32ch'
            parameters.get_channels         = channel_order_32ch;
    end
    
    % all variables that are only relevant in the 'parameters' struct for
    % extract_ephys_data_function can go directly into the struct
    
    parameters.digital_events     	= metadata{a,dig_events_col};
    parameters.data_prefix         	= metadata{a,prefix_col};
    parameters.channelmap         	= metadata{a,chanmap_col};
    parameters.get_LFP              = get_LFP;
    parameters.LED_power            = this_LEDpw;
    
    parameters.spike_thresh        	= metadata{a,spkthresh_col};
    parameters.spike_smoothwin    	= metadata{a,smoothw_col};
    parameters.LED_conditions_res	= metadata{a,LED_res_col};
    parameters.whisk_conditions_res = metadata{a,whisk_res_col};
    parameters.spontwin          	= [0 metadata{a,spontwin_col}]; % start of trial, 0 to X ms
    parameters.morse                = metadata{a,morse_col};
    
    parameters.trial_channel      	= metadata{a,trialTTL_col};         % Which input channel has the episode TTL
    parameters.whisk_channel     	= metadata{a,whiskTTL_col};         % Which input channel has the piezo / whisk TTL
    parameters.LED_channel          = metadata{a,LEDTTL_col};           % Which input channel has the LED TTL    
    parameters.stim_switch_channel  = metadata{a,switchTTL_col};        % Which input channel has the stimulator switch TTL
    parameters.experiment_type      = metadata{a,expt_col};             % What type of experiment is this?
    
    parameters.override_conds       = metadata{a,override_col};
    parameters.n_conds              = metadata{a,n_cond_col};
    
    parameters.animal_type          = metadata(a,type_col);
    parameters.target_whisker       = this_whisk;
    
    parameters.data_output          = data_output;
    
    %% time to extract the data
    this_data_folder                = [data_folder filesep this_date];
    
    % Using a 'try' and 'catch' structure so that data processing will
    % continue even if there is an error with one particular data file
    try 
        % data extraction function (where it all really happens)
        ephys_data              	= extract_ephys_data_function(this_data_folder, this_file, q_spike_detection, parameters);
    catch
        % Report any errors in processing
        warning(['Error processing file ' this_data_folder ' #' num2str(this_file) ])
        disp('Error message: ')
        errhandle = lasterror;
        disp(errhandle.message)
        for a = 1:length(errhandle.stack)
            disp('in')
            disp(errhandle.stack(a))
        end
        continue
    end
    
    % organise output in folders by experiment name (expt) - timing, drive,
    % velocity, frequency, etc.
    this_save_file_name             = [this_date '-' num2str(this_file) '-' this_expt];
    this_save_folder                = [save_folder filesep this_expt filesep this_date];
    
    % If folder of that name does not exist, create it
    if ~isdir(this_save_folder)
        mkdir(this_save_folder)
    end
    
    % save channels struct in a folder 
    save([this_save_folder filesep this_save_file_name],'ephys_data')
    
end
