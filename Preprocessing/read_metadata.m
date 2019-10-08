function metadata_info = read_metadata(metadata,headers,date_string,rec_number)
% METADATA_INFO = read_metadata(METADATA,HEADERS,DATE_STRING,REC_NUMBER)
%
% Extract METADATA_INFO from METADATA and HEADERS (obtained with LOAD_METADATA 
% function).
% 
% INPUTS:
% 
% METADATA: Metadata loaded using LOAD_METADATA
%
% HEADERS: Headers for METADATA, obtained from LOAD_METADATA
%
% DATE_STRING: Target date in format 'yyyy_mm_dd' (as entered in metadata file!)
%
% REC_NUMBER: number indicating which recording number from the given date
% you want to get metadata for.
%
% OUTPUT:
% 
% METADATA_INFO: Information from the metadata file (METADATA_FILE) regarding
% the target recording number (REC_NUMBER) from the target date (DATE_STRING).
%
%



%% global analysis metadata_info

% set this up for different probe types:
channel_order_16ch  = [13 29 4 20 15 31 3 19 16 32 1 17 2 18 14 30]; %% Corrected 08/02/2019 % 16ch linear silicon A16 probe neuronexus:
channel_order_32ch  = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31]; % Corrected 08/02/2019

%% running code starts here


% Find relevant columns in excel file; note that by using strcmpi the column
% order of the excel file is not fixed and it is easy to add extra columns,
% however the script will need updating if you change the text of a header
date_col        = find(strcmpi('Date',headers)); % Date
ID_col          = find(strcmpi('Animal ID',headers)); % Animal ID (Sanger system preferred, to trace animal / electroporation notes)
type_col        = find(strcmpi('Animal type',headers)); % EG 'IUE ChR2-YFP-ON' or 'CAMKII-GCAMP Transgenic'
age_col         = find(strcmpi('Age (w)',headers)); % Age of animal in weeks
sex_col         = find(strcmpi('Sex',headers)); % Sex of animal

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

trials_from_whisk_col   = find(strcmpi('Trials from whisk',headers));
whisk_buffer_col        = find(strcmpi('Whisk buffer',headers));

% find experiment date; see if this date is within the processing range
all_dates       = metadata(:,date_col);
all_rec_nrs     = [metadata{:,file_col}]';

q_date          = strcmp(date_string,all_dates);
q_rec_nr        = all_rec_nrs == rec_number;

q_select        = q_date & q_rec_nr;

if sum(q_select) == 0 
    warning(['Target recording info not found for date ' date_string ' rec nr ' num2str(rec_number) ' in metadata file ' metadata_file '.'])
    metadata_info = [];
    return
elseif sum(q_select) > 1
    warning(['Multiple matches found for date ' date_string ' rec nr ' num2str(rec_number) ' in metadata file ' metadata_file '.'])
    metadata_info = [];
    return
end

% Extract other relevant data
metadata_info.experiment_date       = date_string;
metadata_info.animal_ID             = metadata{q_select,ID_col};
metadata_info.animal_type           = metadata{q_select,type_col};
metadata_info.probe_type            = metadata{q_select,probe_col};
metadata_info.animal_age            = metadata{q_select,age_col};
if ~isempty(sex_col)
    metadata_info.animal_sex      	= metadata{q_select,sex_col};
else
    metadata_info.animal_sex      	= 'U'; % U for 'Unknown'
end
metadata_info.implantation_nr       = metadata{q_select,pen_col};
metadata_info.caudal_coordinate     = metadata{q_select,caud_col};
metadata_info.lateral_coordinate    = metadata{q_select,lat_col};
metadata_info.depth_coordinate      = metadata{q_select,depth_col};
metadata_info.implantation_angle  	= metadata{q_select,angle_col};
metadata_info.recording_nr          = metadata{q_select,file_col};

metadata_info.stim_nr               = metadata{q_select,stim_col};

switch metadata_info.probe_type
    case '16ch'
        metadata_info.get_channels	= channel_order_16ch;
    case '32ch'
        metadata_info.get_channels	= channel_order_32ch;
end

% all variables that are only relevant in the 'metadata_info' struct for
% extract_ephys_data_function can go directly into the struct

metadata_info.digital_events     	= metadata{q_select,dig_events_col};
metadata_info.data_prefix         	= metadata{q_select,prefix_col};
metadata_info.channelmap         	= metadata{q_select,chanmap_col};
metadata_info.LED_power             = metadata{q_select,LEDpw_col};

metadata_info.spike_thresh        	= metadata{q_select,spkthresh_col};
metadata_info.spike_smoothwin    	= metadata{q_select,smoothw_col};
metadata_info.LED_conditions_res	= metadata{q_select,LED_res_col};
metadata_info.whisk_conditions_res  = metadata{q_select,whisk_res_col};
metadata_info.spontwin          	= [0 metadata{q_select,spontwin_col}]; % start of trial, 0 to X ms
metadata_info.morse                 = metadata{q_select,morse_col};

metadata_info.trial_channel      	= metadata{q_select,trialTTL_col};         % Which input channel has the episode TTL
metadata_info.whisk_channel     	= metadata{q_select,whiskTTL_col};         % Which input channel has the piezo / whisk TTL
metadata_info.LED_channel           = metadata{q_select,LEDTTL_col};           % Which input channel has the LED TTL
metadata_info.stim_switch_channel   = metadata{q_select,switchTTL_col};        % Which input channel has the stimulator switch TTL
metadata_info.experiment_type       = metadata{q_select,expt_col};             % What type of experiment is this?

metadata_info.override_conds        = metadata{q_select,override_col};
metadata_info.n_conds               = metadata{q_select,n_cond_col};

% New additions to metadata file - if not present, default to 0 / false
if ~isempty(trials_from_whisk_col)
    metadata_info.trials_from_whisk     = metadata{q_select,trials_from_whisk_col};
    metadata_info.whisk_buffer          = metadata{q_select,whisk_buffer_col};
else
    metadata_info.trials_from_whisk     = 0;
    metadata_info.whisk_buffer          = 0;
end

metadata_info.animal_type           = metadata(q_select,type_col);
metadata_info.target_whisker        = metadata(q_select,whisk_col);


