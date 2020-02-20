function preprocess_multiunit(metadata_file, data_folder, save_folder, start_date, end_date, process_expts, baseline_moving_window)
% function preprocess_multiunit(METADATA_FILE, DATA_FOLDER, SAVE_FOLDER, START_DATE, END_DATE, PROCESS_EXPTS,BASELINE_MOVING_WINDOW)
% OR
% function preprocess_multiunit(METADATA_FILE, DATA_FOLDER, SAVE_FOLDER, START_DATE, FILE_NRS,BASELINE_MOVING_WINDOW)
% 
% First step in the OpenEphys analysis pipeline in the Akerman lab. This will
% read in the raw data in DATA_FOLDER based on the info in metadata spreadsheet
% METADATA_FILE and distribute per-experiment files of spike times and LFP traces
% relative to trial onset organised by condition in SAVE_FOLDER. START_DATE sets the
% date from which to preprocess data - if only START_DATE is specified the function
% will only preprocess data from this date, if an END_DATE is specified the function
% will preprocess data from START_DATE up to and including END_DATE.
% 
% If FILE_NRS are specified instead of END_DATE, the function will preprocess only
% the specified file numbers from START_DATE (and preprocess for this date only).
% 
% You can preprocess only target experiment types using PROCESS_EXPTS.
% 
% REQUIRED INPUTS:
% 
% METADATA_FILE: Full file path of metadata file containing metadata for this experiment.
% See /Metadata/metadata_template.xlsx for an example.
%
% DATA_FOLDER: Full path to the folder containing the data (the parent folder
% that contains the folders for individual dates).
%
% SAVE_FOLDER: The folder where you want the preprocessed data to be saved.
%
% START_DATE: Preprocess experiments from this date onwards - format 'yyyy_mm_dd'
% e.g. '2019_01_01'.
% 
% OPTIONAL INPUTS:
% 
% END_DATE: Preprocess experiments until this date. Default is the same as START_DATE
% so the function will only process data for one target day.
% OR
% FILE_NRS: Preprocess only the target file numbers on the specified START_DATE
% (if FILE_NRS is specified the function will only preprocess a single date).
% 
% PROCESS_EXPTS: Enter experiment types that you want to preprocess as a list of strings
% in a cell, e.g. {'Drive' 'Timing' 'Laser_pulse'}; defaults to {'All'}, i.e. all experiment
% types in the date range will be preprocessed.
% 
% BASELINE_MOVING_WINDOW: Moving minimum window for fixing baseline in milliseconds. 
% Baseline fixing is disabled by default, but is enabled when this input is provided. 
% For the moving minimum window, the best value is 'a few milliseconds longer than 
% your longest LED stimulus', and in general the baseline fix will give better results
% for recordings with shorter LED stimuli.
%
% TO DO - Allow trials from whisk & whisk_buffer
% 

if nargin < 5
    end_date            = start_date;
end

if isnumeric(end_date)
    file_nrs    = end_date;
    end_date    = start_date;
    if exist('process_expts','var')
        baseline_moving_window  = process_expts;
        process_expts   = 'All';
    end
end

if nargin < 6
    process_expts       = 'All';
end

if exist('baseline_moving_window','var')
    q_fix_baseline          = true;
else
    q_fix_baseline          = false;
    baseline_moving_window  = NaN;
end

    
if nargin < 8
    %% Change these only in exceptional cases
    get_LFP             = true;             % get LFP traces? This does increase output data size
    trials_from_whisk   = false;          	% discard trial information from ADC channels and determine trials based on whisker instead?
    whisk_buffer        = 3;              % 0.0625 	% if using whisk stim to divide recording into trials (above), trials start whisk_buffer (in seconds) before the whisker stim onset, and end 2*whisk buffer after whisker stim ONSET
end

%% These things shouldn't really change
do_CAR              = true;             % Do common average referencing?
save_sync_chans     = true;
sync_chans_res      = 1000;

%% HARDCODED CHANNEL MAPS!!!
channel_order_16ch  = [13 29 4 20 15 31 3 19 16 32 1 17 2 18 14 30]; %% Corrected 08/02/2019 % 16ch linear silicon A16 probe neuronexus: 
channel_order_32ch  = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31]; % Corrected 08/02/2019


%% running code starts here

tic

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
animal_nr_col   = find(strcmpi('Animal number',headers));
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

spkthresh_col   = find(strcmpi('Spike threshold',headers)); % Threshold for spike detection
smoothw_col     = find(strcmpi('Smoothing',headers)); % Smoothing of the signal (in number of samples) before applying threshold for spike detection
prefix_col      = find(strcmpi('File prefix',headers)); % Prefix for data files
whisk_res_col   = find(strcmpi('Whisker timing minimal increment',headers)); % Minimal increment for automatic extraction of whisker trigger onset. For old data (2017), recommended to keep this at 250 ms
LED_res_col     = find(strcmpi('LED timing minimal increment',headers)); % Minimal increment for automatic extraction of LED timing conditions

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
    
    % find experiment number; see if it needs processing
    this_file_nr    = metadata{a,file_col};
    
    if exist('file_nrs','var') && ~ismember(this_file_nr, file_nrs)
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

    parameters.data_prefix         	= metadata{a,prefix_col};
    parameters.get_LFP              = get_LFP;
    parameters.LED_power            = this_LEDpw;
    
    parameters.spike_thresh        	= metadata{a,spkthresh_col};
    parameters.spike_smoothwin    	= metadata{a,smoothw_col};
    parameters.LED_conditions_res	= metadata{a,LED_res_col};
    parameters.whisk_conditions_res = metadata{a,whisk_res_col};
    
    parameters.trial_channel      	= metadata{a,trialTTL_col};         % Which input channel has the episode TTL
    parameters.whisk_channel     	= metadata{a,whiskTTL_col};         % Which input channel has the piezo / whisk TTL
    parameters.LED_channel          = metadata{a,LEDTTL_col};           % Which input channel has the LED TTL
    parameters.stim_switch_channel  = metadata{a,switchTTL_col};        % Which input channel has the stimulator switch TTL
    parameters.experiment_type      = metadata{a,expt_col};             % What type of experiment is this?
    
    parameters.animal_type          = metadata(a,type_col);
    parameters.animal_nr            = metadata{a,animal_nr_col};
    parameters.target_whisker       = this_whisk;
    
    parameters.trials_from_whisk    = trials_from_whisk;
    parameters.whisk_buffer         = whisk_buffer;
    parameters.do_CAR               = do_CAR;
    parameters.save_sync_chans      = save_sync_chans;
    parameters.sync_chans_res       = sync_chans_res;
    
    % set moving minimum window for baseline fixing
    switch parameters.experiment_type
        case 'LED_power'
            baseline_moving_win = 12;
        case 'Laser_pulse'
            baseline_moving_win = 12;
        case 'Timing'
            baseline_moving_win = 12;
        case 'Drive'
            baseline_moving_win = 1020;
        case 'Drive_early'
            baseline_moving_win = 1020;
        case 'Ramp'
            baseline_moving_win = 10020;
        case 'Short_ramp'
            baseline_moving_win = 10020;
        case 'Long_opto'
            baseline_moving_win = 10020;
        case 'Dual_delay'
            baseline_moving_win = 12;
        case 'Multiwhisk_ramp'
            baseline_moving_win = 10020;
        otherwise
            baseline_moving_win = 10020;
    end
    
    parameters.baseline_moving_window   = baseline_moving_win;
    
    %% time to extract the data
    this_data_folder                = [data_folder filesep this_date];
    
    % Using a 'try' and 'catch' structure so that data processing will
    % continue even if there is an error with one particular data file
    try 
        % data extraction function (where it all really happens)
        ephys_data              	= extract_ephys_data(this_data_folder, this_file, parameters);
        toc
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
