function run_Kilosort_AkermanLab(data_folder, sort_folder, target_rec_nrs, metadata_file, varargin)
% function run_Kilosort_AkermanLab(DATA_FOLDER, SORT_FOLDER, TARGET_REC_NRS, METADATA_FILE, OPTIONS)
% 
% Your one-stop-shop for running Kilosort in the Akerman lab pipeline.
% 
% Takes openepys .continuous data, reformats it as Kilosort-compatible raw
% binary, applies common average referencing, copies relevant Kilosort running
% and config files, retrieves relevant metadata and stimulus synchronisation 
% TTL pulses on the OpenEphys ADC channels, and runs the Kilosort algorithm 
% in the target folder.
% 
% REQUIRED ARGUMENTS:
% DATA_FOLDER: The data folder containing the OpenEphys .continuous files.
% SORT_FOLDER: The folder where the spike sorting files will be kept. Will 
% be created if it does not yet exist.
% TARGET_REC_NRS: Recording nrs to be concatenated and sorted.
% METADATA_FILE: Full path to the metadata file with the info for the relevant
% recordings.
% 
% Optional:
% 
% OPTIONS:
% various optional input arguments:
% [32] / [16]: sets number of channels to specified number
% '-cat': do concatenation of openephys data
% '-CAR': do common average referencing
% '-copy': copy relevant Kilosort files to target folder
% '-sync': get stimulus synchronisation data
% '-sort': run kilosort in the folder
% The function will run ALL OF THE ABOVE by default, currently functionality
% for running individual sections is not yet available because of dependence
% of some bits on others.
% 

if nargin < 5
    do_CAR              = true;
    get_KS_copyfiles  	= true;
    get_sync_data       = true;
    do_dat_and_concat  	= true;
    run_kilosort        = true;
    n_channels          = 32;
else
    do_CAR              = false;
    get_KS_copyfiles  	= false;
    get_sync_data       = false;
    do_dat_and_concat  	= false;
    run_kilosort        = false;
 
    if any(strcmpi('-CAR',varargin))
        do_CAR              = true;
    end
    if any(strcmpi('-cat',varargin))
        do_dat_and_concat  	= true;
    end
    if any(strcmpi('-sync',varargin))
        get_sync_data       = true;
    end
    if any(strcmpi('-copy',varargin))
        get_KS_copyfiles  	= true;
    end
    if any(strcmpi('-sort',varargin))
        run_kilosort        = true;
    end
    numeric_ind    = cellfun(@isnumeric,varargin);
    if any(numeric_ind)
        n_channels = varargin{numeric_ind};
    end
%     if any(~ismember(varargin,{'-CAR' '-cat' '-sync' '-copy' '-sort' '' '' ''}))
%         error('Unknown option provided. Valid options are: ''-copy''  ''-sync''  ''-cat''  ''-CAR''  ''-sort''')
%     end
end

% Hardcoded concatenated data file name
dat_file_name       = 'concatenated_data.dat';

% Use current directory in repository to reconstruct file path of Kilosort copyfiles
function_folder  	= fileparts(mfilename('fullpath'));
root_folder         = fileparts(function_folder);
KS_copy_file_folder = [root_folder filesep 'Kilosort copyfiles/32Ch_standard'];

%% Copy Kilosoft copy files to new sorting folder

% Always make sort folder
if ~isdir(sort_folder)
	mkdir(sort_folder);
end

%% If requested, copy relevant Kilosort files over
if get_KS_copyfiles
    Kilosort_files        = dir([KS_copy_file_folder]);
    Kilosort_files        = {Kilosort_files.name};
    Kilosort_files(ismember(Kilosort_files,{'.' '..' '.DS_Store'}))     = []; % get rid of '.' and '..'
    
    for i = 1:length(Kilosort_files)
        [success, message, messageid] = copyfile([KS_copy_file_folder filesep Kilosort_files{i}],sort_folder);
    end
end

%% Find the data folders for all the requested file numbers

protocol_folders        = dir([data_folder]);
protocol_folders        = {protocol_folders.name};
protocol_folders(ismember(protocol_folders,{'.' '..' '.DS_Store'}))     = []; % get rid of '.' and '..'

% find the part of the file name where the protocol nr is mentioned (it is the number at the end, after the underscore)
rec_number_start_ind    = cell2mat(regexp(protocol_folders,repmat({'\d+$'},size(protocol_folders))));
rec_number_end_ind      = cell2mat(regexp(protocol_folders,repmat({'\d$'},size(protocol_folders))));

protocol_nrs            = NaN(size(protocol_folders));
for a = 1:length(protocol_folders)
    this_folder         = protocol_folders{a};
    
    protocol_nrs(a)     = str2num(this_folder(rec_number_start_ind(a):rec_number_end_ind(a)));
end

q_protocols             = ismember(protocol_nrs,target_rec_nrs);

concatenate_folders     = protocol_folders(q_protocols);
rec_nrs                 = protocol_nrs(q_protocols);

% loop to add full file path to the folder names
for c = 1:length(concatenate_folders)
    concatenate_folders{c} = [data_folder filesep concatenate_folders{c}];
end

%% Stimulus synchronisation data (from TTL inputs and Metadata file)
if get_sync_data
    
    [~, date_string]        = fileparts(data_folder);
    
    [metadata, headers]     = load_metadata(metadata_file);
    
    for i = 1:length(concatenate_folders)
        
        metadata_info           = read_metadata(metadata, headers, date_string, rec_nrs(i));
        
        expt_type               = metadata_info.experiment_type;
        
        % set moving minimum window for baseline fixing
        switch expt_type
            case 'LED_power'
                baseline_moving_win = 12;
                LED_conditions_res  = 5;
            case 'Laser_pulse'
                baseline_moving_win = 12;
                LED_conditions_res  = 5;
            case 'Timing'
                baseline_moving_win = 12;
                LED_conditions_res  = 5;
            case 'Drive'
                baseline_moving_win = 1020;
                LED_conditions_res  = 100;
            case 'Drive_early'
                baseline_moving_win = 1020;
                LED_conditions_res  = 100;
            case 'Ramp'
                baseline_moving_win = 10020;
                LED_conditions_res  = 1000;
            case 'Short_ramp'
                baseline_moving_win = 10020;
                LED_conditions_res  = 100;
            case 'Long_opto'
                baseline_moving_win = 10020;
                LED_conditions_res  = 100;
            case 'Dual_delay'
                baseline_moving_win = 12;
                LED_conditions_res  = 5;
            case 'Multiwhisk_ramp'
                baseline_moving_win = 10020;
                LED_conditions_res  = 100;
            otherwise
                baseline_moving_win = 10020;
        end
        
        sync_data(i)            = get_stim_sync_data(concatenate_folders{i},metadata_info,0,[],baseline_moving_win);
        
    end
    save(fullfile(sort_folder,'sync_data'),'sync_data')
end


%% convert to raw binary format and concatenate in one step:
full_dat_file_name   = fullfile(sort_folder, dat_file_name);

if do_dat_and_concat
    % full file path of the .dat file
    
    % function that takes all individual channel .continuous files and merges
    % them into a single .dat file; it will work on single data sets or
    % concatenate multiple data sets
    concatenate_continuous_as_dat(full_dat_file_name,concatenate_folders, n_channels);
end

%% Common average referencing (not optional)

if do_CAR
    % do common average referencing:
    disp('Applying common average referencing...')
    applyCARtoDat(full_dat_file_name,n_channels,sort_folder);
end


%% Run kilosort?

if run_kilosort
    disp('Starting Kilosort sorting...')
    current_dir = cd;
    
    % Move current directory to relevant sort folder
    cd(sort_folder)
    Master_spike_sort
    cd(current_dir);
end