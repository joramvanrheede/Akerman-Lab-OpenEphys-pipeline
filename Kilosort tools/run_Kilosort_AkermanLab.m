function run_Kilosort_AkermanLab(data_folder, sort_folder, target_rec_nrs, events_chans, n_channels, varargin)
% function run_Kilosort_AkermanLab(data_folder, sort_folder, target_rec_nrs, events_chans, n_channels, varargin)
% 
% Your one-stop-shop for running Kilosort in the Akerman lab pipeline.
% 
% Takes openepys .continuous data, reformats it as Kilosort-compatible raw
% binary, applies common average referencing, copies relevant Kilosort running
% and config files, links data up with stimulus synchronisation TTL pulses
% on the OpenEphys ADC channels, and runs the Kilosort algorithm in the
% target folder.
% 
% REQUIRED ARGUMENTS:
% DATA_FOLDER: The data folder containing the OpenEphys .continuous files
% SORT_FOLDER: The folder where the spike sorting files will be kept. Will 
% be created if it does not yet exist.
% TARGET_REC_NRS: Recording nrs to be concatenated and sorted
% 
% OPTIONAL:
% 
% EVENTS_CHANS: ADC hannel numbers for synchronisation:
% [Trials Whisk Opto Stim_nr] --> defaults to [1 1 3 2]
% 
% N_CHANNELS: Number of channels of data. Defaults to 32.
% 
% VARARGIN (CURRENTLY NON-FUNCTIONAL):
% various optional input arguments:
% '-CAR': do common average referencing
% '-copy': copy relevant Kilosort files
% '-sync': get stimulus synchronisation data
% '-sort': run kilosort in the folder
% The function will run ALL OF THE ABOVE by default, currently functionality
% for running individual sections is not yet available because of dependence
% of some bits on others.
% 
% 

if nargin < 4
    events_chans    = [1 1 3 2]; % [Trials Whisk Opto Stim_nr]
end

if nargin < 5
    n_channels      = 32; % Number of channels
end

if nargin < 6
    do_CAR              = true;
    get_KS_copyfiles  	= true;
    get_sync_data       = true;
    do_dat_and_concat  	= true;
    run_kilosort        = true;
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

%% convert to raw binary format and concatenate in one step:

% loop to add full file path to the folder names
for c = 1:length(concatenate_folders)
    concatenate_folders{c} = [data_folder filesep concatenate_folders{c}];
end

%% Stimulus synchronisation data (from TTL inputs)
if get_sync_data
    for i = 1:length(concatenate_folders)
        sync_data(i)    = get_stim_sync_data(concatenate_folders{i},events_chans);
    end
    save(fullfile(sort_folder,'sync_data'),'sync_data')
end

full_dat_file_name   = fullfile(sort_folder, dat_file_name);

%% 
if do_dat_and_concat    
    % full file path of the .dat file
    
    % function that takes all individual channel .continuous files and merges
    % them into a single .dat file; it will work on single data sets or
    % concatenate multiple data sets
    concatenate_continuous_as_dat(full_dat_file_name,concatenate_folders, n_channels);
end

%% Optional common average referencing?

if do_CAR
    % do common average referencing:
    disp('Applying common average referencing...')
    applyCARtoDat(full_dat_file_name,n_channels,sort_folder);
end


%% Run kilosort?

if run_kilosort
    disp('Starting Kilosort sorting...')
    current_dir = cd;
    cd(sort_folder)
    Master_spike_sort
    cd(current_dir);
end