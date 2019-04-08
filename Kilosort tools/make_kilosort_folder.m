% make_kilosort_folder

% This process should do the following:
% 1) convert target data set to binary format - check
% 2) concatenate the data for the recordings in the data set so they can
% all be sorted in one go and units can be compared across these recordings
% - check
% 3) Optional common average referencing - NEEDS WORK
% 4) Put the converted, concatenated binary data file in a new folder that
% will be the root directory for running kilosort - check
% 5) Copy the relevant kilosort files into the root folder for running
% Kilosort as well - check

% KILOSORT COPYFILES
KS_copy_file_folder	= '/Users/Joram/Dropbox/Akerman Postdoc/Code/OpenEphys_pipeline/Kilosort copyfiles/32Ch_standard';

% DATA IN
data_folder         = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data';
target_date         = '2019_03_09';
target_file_nrs     = [1 3 4 5 6];

% DATA OUT
sorting_folder      = '/Volumes/Akermanlab/Joram/Spike_sorting'; %
dat_file_name       = 'concatenated_data.dat';

% OPTIONS
do_CAR              = false;

% currently non-functional:
get_KS_copyfiles  	= false;
get_sync_data       = true;
do_dat_and_concat  	= false;

events_chans        = [1 1 3 2]; % sync data ADC channels order - [Trials Whisk Opto Stim_nr]

n_channels          = 32; % number of channels of data

%% Copy Kilosoft copy files to new sorting folder

% Create target folder 
folder_postfix  = '-';
for b = 1:length(target_file_nrs)
    folder_postfix  = [folder_postfix '_' num2str(target_file_nrs(b))];
end

sort_folder     = [sorting_folder filesep target_date folder_postfix];

if ~isdir(sort_folder)
	mkdir(sort_folder);
end

%%
if get_KS_copyfiles
    Kilosort_files        = dir([KS_copy_file_folder]);
    Kilosort_files        = {Kilosort_files.name};
    Kilosort_files(ismember(Kilosort_files,{'.' '..' '.DS_Store'}))     = []; % get rid of '.' and '..'
    
    for i = 1:length(Kilosort_files)
        [success, message, messageid] = copyfile([KS_copy_file_folder filesep Kilosort_files{i}],sort_folder);
    end
end

%% Find the data folders

data_date_folder        = [data_folder filesep target_date filesep];

protocol_folders        = dir([data_date_folder]);
protocol_folders        = {protocol_folders.name};
protocol_folders(ismember(protocol_folders,{'.' '..'}))     = []; % get rid of '.' and '..'

% find the part of the file name where the protocol nr is mentioned (it is the number at the end, after the underscore)
rec_number_start_ind    = cell2mat(regexp(protocol_folders,repmat({'\d+$'},size(protocol_folders))));
rec_number_end_ind      = cell2mat(regexp(protocol_folders,repmat({'\d$'},size(protocol_folders))));

protocol_nrs            = NaN(size(protocol_folders));
for a = 1:length(protocol_folders)
    this_folder         = protocol_folders{a};
    
    protocol_nrs(a)     = str2num(this_folder(rec_number_start_ind(a):rec_number_end_ind(a)));
end

q_protocols             = ismember(protocol_nrs,target_file_nrs);

concatenate_folders     = protocol_folders(q_protocols);

%% convert to binary and concatenate in one step:

% loop to add full file path to the folder names
for c = 1:length(concatenate_folders)
    concatenate_folders{c} = [data_date_folder filesep concatenate_folders{c}];
end

if get_sync_data
    clear sync_data
    for i = 1:length(concatenate_folders)
        sync_data(i)    = get_stim_sync_data(concatenate_folders{i},events_chans);
    end
    save(fullfile(sort_folder,'sync_data'),'sync_data')
end
    
if do_dat_and_concat    
    % full file path of the .dat file
    dat_file_name   = [sort_folder filesep dat_file_name];
    
    % function that takes all individual channel .continuous files and merges
    % them into a single .dat file; it will work on single data sets or
    % concatenate multiple data sets
    concatenate_continuous_as_dat(dat_file_name,concatenate_folders, n_channels);
end

%% Optional common average referencing?

if do_CAR
    % do common average referencing:
    disp('Applying common average referencing')
    applyCARtoDat([dat_file_name],n_channels,sort_folder);
end

%% 