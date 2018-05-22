
function oe2mda(oe_data_loc,save_loc,channels,detect_threshold)

% detect_threshold            = 5;
% channels                    = [1:4 13:20 29:32];
% save_loc                    = '/Users/Joram/Mountainsort_test/mda';
% file_loc                    = '/Users/Joram/Documents/In vivo data/2018_04_24/ChR2-YFP-ON_2018-04-24_13-32-29_1';

file_prefix               	= '100_CH'; % hardcoded in as this doesn't change in our data set
file_ext                    = '.continuous'; % as above

% Loop over channels
for a = 1:length(channels)
    
    % find which channel nr we're on
    this_channel                = channels(a);
    chan_nr                     = num2str(this_channel);
    
    % make file name string of openephys .continuous file to be read
    read_file                   = [oe_data_loc filesep file_prefix chan_nr file_ext];
    
    % load openephys data file
    [data, timestamps, info]    = load_open_ephys_data(read_file);
    
    % make file name string
    write_file                  = [save_loc filesep file_prefix chan_nr '.mda'];
    
    % make directory for mda files if it does not yet exist
    if ~isdir(save_loc)
        mkdir(save_loc)
    end
    
    % write file as mda (note that data has to go in as a row)
    disp('writing mda file')
    writemda(data',write_file)
    
    % Run mountainlab processes using a string of system commands; all as
    % one long command to system because otherwise matlab executes each
    % statement in isolation rather than letting them build on one another
    % (i.e. it forgets the path and virtual environment).
    % Any relevant .bash_profile statements need to be included too.
    system(['export PATH=/Applications/mountainlab-js/bin/:$PATH ; ' ...    % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/ephys-viz/bin/:$PATH ; ' ...   % Duplicate from ~/.bash_profile
        'export PATH=/usr/local/bin/:$PATH ; ' ...                          % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/:$PATH ; ' ...                 % Duplicate from ~/.bash_profile
        'export PATH=/usr/local/lib/python3.6/site-packages/:$PATH ; ' ...  % Duplicate from ~/.bash_profile
        'alias python=''python3'' ; ' ...                                   % Duplicate from ~/.bash_profile
        'source ~/ml_venv/bin/activate ; ' ...                              % Start previously created ml_venv virtual environment for running mountainlab
        'cd ' save_loc ' ; ' ...                                            % cd to working directory
        
        % Use ephys.bandpass_filter to filter data
        'ml-run-process ephys.bandpass_filter --inputs timeseries:' write_file ' --outputs timeseries_out:./filt_' chan_nr '.mda.prv --parameters samplerate:30000 freq_min:300 freq_max:6000 ; ' ...
        % Use ephys.whiten to whiten the data
        'ml-run-process ephys.whiten --inputs timeseries:./filt_' chan_nr '.mda.prv --outputs timeseries_out:./pre' chan_nr '.mda.prv ;' ...
        % Use ms4alg.sort for spike detection, clustering and sorting
        'ml-run-process ms4alg.sort --inputs timeseries:./pre' chan_nr '.mda.prv  --outputs firings_out:./firings' chan_nr '.mda.prv --parameters adjacency_radius:0 detect_sign:-1 detect_threshold:' num2str(detect_threshold)]);
        
end


