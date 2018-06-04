

channels            = [1:4 13:20 29:32];
detect_threshold    = 4;

range               = 1:(30000*600);

oe_data_loc         = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data/2018_04_05/ChR2-YFP-ON_2018-04-05_16-55-04_4';
mda_loc             = [oe_data_loc filesep 'mda'];
geom_file           = [cd filesep '16ch_site_geometry.csv'];

file_prefix         = '100_CH';

qconvert            = 0; % (re)do conversion to mda?
qsort               = 1; % (re)do sorting?

view_channel        = 1; % 2018/04/05: 16,29,13

qinspect            = 1;
qmountainview       = 0;

%% File conversion from openephys .continuous to .mda
if qconvert
    oe2mda_multi(oe_data_loc,mda_loc,channels,range) % convert data for this channel to mda
    copyfile(geom_file,[mda_loc filesep 'geom.csv']) % copy geometry file to mda directory
end

%% Spike sorting
if qsort

        
        % make file name string
        write_file  	= [mda_loc filesep 'All_16_chan.mda'];
        
        disp(' ')
        disp(' ')
        disp(['SORTING DATA FOR CHANNELS'])
        disp(' ')
        disp(' ')
        
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
            'cd ' mda_loc ' ; ' ...                                            % cd to working directory
            ...
            ... % Use ephys.bandpass_filter to filter data
            'ml-run-process ephys.bandpass_filter --inputs timeseries:' write_file ' --outputs timeseries_out:' mda_loc '/filt' '_ALL' '.mda.prv --parameters samplerate:30000 freq_min:300 freq_max:6000 ; ' ...
            ... % Use ephys.whiten to whiten the data
            'ml-run-process ephys.whiten --inputs timeseries:./filt_ALL.mda.prv --outputs timeseries_out:' mda_loc '/pre' '_ALL' '.mda.prv ;' ...
            ... % Use ms4alg.sort for spike detection, clustering and sorting
            'ml-run-process ms4alg.sort --inputs timeseries:./All_16_chan.mda geom:' mda_loc '/geom.csv  --outputs firings_out:' mda_loc '/firings' '_ALL' '.mda.prv --parameters adjacency_radius:-1 detect_sign:-1 detect_threshold:' num2str(detect_threshold)]);
        %
        disp(' ')
        disp(' ')
        disp(['SORTED DATA FOR CHANNELS'])
        disp(' ')
        disp(' ')
end

%%
if qinspect

    system(['export PATH=/Applications/mountainlab-js/bin/:$PATH ; ' ...    % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/ephys-viz/bin/:$PATH ; ' ...   % Duplicate from ~/.bash_profile
        'export PATH=/usr/local/bin/:$PATH ; ' ...                          % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/:$PATH ; ' ...                 % Duplicate from ~/.bash_profile
        'export PATH=/usr/local/lib/python3.6/site-packages/:$PATH ; ' ...  % Duplicate from ~/.bash_profile
        'alias python=''python3'' ; ' ...                                   % Duplicate from ~/.bash_profile
        'source ~/ml_venv/bin/activate ; ' ...                              % Start previously created ml_venv virtual environment for running mountainlab
        'cd ~/Mountainsort_test ; ' ...
        './quick_inspect_sort.py ' '_ALL' ' ' mda_loc])
    
end

%% 

if qmountainview
    chan_nr     = num2str(view_channel);
    
    system(['export PATH=/Applications/mountainlab-js/bin/:$PATH ; ' ...    % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/ephys-viz/bin/:$PATH ; ' ...   % Duplicate from ~/.bash_profile
        'export PATH=/usr/local/bin/:$PATH ; ' ...                          % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/:$PATH ; ' ...                 % Duplicate from ~/.bash_profile
        'export PATH=/usr/local/lib/python3.6/site-packages/:$PATH ; ' ...  % Duplicate from ~/.bash_profile
        'export PATH=~/.mountainlab/packages/qt-mountainview/bin/:$PATH ; ' ...
        'alias python=''python3'' ; ' ...                                   % Duplicate from ~/.bash_profile
        'source ~/ml_venv/bin/activate ; ' ...                              % Start previously created ml_venv virtual environment for running mountainlab/view
        'qt-mountainview --raw=' mda_loc filesep file_prefix '_ALL' '.mda --firings=' mda_loc filesep 'firings' chan_nr '.mda.prv --samplerate=30000']);
    
end

