

% channels            = [1:4 13:20 29:32];
detect_threshold    = 4;

range               = 1:(30000*600);

oe_data_loc         = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data/2018_09_18/ChR2-YFP-OFF_2018-09-18_15-54-30_6';
mda_loc             = [oe_data_loc filesep 'mda'];
geom_file           = [cd filesep '32ch_site_geometry.csv'];

file_prefix         = '100_CH';

qconvert            = 0; % (re)do conversion to mda?
qsort               = 0; % (re)do sorting?

view_channel        = 20; % 2018/04/05: 16,29,13

qinspect            = 0;
qmountainview       = 1;

probe_type          = '32ch'; % '16ch' or '32ch'

layer               = 'L5';

L2_3_chans          = [1:8];
L4_chans            = [10:14];
L5_chans            = [16:28];

% set this up for different probe types:
channel_order_16ch  = [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]; % 16ch linear silicon A16 probe neuronexus: [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]
channel_order_32ch  = [16 32 1 17 14 30 3 19 8 24 7 29 9 25 15 20 10 23 2 28 6 26 5 21 11 31 4 27 12 22 13 18]; % 32Ch linear A32 probe neuronexus [18 2 31 15 19 3 30 14 17 1 32 16 24 8 25 9 23 7 26 10 22 6 27 11 21 5 28 12 20 4 29 13];


if ~isdir(mda_loc)
    mkdir(mda_loc)
end
    
%% File conversion from openephys .continuous to .mda
if qconvert
    
    switch layer
        case 'L2_3'
            channels        = L2_3_chans;
        case 'L4'
            channels        = L4_chans;
        case 'L5'
            channels        = L5_chans;
    end
    
    geom_filenm     = [layer '_geom.csv'];
    geom_mat        = [[1:length(channels)]; zeros(1,length(channels))]';
    % determine which channel map to use
    switch probe_type
        case '16ch'
            channel_map     = channel_order_16ch;
            geom_mat        = geom_mat * 50;
        case '32ch'
            channel_map     = channel_order_32ch;
            geom_mat        = geom_mat * 25;
    end
    
    csvwrite([mda_loc filesep geom_filenm],geom_mat);
    
    oe2mda_multi(oe_data_loc,mda_loc,channel_map(channels),range,layer) % convert data for this mapped channel range to mda
    
end

%% Spike sorting
if qsort

        % make file name string
        switch layer
            case 'L2_3'
                write_file  	= [mda_loc filesep 'L2_3.mda'];
            case 'L4'
                write_file  	= [mda_loc filesep 'L4.mda'];
            case 'L5'
                write_file  	= [mda_loc filesep 'L5.mda'];
        end
        
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
            '__conda_setup="$(CONDA_REPORT_ERRORS=false ''/Users/Joram/anaconda3/bin/conda'' shell.bash hook 2> /dev/null)" ' ...
            'if [ $? -eq 0 ]; then ' ... 
            '    \eval "$__conda_setup"' ...
            'else' ... 
            '    if [ -f "/Users/Joram/anaconda3/etc/profile.d/conda.sh" ]; then ' ...
            '        . "/Users/Joram/anaconda3/etc/profile.d/conda.sh" ' ...
            '        CONDA_CHANGEPS1=false conda activate base ' ...
            '    else ' ...
            '        \export PATH="/Users/Joram/anaconda3/bin:$PATH" ' ...
            '    fi ' ...
            'fi ;' ...
            'unset __conda_setup ;' ...
            'conda activate mountainlab ; ' ...                                 % Start conda mountainlab virtual environment for running mountainlab
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
        './quick_inspect_sort.py ' '_' layer ' ' mda_loc])
    
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

