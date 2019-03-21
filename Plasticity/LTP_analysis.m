% Plasticity 

LTP_expt_dir                = '/Volumes/Akermanlab/Joram/Plasticity_8Hz_PW'; % where are the data files? 

LTP_expt_day                = '2019_03_08'; % This script only deals with data from 1 day; it will read in pre and post files

% control: 01/31, 02/01

%%  parameters 
whisk_resp_win              = [0.004 0.025];    % Window for determining spiking responses (peak is determined within this time window; binned response is area under the psth curve in this window)
whisk_trace_win             = [0.004 0.100];    % Window for what time window of the PSTH to plot

whisk_P1_win                = [0.004 0.200];    % Window for determining the P1 LFP response (= maximum within this time window)
whisk_N1_win                = [0.004 0.200];    % Window for determining the N1 LFP response (= minimum within this time window)

n_per_point                 = [1];              % How many repeats to average per data point; 1 = plot every data point individually
target_channels             = [1:32];           % Which channels to analyse (1 = most superficial, 32 = deepest)

hist_binsize                = 0.001;            % Bin size for psths
rate_window                 = [0.015];          % Window for gaussian smoothing for determining peak spike rate; currently not used apart from as input for ephys_data_psths function

peak_smooth_win             = [3];              % Smoothing window for determining peak response
trace_smooth_win            = [1];              % Smoothing window for plotting PSTHs



%% Make single LFP measure of peak-to-peak value
%% Make separate subplots for traces
%% Merge plots for different measures that use same trace:
    %% 1: spikes: Spike peak, Spikes binned
    %% 2: LFP: Peak-to-peak, N1, P1

%% Add standard error to plots when npoints > 1

%% Add mean and standard deviation for ALL points in protocol

%% Add figure splitting this out by channel group?

%% Printing
qprint              = 1;
figure_dir          = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/2019_Lab_meeting';
anaesth             = '100Hz_L23';

qrestart_count      = 0; 
%% To do - correct for spontaneous (subtract)?


%%
close all % Close any remaining open figure windows
%%
% get contents of the folder where the experimental data is kept
expt_folders                = dir(LTP_expt_dir); 

% remove some of the folder contents that we don't need
qremove                     = ismember({expt_folders.name},{'.','..','.DS_Store'}); 
expt_folders(qremove)       = []; 

% Look for folders with these names; these will contain the relevant data
expt_type_names             = {'Test_12_pre','Test_12_post_1','Test_12_post_2','Test_12_post_12','Test_1_post_1','Test_1_pre', 'Test_1_pre_1','Test_1_pre_0', 'Test_1_post_0'};
[q_test,which_type]       	= ismember({expt_folders.name},expt_type_names);
test_folders                = {expt_folders(q_test).name};
types                       = which_type(q_test);

% Set up some empty variables that we will use in the loop below
experiment_files            = [];
expt_files_order            = [];
expt_file_type              = [];

% Loop through the folders of different experiment types
for a = 1:length(test_folders)
    this_folder                 = test_folders{a};
    this_path                   = [LTP_expt_dir filesep this_folder];
    day_folders                 = dir(this_path);
    
    % Find the folder that corresponds to the selected experimental day
    q_day                       = ismember({day_folders.name},{LTP_expt_day});
    day_folder                  = day_folders(q_day);
    if ~isempty(day_folder)
        day_path                    = [this_path filesep day_folder.name];
        day_files               	= dir(day_path);
        
        % remove clutter from the list of folder contents
        qremove                     = ismember({day_files.name},{'.','..','.DS_Store'});
        day_files(qremove)          = [];
        
        % Loop over the files in the folder for this date
        for b = 1:length(day_files)
            % add the files to the variables that keep track of
            experiment_files            = [experiment_files; {[this_path filesep day_folder.name filesep day_files(b).name]}];
            expt_files_order            = [expt_files_order; {day_files(b).name}];
            expt_file_type              = [expt_file_type; expt_type_names(types(a))];
        end
    end
end

% Sort the list of files in order of recording number
[sort_temp sort_inds]   = sort(expt_files_order);
experiment_files        = experiment_files(sort_inds);
expt_file_type          = expt_file_type(sort_inds);

% Which whisker stimulator was used for induction / whick whisker was
% induced?
if any(ismember(expt_file_type,{'Test_12_post_1'}))
    RWS_whisker     = 1;
elseif any(ismember(expt_file_type,{'Test_12_post_2'}))
    RWS_whisker     = 2;
elseif any(ismember(expt_file_type,{'Test_1_post_1'}))
    RWS_whisker     = 1;
else
    RWS_whisker     = 0;
    warning('No whisker stimulated - correct?')
end

% Set up some empty variables that will be used in the loop below
start_time          	= [];

cond_stims            	= [];
cond_amps             	= [];
expt_start_times      	= [];
expt_end_times        	= [];

multi_expt_whisk_binned = [];
multi_expt_whisk_peak   = [];
multi_expt_whisk_P1     = [];
multi_expt_whisk_N1     = [];
multi_expt_whisk_psth   = [];
multi_expt_whisk_LFP    = [];
multi_expt_LFP_vect     = [];
multi_expt_psth_vect    = [];

time_vect               = [];

N1_pre_y_stim           = [];
N1_pre_n_stim           = [];

spk_pre_y_stim          = [];
spk_pre_n_stim          = [];

N1_post_y_stim          = [];
N1_post_n_stim          = [];

spk_post_y_stim         = [];
spk_post_n_stim         = [];

N1err_pre_y_stim        = [];
N1err_pre_n_stim        = [];

spkerr_pre_y_stim       = [];
spkerr_pre_n_stim       = [];

N1err_post_y_stim       = [];
N1err_post_n_stim       = [];

spkerr_post_y_stim      = [];
spkerr_post_n_stim      = [];

%% Loop over the relevant experiment files
for b = 1:length(experiment_files)
    
    % Select file and recording type from the list
    this_file           = experiment_files{b};
    this_type           = expt_file_type{b};
    
    % Load file
    load(this_file)
    
    % Determine length of block of trials (used to determine spacing of plot points)
    ephys_data.block_length     = ephys_data.trial_interval * length(ephys_data.conditions);
    
    % for the first experimental file
    if isempty(cond_stims)
        cond_stims  = [ephys_data.conditions.whisk_stimulator];
        cond_amps   = [ephys_data.conditions.whisk_amplitude];
    end
    
    % Get the name of the folder with the experimental data; this has
    % information about the date and time of the recordings which is used
    % to plot the timeline of the experiment
    folder_name       	= ephys_data.data_folder;
    
    % Extract date and time information from the folder name:
    % 1 - find the string that has the date and time
    date_ind            = strfind(folder_name,'_20');
    date_start          = date_ind(1) + 1;
    date_str            = folder_name(date_start:date_start+18);
    time_str            = date_str(12:19);
    
    % 2 - break down the date and time into individual values
    this_year           = str2num(date_str(1:4));
    this_month          = str2num(date_str(6:7));
    this_day          	= str2num(date_str(9:10));
    
    this_hour           = str2num(time_str(1:2));
    this_min            = str2num(time_str(4:5));
    this_sec            = str2num(time_str(7:8));
    
    % Make a time variable that can be compared to the start time of the
    % first protocol that was run
    protocol_time       = datetime(this_year,this_month,this_day,this_hour,this_min,this_sec);
    
    % If this is the first protocol, use this as our 'start time'
    % reference, t = 0
    if isempty(start_time)
        start_time  = protocol_time;
    end
    
    % Get elapsed time in minutes between start of this protocol and start
    % time
    time_counter        = minutes(protocol_time - start_time);
    
    expt_start_times    = [expt_start_times; time_counter];
    expt_end_times      = [expt_end_times; time_counter + ephys_data.protocol_duration/60];
    
    %% 
    
    % use the spike time data in ephys_data to generate post stimulus time
    % histograms for all trials and conditions so that we can easily
    % extract information about spike counts / rates
    ephys_data              = ephys_data_psths(ephys_data, hist_binsize, rate_window); % this function adds psths to the ephys_data struct
    
    % What amplitudes where the whisker stimuli?
    amplitudes              = ephys_data.condition_values(:,7);
    
    % What is the maximal whisker stimulus amplitude?
    max_amp                 = max(amplitudes);
    
    % Set up empty variables / clear variables that we will use in the loop
    % below
    expt_whisk_binned       = [];
    expt_whisk_peak         = [];
    expt_whisk_P1           = [];
    expt_whisk_N1           = [];
    expt_whisk_P2P          = [];
    
    expt_whisk_P1_serr      = [];
    expt_whisk_N1_serr      = [];
    expt_whisk_P2P_serr     = [];
    
    expt_whisk_psth         = [];
    expt_whisk_LFP          = [];
    expt_whisk_psth_vect    = [];
    expt_whisk_LFP_vect     = [];
    
    % Loop over the different experimental conditions (whisker stimulator 
    % nr, whisker stimulus amplitude, etc)
    for c = 1:length(ephys_data.conditions)
        
        % Make variable for data from this condition only
        cond_data           = ephys_data.conditions(c);
        
        % subtract whisker time so that psth bins are relative to whisker
        % onset (whisker onset == time 0)
        whisk_bins          = cond_data.psth_bins - cond_data.whisk_onset;
        
        % Make a boolean for selecting data from the relevant time window
        % relative to the whisker stimulus
        q_whisk_win         = whisk_bins > min(whisk_resp_win) & whisk_bins < max(whisk_resp_win);
        q_whisk_trace       = whisk_bins > min(whisk_trace_win) & whisk_bins < max(whisk_trace_win);
        
        
        q_P1_win            = whisk_bins > min(whisk_P1_win) & whisk_bins < max(whisk_P1_win);
        q_N1_win            = whisk_bins > min(whisk_N1_win) & whisk_bins < max(whisk_N1_win);
        
        whisk_rate_bins     = cond_data.spike_rate_times - cond_data.whisk_onset;
        q_whisk_rate_win    = whisk_rate_bins > min(whisk_resp_win) & whisk_rate_bins < max(whisk_resp_win);
        
        this_psth           = cond_data.trial_psths;
        
        n_points            = ceil(cond_data.n_trials / n_per_point);
        
        whisk_resp_binned   = zeros(n_points,1);
        whisk_peak_rate     = zeros(n_points,1);
        whisk_P1            = zeros(n_points,1);
        whisk_N1            = zeros(n_points,1);
        
        cond_stim           = repmat(cond_data.whisk_stimulator,n_points,1);
        cond_amp            = repmat(cond_data.whisk_amplitude,n_points,1);
        
        %% for comparison with other experiments
    	whisk_resp_selection        = cond_data.trial_psths(target_channels,:,q_whisk_win);
        whisk_resp_selection        = squeeze(mean(whisk_resp_selection,1));
        whisk_resp_selection_binned = mean(whisk_resp_selection,2);
        
        all_whisk_resp_binned       = mean(whisk_resp_selection_binned(:));
        all_whisk_resp_serr         = serr(whisk_resp_selection_binned(:));

        N1_LFP_select               = cond_data.LFP_trace(target_channels,:,q_N1_win);
        N1_LFP_select               = squeeze(mean(N1_LFP_select,1));
        
        N1_LFP_peaks                = min(N1_LFP_select,[],2);
        P1_LFP_peaks                = max(N1_LFP_select,[],2);
        
        N1_LFP_peaks                = P1_LFP_peaks - N1_LFP_peaks;
        
        all_whisk_N1                = mean(N1_LFP_peaks);
        all_whisk_N1_serr           = serr(N1_LFP_peaks);
        
        %% Find out whether this data point is pre / post, and which stimulator was used for induction
        
        if (strcmpi(this_type,'Test_12_pre') | strcmpi(this_type,'Test_1_pre') | strcmpi(this_type,'Test_1_pre_1') | strcmpi(this_type,'Test_1_pre_0')) & cond_amp == max_amp
            
            if cond_data.whisk_stimulator == RWS_whisker
                N1_pre_y_stim       = [N1_pre_y_stim all_whisk_N1];
                N1err_pre_y_stim    = [N1err_pre_y_stim all_whisk_N1_serr];
                
                spk_pre_y_stim      = [spk_pre_y_stim all_whisk_resp_binned];
                spkerr_pre_y_stim 	= [spkerr_pre_y_stim all_whisk_resp_serr];
            else
                N1_pre_n_stim       = [N1_pre_n_stim all_whisk_N1];
                N1err_pre_n_stim    = [N1err_pre_n_stim all_whisk_N1_serr];
                
                spk_pre_n_stim      = [spk_pre_n_stim all_whisk_resp_binned];
                spkerr_pre_n_stim 	= [spkerr_pre_n_stim all_whisk_resp_serr];
            end
            
        elseif (strcmpi(this_type,'Test_12_post_1') | strcmpi(this_type,'Test_12_post_2') | strcmpi(this_type,'Test_1_post_0') | strcmpi(this_type,'Test_1_post_1')) & cond_amp == max_amp
            
            if cond_data.whisk_stimulator == RWS_whisker
                N1_post_y_stim    	= [N1_post_y_stim all_whisk_N1];
                N1err_post_y_stim  	= [N1err_post_y_stim all_whisk_N1_serr];
                
                spk_post_y_stim   	= [spk_post_y_stim all_whisk_resp_binned];
                spkerr_post_y_stim 	= [spkerr_post_y_stim all_whisk_resp_serr];
            else
                N1_post_n_stim     	= [N1_post_n_stim all_whisk_N1];
                N1err_post_n_stim 	= [N1err_post_n_stim all_whisk_N1_serr];
                
                spk_post_n_stim  	= [spk_post_n_stim all_whisk_resp_binned];
                spkerr_post_n_stim	= [spkerr_post_n_stim all_whisk_resp_serr];
            end
            
        end
        
        
        %% Get data for individual points to plot. n_points is total number
        % of points to plot; 
        
        for d = 1:n_points
            
            target_points        	= [1:n_per_point] + (d - 1) * n_per_point;
            
            %
            if max(target_points) > size(cond_data.trial_psths,2)
                target_points   = [min(target_points):size(cond_data.trial_psths,2)];
            end
            
            whisk_resp_selection        = smooth(squeeze(mean(mean(cond_data.trial_psths(target_channels,target_points,q_whisk_win),1),2)),peak_smooth_win,'loess');
            whisk_resp_binned(d)        = mean(whisk_resp_selection(:));
            whisk_peak_rate(d)          = max(whisk_resp_selection(:));
            
            P1_LFP_mean                 = squeeze(mean(mean(cond_data.LFP_trace(target_channels,target_points,q_P1_win),1),2));
            whisk_P1(d)                 = max(P1_LFP_mean);
            
            N1_LFP_mean                 = squeeze(mean(mean(cond_data.LFP_trace(target_channels,target_points,q_N1_win),1),2));
            whisk_N1(d)                 = min(N1_LFP_mean);
            
            whisk_P2P(d)                = whisk_P1(d) - whisk_N1(d);
            
            % create some variables of 'NaN' values of the right size
            whisk_resp_binned_single    = NaN(length(target_points),1);
            whisk_peak_rate_single      = NaN(length(target_points),1);
            whisk_P1_single             = NaN(length(target_points),1);
            whisk_N1_single             = NaN(length(target_points),1);
            whisk_P2P_single            = NaN(length(target_points),1);
            for e = 1:length(target_points)
                % Which point are we dealing with?
                this_point                      = target_points(e);
                
                % Select the psth information for these channels, these
                % trials, and over the relevant time window
                whisk_resp_selection_single     = smooth(squeeze(mean(mean(cond_data.trial_psths(target_channels,this_point,q_whisk_win),1),2)),peak_smooth_win,'loess');
                
                % Take the area under the histogram in this window (total spike count) and save it in a variable
                whisk_resp_binned_single(e)     = mean(whisk_resp_selection_single(:));
                
                % Determine the peak of the histogram (approx peak spike rate)
                whisk_peak_rate_single(e)       = max(whisk_resp_selection_single(:));
                
                % Select LFP data for these channels, these trials, and
                % over the relevant time window, for first positive peak of
                % LFP (P1)
                P1_LFP_mean_single              = squeeze(mean(mean(cond_data.LFP_trace(target_channels,this_point,q_P1_win),1),2));
                whisk_P1_single(e)              = max(P1_LFP_mean_single);
                
                % Select LFP data for these channels, these trials, and
                % over the relevant time window, for first negative peak of
                % LFP (N1)
                N1_LFP_mean_single              = squeeze(mean(mean(cond_data.LFP_trace(target_channels,this_point,q_N1_win),1),2));
                whisk_N1_single(e)              = min(N1_LFP_mean_single);
                
                % Determine peak to peak size of LFP response (from min to 
                % max)
                whisk_P2P_single(e)             = whisk_P1_single(e) - whisk_N1_single(e);
                
            end
            
            % Determine standard errors for each point
            P1_serr(d)      = serr(whisk_P1_single);
            N1_serr(d)      = serr(whisk_N1_single);
            P2P_serr(d)    	= serr(whisk_P2P_single);
            
            % Make a time vector for each point; only needs to be done for
            % one condition (as the other conditions will use the same time 
            % points)
            if c == 1
                time_counter        = time_counter + n_per_point * ephys_data.block_length / 60;
                time_vect           = [time_vect; time_counter];
            end
        end
        
        expt_whisk_binned       = [expt_whisk_binned whisk_resp_binned];
        expt_whisk_peak         = [expt_whisk_peak whisk_peak_rate];
        
        expt_whisk_P1           = [expt_whisk_P1 whisk_P1];
        expt_whisk_N1           = [expt_whisk_N1 whisk_N1];
        expt_whisk_P2P          = [expt_whisk_P2P whisk_P2P];
        
        expt_whisk_P1_serr      = [expt_whisk_P1_serr P1_serr];
        expt_whisk_N1_serr      = [expt_whisk_N1_serr N1_serr];
        expt_whisk_P2P_serr     = [expt_whisk_P2P_serr P2P_serr];
        
        expt_whisk_psth         = [expt_whisk_psth smooth(squeeze(mean(mean(cond_data.trial_psths(target_channels,:,q_whisk_trace),1),2)),trace_smooth_win)];
        expt_whisk_LFP          = [expt_whisk_LFP squeeze(mean(mean(cond_data.LFP_trace(target_channels,:,q_P1_win),1),2));];
        
        expt_whisk_psth_vect    = [expt_whisk_psth_vect linspace(expt_start_times(b),expt_end_times(b),size(expt_whisk_psth,1))'];
        expt_whisk_LFP_vect     = [expt_whisk_LFP_vect linspace(expt_start_times(b),expt_end_times(b),size(expt_whisk_LFP,1))'];
        
    end
    
    multi_expt_whisk_binned = [multi_expt_whisk_binned; expt_whisk_binned];
    multi_expt_whisk_peak   = [multi_expt_whisk_peak; expt_whisk_peak];
    multi_expt_whisk_P1     = [multi_expt_whisk_P1; expt_whisk_P1];
    multi_expt_whisk_N1     = [multi_expt_whisk_N1; expt_whisk_N1];

    multi_expt_whisk_psth 	= cat(3,multi_expt_whisk_psth,expt_whisk_psth);
    multi_expt_psth_vect    = cat(3,multi_expt_psth_vect,expt_whisk_psth_vect);
    multi_expt_whisk_LFP	= cat(3,multi_expt_whisk_LFP,expt_whisk_LFP);
    multi_expt_LFP_vect     = cat(3,multi_expt_LFP_vect,expt_whisk_LFP_vect);
    
end


psth_mat_size                   = size(multi_expt_whisk_psth);
multi_expt_whisk_psth           = [multi_expt_whisk_psth; NaN(1,psth_mat_size(2),psth_mat_size(3))];
    
LFP_mat_size                    = size(multi_expt_whisk_LFP);
multi_expt_whisk_LFP            = [multi_expt_whisk_LFP; NaN(1,LFP_mat_size(2),LFP_mat_size(3))];

multi_expt_psth_vect            = [multi_expt_psth_vect; NaN(1,psth_mat_size(2),psth_mat_size(3))];

multi_expt_LFP_vect             = [multi_expt_LFP_vect; NaN(1,LFP_mat_size(2),LFP_mat_size(3))];

time_vect                       = time_vect - .5 * (n_per_point * ephys_data.block_length / 60);

uniq_stims  = unique(cond_stims);
n_stims     = length(uniq_stims);
uniq_amps   = unique(cond_amps);
for a = 1:n_stims
    q_stim      = cond_stims == uniq_stims(a);
    for b = 1:length(uniq_amps)
        q_amp   = cond_amps == uniq_amps(b);
        
        this_whisk_binned   = multi_expt_whisk_binned(:,q_amp & q_stim);
        this_whisk_peak     = multi_expt_whisk_peak(:,q_amp & q_stim);
        
        this_whisk_P1       = multi_expt_whisk_P1(:,q_amp & q_stim);
        this_whisk_N1       = multi_expt_whisk_N1(:,q_amp & q_stim);
        
        this_whisk_psth     = multi_expt_whisk_psth(:,q_amp & q_stim,:);
        this_psth_vect      = multi_expt_psth_vect(:,q_amp & q_stim,:);
        
        this_whisk_LFP      = multi_expt_whisk_LFP(:,q_amp & q_stim,:);
        this_LFP_vect       = multi_expt_LFP_vect(:,q_amp & q_stim,:);
        
        max_binned          = max(max(multi_expt_whisk_binned(:,q_stim)));
        max_peak            = max(max(multi_expt_whisk_peak(:,q_stim)));
        max_P1              = max(max(multi_expt_whisk_P1(:,q_stim)));
        max_N1              = max(max(-multi_expt_whisk_N1(:,q_stim)));
        
        
        figure(1)
        
        subplot(3,n_stims,a)
        plot(this_psth_vect(:),this_whisk_psth(:),'Color',[1 1 1] -uniq_amps(b)/max(uniq_amps))
        %barh = bar(this_psth_vect(:),this_whisk_psth(:),'histc','k');
        ylabel(['PSTH stimulator nr ' num2str(a)])
        
        subplot(3,n_stims,n_stims+a)
        plot(time_vect,this_whisk_binned,'.','MarkerSize',20,'Color',[1 1 1] - uniq_amps(b)/max(uniq_amps))
        ylabel(['Spike count stimulator ' num2str(a)])

        subplot(3,n_stims,2*n_stims+a)
        plot(time_vect,this_whisk_peak,'.','MarkerSize',20,'Color',[1 1 1] -uniq_amps(b)/max(uniq_amps))
        ylabel(['Peak spike rate stimulator' num2str(a)])
        
        figure(2)
        
        subplot(4,n_stims,a)
        plot(this_LFP_vect(:),this_whisk_LFP(:),'Color',[1 1 1]-uniq_amps(b)/max(uniq_amps))
        ylabel(['LFP simulator nr ' num2str(a)])
        
        subplot(4,n_stims,n_stims+a)
        plot(time_vect,this_whisk_P1 - this_whisk_N1,'.','MarkerSize',20,'Color',[1 1 1]-uniq_amps(b)/max(uniq_amps))
        ylabel(['LFP peak to peak simulator nr ' num2str(a)])
        
        subplot(4,n_stims,2*n_stims+a)
        plot(time_vect,this_whisk_P1,'.','MarkerSize',20,'Color',[1 1 1]-uniq_amps(b)/max(uniq_amps))
        ylabel(['LFP P1 simulator nr ' num2str(a)])
        
        subplot(4,n_stims,3*n_stims+a)
        plot(time_vect,this_whisk_N1,'.','MarkerSize',20,'Color',[1 1 1]-uniq_amps(b)/max(uniq_amps))
        ylabel(['LFP N1 simulator nr ' num2str(a)])
    end
end

% some figure things
for a = 1:2
    figure(a)
    set(gcf,'Units','Normalized')
    set(gcf,'Position',[0 .4 1 .3])
    % Set all y axes to the same range (based on the largest range)
    plotaxes    = get(gcf,'Children');
    % maxy        = cellfun(@max,get(plotaxes,'Ylim'));
    
%     set(plotaxes,'XLim',[ -10, max(time_vect) + time_vect(1) + 10]) % add margin 
%     
    hold off

end


if qprint
    
    figure(1)
    print(gcf,'-depsc','-r300',[figure_dir filesep 'Plasticity Spikes Ch' num2str(min(target_channels)) '-' num2str(max(target_channels)) ' ' LTP_expt_day ' ' anaesth])
    figure(2)
    print(gcf,'-depsc','-r300',[figure_dir filesep 'Plasticity LFP Ch' num2str(min(target_channels)) '-' num2str(max(target_channels)) ' ' LTP_expt_day ' ' anaesth])

end
