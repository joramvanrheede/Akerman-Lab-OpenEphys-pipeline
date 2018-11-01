% Plasticity 

LTP_expt_dir                = '/Volumes/Akermanlab/Joram/LTPtest'; % where are the data files? 

LTP_expt_day                = '2018_10_11'; % This script only deals with data from 1 day; it will read in pre and post files

%%  parameters 
whisk_resp_win    	= [0.005 0.200];    % Window for determining spiking responses (peak is determined within this time window; binned response is area under the psth curve in this window)
whisk_trace_win     = [0.005 0.300];    % Window for what time window of the PSTH to plot

whisk_P1_win        = [0.002 0.300];    % Window for determining the P1 LFP response (= maximum within this time window)
whisk_N1_win        = [0.002 0.150];    % Window for determining the N1 LFP response (= minimum within this time window)
n_per_point         = [10];              % How many repeats to average per data point; 1 = plot every data point individually
target_channels     = [1:32];           % Which channels to analyse

hist_binsize        = 0.001;            % Bin size for psths
rate_window         = [0.015];          % Window for gaussian smoothing for determining peak spike rate; currently not used apart from as input for ephys_data_psths function

peak_smooth_win   	= [7];              % Smoothing window for determining peak response
trace_smooth_win    = [20];              % Smoothing window for plotting PSTHs

%% To do - correct for spontaneous (subtract)?


%%
close all
%%
expt_folders                = dir(LTP_expt_dir);
qremove                     = ismember({expt_folders.name},{'.','..','.DS_Store'});
expt_folders(qremove)       = [];

q_test                      = ismember({expt_folders.name},{'Test_12_pre','Test_12_post_1','Test_12_post_2'});

test_folders                = {expt_folders(q_test).name};

experiment_files            = [];
expt_files_order            = [];
for a = 1:length(test_folders)
    this_folder                 = test_folders{a};
    this_path                   = [LTP_expt_dir filesep this_folder];
    day_folders                 = dir(this_path);
    q_day                       = ismember({day_folders.name},{LTP_expt_day});
    day_folder                  = day_folders(q_day);
    if ~isempty(day_folder)
        day_path                    = [this_path filesep day_folder.name];
        day_files               	= dir(day_path);
        qremove                     = ismember({day_files.name},{'.','..','.DS_Store'});
        day_files(qremove)          = [];
        for b = 1:length(day_files)
            experiment_files            = [experiment_files; {[this_path filesep day_folder.name filesep day_files(b).name]}];
            expt_files_order            = [expt_files_order; {day_files(b).name}];
        end
    end
end

[sort_temp sort_inds]   = sort(expt_files_order);
experiment_files        = experiment_files(sort_inds); % arrange by order of recording number

start_time                  = [];

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
for b = 1:length(experiment_files)
    this_file           = experiment_files{b};
    
    load(this_file)
    
    ephys_data.block_length     = ephys_data.trial_interval * length(ephys_data.conditions);
    
    % for the first experimental file
    if isempty(cond_stims)
        cond_stims  = [ephys_data.conditions.whisk_stimulator];
        cond_amps   = [ephys_data.conditions.whisk_amplitude];
    end
    
    file_name           = ephys_data.data_folder;
    
    date_ind            = strfind(file_name,'_20');
    date_start          = date_ind(1) + 1;
    date_str            = file_name(date_start:date_start+18);
    time_str            = date_str(12:19);
    
    this_year           = str2num(date_str(1:4));
    this_month          = str2num(date_str(6:7));
    this_day          	= str2num(date_str(9:10));
    
    this_hour           = str2num(time_str(1:2));
    this_min            = str2num(time_str(4:5));
    this_sec            = str2num(time_str(7:8));
    
    protocol_time       = datetime(this_year,this_month,this_day,this_hour,this_min,this_sec);
    
    
    %
    if isempty(start_time)
        start_time  = protocol_time;
    end
    
    time_counter    = minutes(protocol_time - start_time);
    
    expt_start_times    = [expt_start_times; time_counter];
    expt_end_times      = [expt_end_times; time_counter + ephys_data.protocol_duration/60];
    
    
    ephys_data          = ephys_data_psths(ephys_data, hist_binsize, rate_window); % this function adds psths to the ephys_data struct
    
    
    expt_whisk_binned       = [];
    expt_whisk_peak         = [];
    expt_whisk_P1           = [];
    expt_whisk_N1           = [];
    
    expt_whisk_psth         = [];
    expt_whisk_LFP          = [];
    expt_whisk_psth_vect    = [];
    expt_whisk_LFP_vect     = [];
    
    for c = 1:length(ephys_data.conditions)
        cond_data           = ephys_data.conditions(c);
        
        whisk_bins          = cond_data.psth_bins - cond_data.whisk_onset;
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
        
        for d = 1:n_points
            
            target_points        	= [1:n_per_point] + (d - 1) * n_per_point;
            
            %
            if max(target_points) > size(cond_data.trial_rate_psths,2)
                target_points   = [min(target_points):size(cond_data.trial_rate_psths,2)];
            end
            
            whisk_resp_selection    = smooth(squeeze(mean(mean(cond_data.trial_rate_psths(target_channels,target_points,q_whisk_win),1),2)),peak_smooth_win,'loess');
            whisk_resp_binned(d)    = mean(whisk_resp_selection(:));
            whisk_peak_rate(d)      = max(whisk_resp_selection(:));
            
            P1_LFP_mean             = squeeze(mean(mean(cond_data.LFP_trace(target_channels,target_points,q_P1_win),1),2));
            whisk_P1(d)             = max(P1_LFP_mean);
            
            N1_LFP_mean             = squeeze(mean(mean(cond_data.LFP_trace(target_channels,target_points,q_N1_win),1),2));
            whisk_N1(d)             = min(N1_LFP_mean);
            
            if c == 1
                time_counter        = time_counter + n_per_point * ephys_data.block_length / 60;
                time_vect           = [time_vect; time_counter];
            end
        end
        
        expt_whisk_binned       = [expt_whisk_binned whisk_resp_binned];
        expt_whisk_peak         = [expt_whisk_peak whisk_peak_rate];
        expt_whisk_P1           = [expt_whisk_P1 whisk_P1];
        expt_whisk_N1           = [expt_whisk_N1 whisk_N1];
        
        expt_whisk_psth         = [expt_whisk_psth smooth(squeeze(mean(mean(cond_data.trial_rate_psths(target_channels,:,q_whisk_trace),1),2)),trace_smooth_win)];
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
uniq_amps   = unique(cond_amps);
for a = 1:length(uniq_stims)
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
        max_N1              = max(max(multi_expt_whisk_N1(:,q_stim)));
        
        
        figure(1)
        subplot(1,2,a)
        title(['Binned whisk response stim ' num2str(a)])
        plot(time_vect,this_whisk_binned,'.','MarkerSize',20,'Color',1-[0.2 0.2 0.2]*b)
        hold on
        plot(this_psth_vect(:),1.1*max_binned + this_whisk_psth(:)/6,'Color',1-[0.2 0.2 0.2]*b)
        
        figure(2)
        subplot(1,2,a)
        title(['Peak whisk response stim ' num2str(a)])
        plot(time_vect,this_whisk_peak,'.','MarkerSize',20,'Color',1-[0.2 0.2 0.2]*b)
        hold on
        plot(this_psth_vect(:),1.1*max_peak + this_whisk_psth(:)/2,'Color',1-[0.2 0.2 0.2]*b)
        
        figure(3)
        subplot(1,2,a)
        title(['P1 whisk response stim ' num2str(a)])
        plot(time_vect,this_whisk_P1,'.','MarkerSize',20,'Color',1-[0.2 0.2 0.2]*b)
        hold on
        plot(this_LFP_vect(:),1.1*max_peak + this_whisk_LFP(:)/2,'Color',1-[0.2 0.2 0.2]*b)
        
        figure(4)
        subplot(1,2,a)
        title(['N1 whisk response stim ' num2str(a)])
        plot(time_vect,-this_whisk_N1,'.','MarkerSize',20,'Color',1-[0.2 0.2 0.2]*b)
        hold on
        plot(this_LFP_vect(:),1.1*max_peak + this_whisk_LFP(:)/2,'Color',1-[0.2 0.2 0.2]*b)
        
    end
end

% some figure 
for a = 1:4
    figure(a)
    set(gcf,'Units','Normalized')
    set(gcf,'Position',[0 .4 1 .3])
    % Set all y axes to the same range (based on the largest range)
    plotaxes    = get(gcf,'Children');
    % maxy        = cellfun(@max,get(plotaxes,'Ylim'));
    
    set(plotaxes,'XLim',[ -10, max(time_vect) + time_vect(1) + 10]) % add margin 
    
    hold off
end




