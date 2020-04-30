% 30th of April 2020 14:17

% The overarching folders containing the group folders which in turn contain
% the sorted data as generated by SYNC_CURATED_DATA and the preprocessed 
% multiunit and LFP data as generated by PREPROCESS_MULTIUNIT.
% e.g.
% EXPERIMENTAL_GROUP/EXPERIMENTAL_PREPS/EXPERIMENT_TITLES/EXPERIMENT_FILES.MAT
sorted_folder       = '/Users/Joram/Data/Sorting/Synced_individual';
multiunit_folder    = '/Users/Joram/Data/Preprocessed';

% The names of the group folders (also used in plotting as group names)
group_folders     	= {'POM' 'M1' 'S1'}; %

% Name of target experiment title folder
expt_folder_name  	= 'Timing';

% Reload and re-merge the single and multiunit data?
q_reload            = true; % Reload all data?

%% Inclusion settings and criteria
target_opto_rate  	= 60;   % Target opto response in Hz % 60? 100?
opto_p_threshold    = 0.1; % p value for opto response to include
min_opto_rate       = 10;  	% minimal opto response rate

% Max and min depth for units to include
max_depth           = -25;  % Depth relative to L4 sink (negative numbers = towards pia)
min_depth           = -400; % Depth relative to L4 sink (negative numbers = towards pia)

whisk_resp_p_thresh = 1;  % p-value threshold for increase in spike probability to whisk stimulus

L5_chans            = 5:12; % Layer 5 channels relative to the L4 sink (sink channel + [L5_chans])

% To do:
% Produce mean +/- serr line across units
% Consider normalised delta spike rate
% Consider normalising to max response / to max delta t
% Get response to control opto only - first opto resp
% Normalise data --> just use raw rates and subtract

%% Reload all the data if required
if q_reload
    % Clear old data
    clear merged_groups
    
    full_sorted_folders = [];
    full_multiunit_folders = [];
    for i = 1:length(group_folders)
        full_sorted_folders{i}      = [sorted_folder filesep group_folders{i} filesep expt_folder_name];
        full_multiunit_folders{i}   = [multiunit_folder filesep group_folders{i} filesep expt_folder_name];
    end
    
    % Load single unit data
    disp('Loading sorted spike data')
    sorted_groups = load_sorted_experiments(full_sorted_folders);
    
    % Load multiunit data. This function also adds current source density analysis
    % so it can then discard the LFP data
    disp('Loading multiunit spike data')
    multi_groups = load_multiunit_spikes(full_multiunit_folders);
    
    % A series of loops that uses sorted_groups as a reference and tries to find
    % corresponding multiunit_groups experiment data for each sorted data file
    % and then merges them.
    for a = 1:length(sorted_groups)
        this_group  = group_folders{a};
        for b = 1:length(sorted_groups(a).prep)
            for c = 1:length(sorted_groups(a).prep(b).expt_data)
                sorted_data     = sorted_groups(a).prep(b).expt_data(c).ephys_data;
                
                %%
                prep_nr     = [];
                expt_nr     = [];
                for i = 1:length(multi_groups(a).prep)
                    for j = 1:length(multi_groups(a).prep(i).expt_data)
                        
                        multi_unit_data         = multi_groups(a).prep(i).expt_data(j).ephys_data;
                        
                        multi_unit_folders_win  = split(multi_unit_data.data_folder,'\');
                        multi_unit_folders_mac  = split(multi_unit_data.data_folder,'/');
                        
                        single_unit_folders_win = split(sorted_data.data_folder,'\');
                        single_unit_folders_mac = split(sorted_data.data_folder,'/');
                        
                        if length(single_unit_folders_mac) > length(single_unit_folders_win)
                            single_unit_folder  = single_unit_folders_mac{end};
                        else
                            single_unit_folder  = single_unit_folders_win{end};
                        end
                        
                        if length(multi_unit_folders_mac) > length(multi_unit_folders_win)
                            multi_unit_folder  = multi_unit_folders_mac{end};
                        else
                            multi_unit_folder  = multi_unit_folders_win{end};
                        end
                        
                        if strcmp(single_unit_folder,multi_unit_folder)
                            prep_nr     = i;
                            expt_nr     = j;
                        end
                        
                    end
                end
                
                if isempty(prep_nr)
                    warning(['No multiunit data found for sorted data from ' sorted_data.data_folder])
                    continue
                end
                
                multi_data  = multi_groups(a).prep(prep_nr).expt_data(expt_nr).ephys_data;
                %%
                
                merged_groups(a).prep(b).expt_data(c) = merge_single_and_multi_unit(sorted_data, multi_data);
            end
        end
    end
    
    clear multi_groups
    % clear sorted_groups
end






for a = 1:length(merged_groups)
    group_data(a).unit_depths               = [];
    group_data(a).delta_t                   = [];
    
    group_data(a).spike_rate                = [];
    group_data(a).spike_rate_p              = [];
    group_data(a).spike_rate_ctrl           = [];
    
    group_data(a).spike_prob                = [];
    group_data(a).spike_prob_p              = [];
    group_data(a).spike_prob_ctrl           = [];
    
    group_data(a).first_spike_time          = [];
    group_data(a).first_spike_p             = [];
    group_data(a).first_spike_time_ctrl     = [];
    
    group_data(a).first_spike_jitter        = [];
    group_data(a).first_spike_jitter_ctrl   = [];
    
    group_data(a).peak_spike_rate           = [];
    group_data(a).peak_spike_rate_ctrl    	= [];
    
    group_units(a).spike_rate               = [];
    group_units(a).depths                   = [];
    
    group_data(a).select_power              = [];
    
    cum_n_units                             = 0;
    for b = 1:length(merged_groups(a).prep)
        
        % No multiunit data for this sorted data?
        if isempty(merged_groups(a).prep(b).expt_data)
            continue
        end
        
        opto_power  = [];
        opto_rate   = [];
        opto_stds   = [];
        opto_p      = [];
        for c = 1:length(merged_groups(a).prep(b).expt_data)
            
            ephys_data          = merged_groups(a).prep(b).expt_data(c);
            
            if isempty(ephys_data.conditions)
                opto_rate(c)    = 0;
                opto_stds(c)    = 0;
                opto_p(c)       = 1;
                opto_power(c)   = 0;
                continue
            end
            
            disp(ephys_data.data_folder)
            
            these_L5_chans      = ephys_data.LFP_min_chan + L5_chans;
            these_L5_chans(these_L5_chans>32) = [];
            
            timing_multi_data 	= timing_results(ephys_data,these_L5_chans);
            
               
            opto_rate(c)    = timing_multi_data.control_opto_rate;
            opto_stds(c)    = timing_multi_data.control_opto_stds;
            opto_p(c)       = timing_multi_data.opto_response_p;
            opto_power(c)   = ephys_data.conditions(1).LED_power;
            
        end
        
        if ~any(opto_p <= opto_p_threshold) | max(opto_rate) < min_opto_rate
            disp(['No opto resp in all of ' ephys_data.data_folder])
            continue
        end
        target_opto_diff    = abs(target_opto_rate - opto_rate);
        target_opto_ind    = find(target_opto_diff == min(target_opto_diff),1);
        
        disp(group_folders{a})
        select_opto         = opto_rate(target_opto_ind);
        
        ephys_data          = merged_groups(a).prep(b).expt_data(target_opto_ind);
        
        disp(['Using data from ' ephys_data.data_folder])
        
        L4_sink_microns  	= ephys_data.LFP_min_chan * 25

        delta_t             = [ephys_data.conditions.LED_onset] - [ephys_data.conditions.whisk_onset];
        delta_t             = round(delta_t,3) % round to nearest ms

        timing_unit_data                    = timing_by_unit(ephys_data);
        
        close all
        n_units                             = size(ephys_data.conditions(1).spikes,1);
        next_unit                           = n_units + 1;
        
        unit_inds                           = cum_n_units + [1:n_units];
        cum_n_units                         = cum_n_units + n_units;
        
        n_dts                               = size(timing_unit_data.spike_rate,2);
        max_dts                             = max(size(timing_unit_data.spike_rate,2),size(group_data(a).spike_rate,2));
        max_dts                             = 12;
        
        %% pre-allocate to max with NaNs
        group_data(a).spike_rate(unit_inds,1:max_dts)           = NaN;
        group_data(a).spike_rate_p(unit_inds,1:max_dts)         = NaN;
        
        group_data(a).spike_prob(unit_inds,1:max_dts)           = NaN;
        group_data(a).spike_prob_p(unit_inds,1:max_dts)         = NaN;
        
        group_data(a).first_spike_time(unit_inds,1:max_dts) 	= NaN;
        group_data(a).first_spike_p(unit_inds,1:max_dts)    	= NaN;
        
        group_data(a).first_spike_jitter(unit_inds,1:max_dts)	= NaN;
        
        group_data(a).peak_spike_rate(unit_inds,1:max_dts)  	= NaN;
        
        group_data(a).unit_depths(unit_inds)                    = NaN;
        group_data(a).delta_t(unit_inds,1:max_dts)          	= NaN;
        group_data(a).select_power(unit_inds)                	= NaN;
        
        %%
        group_data(a).spike_rate(unit_inds,1:n_dts)      	= timing_unit_data.spike_rate;
        group_data(a).spike_rate_p(unit_inds,1:n_dts)    	= timing_unit_data.spike_rate_p;
        group_data(a).spike_rate_ctrl(unit_inds)            = timing_unit_data.spike_rate(:,end);
        
        group_data(a).spike_prob(unit_inds,1:n_dts)      	= timing_unit_data.spike_probabilities;
        group_data(a).spike_prob_p(unit_inds,1:n_dts)       = timing_unit_data.spike_prob_p;
        group_data(a).spike_prob_ctrl(unit_inds)            = timing_unit_data.spike_probabilities(:,end);
        
        group_data(a).first_spike_time(unit_inds,1:n_dts) 	= timing_unit_data.first_spike_times;
        group_data(a).first_spike_p(unit_inds,1:n_dts)    	= timing_unit_data.first_spike_p;
        group_data(a).first_spike_time_ctrl(unit_inds)    	= timing_unit_data.first_spike_times(:,end);
        
        group_data(a).first_spike_jitter(unit_inds,1:n_dts)	= timing_unit_data.first_spike_jitter;
        group_data(a).first_spike_jitter_ctrl(unit_inds)    = timing_unit_data.first_spike_jitter(:,end);
        
        group_data(a).peak_spike_rate(unit_inds,1:n_dts)  	= timing_unit_data.peak_spike_rates;
        group_data(a).peak_spike_rate_ctrl(unit_inds)       = timing_unit_data.peak_spike_rates(:,end);
        
        group_data(a).unit_depths(unit_inds)              	= (800 - ephys_data.unit_depths) - L4_sink_microns;
        group_data(a).delta_t(unit_inds,1:n_dts)          	= repmat(delta_t,n_units,1);
        group_data(a).select_opto(unit_inds)             	= select_opto;
        
        group_data(a).spike_resp_prob_p(unit_inds)        	= timing_unit_data.spike_resp_prob_p;
        group_data(a).whisk_resp_stds(unit_inds)            = timing_unit_data.whisk_resp_stds;
        
    end
    
    
    [group_data(a).unit_depths, depth_sort_inds]    = sort(group_data(a).unit_depths);
    
    group_data(a).spike_rate                        = group_data(a).spike_rate(depth_sort_inds,:);
    group_data(a).spike_rate_p                      = group_data(a).spike_rate_p(depth_sort_inds,:);
    group_data(a).spike_rate_ctrl                 	= group_data(a).spike_rate_ctrl(depth_sort_inds);
    
    group_data(a).spike_prob                        = group_data(a).spike_prob(depth_sort_inds,:);
    group_data(a).spike_prob_p                     	= group_data(a).spike_prob_p(depth_sort_inds,:);
    group_data(a).spike_prob_ctrl                 	= group_data(a).spike_prob_ctrl(depth_sort_inds);
    
    group_data(a).first_spike_time              	= group_data(a).first_spike_time(depth_sort_inds,:);
    group_data(a).first_spike_p                     = group_data(a).first_spike_p(depth_sort_inds,:);
    group_data(a).first_spike_time_ctrl            	= group_data(a).first_spike_time_ctrl(depth_sort_inds);
    
    group_data(a).first_spike_jitter                = group_data(a).first_spike_jitter(depth_sort_inds,:);
    group_data(a).first_spike_jitter_ctrl           = group_data(a).first_spike_jitter_ctrl(depth_sort_inds);
    
    group_data(a).peak_spike_rate                   = group_data(a).peak_spike_rate(depth_sort_inds,:);
    group_data(a).peak_spike_rate_ctrl              = group_data(a).peak_spike_rate_ctrl(depth_sort_inds);
    
    group_data(a).delta_t                           = group_data(a).delta_t(depth_sort_inds,:);
    
    group_data(a).spike_resp_prob_p              	= group_data(a).spike_resp_prob_p(depth_sort_inds);
  	group_data(a).whisk_resp_stds               	= group_data(a).whisk_resp_stds(depth_sort_inds);
        
end


figure(1)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(2)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(3)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(4)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(5)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
for a = 1:length(merged_groups)
    unit_depths         = group_data(a).unit_depths;
    q_depth             = unit_depths <= max_depth & unit_depths >= min_depth;
    
    delta_t             = group_data(a).delta_t;
    uniq_delta_ts    	= unique(delta_t(:));
    uniq_delta_ts       = uniq_delta_ts(~isnan(uniq_delta_ts));
    
    %% Spike probability
    
    spike_prob          = group_data(a).spike_prob;
    spike_prob_ctrl     = group_data(a).spike_prob_ctrl;
    delta_spike_prob    = spike_prob - spike_prob_ctrl(:);
    
    spike_prob_p        = group_data(a).spike_prob_p;
    
    q_prob_p                        = spike_prob_p <= 0.05;
    sig_delta_spike_prob            = delta_spike_prob;
    sig_delta_spike_prob(~q_prob_p) = NaN;
    
    spike_prob_up       = delta_spike_prob >= 0;
    
    q_whisk_resp        = group_data(a).spike_resp_prob_p < whisk_resp_p_thresh;
    
    q_all   = q_depth & q_whisk_resp;
    
    if sum(q_all) == 0
        continue
    end
    
    %%
    spike_prob_means    = NaN(size(uniq_delta_ts));
    for i = 1:length(uniq_delta_ts)
        this_delta_t            = uniq_delta_ts(i);
        q_delta_t               = delta_t == this_delta_t;
        spike_prob_means(i)  	= nanmean(delta_spike_prob(q_delta_t & q_all'));
    end
    
    %%
    figure(1)
    subplot(1,length(merged_groups),a)
    plot(delta_t(q_all,:)',delta_spike_prob(q_all,:)','Color',[0 0 0 .2])
    hold on
    
    plot(uniq_delta_ts,spike_prob_means,'Color',[1 0 0],'LineWidth',4)
    
    title([group_folders{a} ', n = ' num2str(sum(q_all)) ' units'])
    plot(delta_t(q_all,:),sig_delta_spike_prob(q_all,:),'k.','MarkerSize',25)
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('Delta spike probability')
    xlabel('Delay')
    fixplot
    
    
    %% Spike rate
    spike_rate          = group_data(a).spike_rate;
    spike_rate_ctrl     = group_data(a).spike_rate_ctrl;
    delta_spike_rate    = spike_rate - spike_rate_ctrl(:);
    
    spike_rate_p        = group_data(a).spike_rate_p;
    
    q_rate_p                        = spike_rate_p <= 0.05;
    sig_delta_spike_rate            = delta_spike_rate;
    sig_delta_spike_rate(~q_rate_p) = NaN;
    
    spike_rate_up       = delta_spike_rate >= 0;
    
    %%
    spike_rate_means   = NaN(size(uniq_delta_ts));
    for i = 1:length(uniq_delta_ts)
        this_delta_t            = uniq_delta_ts(i);
        q_delta_t               = delta_t == this_delta_t;
        spike_rate_means(i)  	= nanmean(delta_spike_rate(q_delta_t & q_all'));
    end
    
    
    %%
    figure(2)
    subplot(1,length(merged_groups),a)
    plot(delta_t(q_all,:)',delta_spike_rate(q_all,:)','Color',[0 0 0 .2])
    hold on
    
    plot(uniq_delta_ts,spike_rate_means,'Color',[1 0 0],'LineWidth',4)
    
    title([group_folders{a} ', n = ' num2str(sum(q_all)) ' units'])
    plot(delta_t(q_all,:),sig_delta_spike_rate(q_all,:),'k.','MarkerSize',25)
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('Delta spike rate')
    xlabel('Delay')
    fixplot
    
    
    %% First spike time
    first_spike_time                = group_data(a).first_spike_time;
    first_spike_time_ctrl           = group_data(a).first_spike_time_ctrl;
    delta_first_spike_time          = first_spike_time - first_spike_time_ctrl(:);
    
    first_spike_p                   = group_data(a).first_spike_p;
    
    q_rate_p                        = first_spike_p <= 0.05;
    sig_delta_first_spike_time            = delta_first_spike_time;
    sig_delta_first_spike_time(~q_rate_p) = NaN;
    
    first_spike_time_up             = delta_first_spike_time >= 0;
    
        
    first_spike_means   = NaN(size(uniq_delta_ts));
    for i = 1:length(uniq_delta_ts)
        this_delta_t            = uniq_delta_ts(i);
        q_delta_t               = delta_t == this_delta_t;
        first_spike_means(i)  	= nanmean(delta_first_spike_time(q_delta_t & q_all'));
    end
    
    
    %%
    figure(3)
    subplot(1,length(merged_groups),a)
    plot(delta_t(q_all,:)',delta_first_spike_time(q_all,:)','Color',[0 0 0 .2])
    hold on
    
    plot(uniq_delta_ts,first_spike_means,'Color',[1 0 0],'LineWidth',4)
    
    title([group_folders{a} ', n = ' num2str(sum(q_all)) ' units'])
    plot(delta_t(q_all,:),sig_delta_first_spike_time(q_all,:),'k.','MarkerSize',25)
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('Delta first spike time')
    xlabel('Delay')
    fixplot
    
    %%
    peak_spike_rate         = group_data(a).peak_spike_rate;
    peak_spike_rate_ctrl    = group_data(a).peak_spike_rate_ctrl(:);
    delta_peak_rate         = peak_spike_rate - peak_spike_rate_ctrl;
    
    peak_rate_means         = NaN(size(uniq_delta_ts));
    for i = 1:length(uniq_delta_ts)
        this_delta_t            = uniq_delta_ts(i);
        q_delta_t               = delta_t == this_delta_t;
        peak_rate_means(i)      = nanmean(delta_peak_rate(q_delta_t & q_all'));
    end
    
    figure(4)
    subplot(1,length(merged_groups),a)
    plot(delta_t(q_all,:)',delta_peak_rate(q_all,:)','Color',[0 0 0 .5])
    hold on
    
    plot(uniq_delta_ts,peak_rate_means,'Color',[1 0 0],'LineWidth',4)
    
    title([group_folders{a} ', n = ' num2str(sum(q_all)) ' units'])
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('Peak spike rate')
    xlabel('Delay')
    fixplot
    
    %%
    first_spike_jitter         = group_data(a).first_spike_jitter;
    first_spike_jitter_ctrl    = group_data(a).first_spike_jitter_ctrl(:);
    delta_spike_jitter        	= first_spike_jitter - first_spike_jitter_ctrl;

    %%
    
    spike_jitter_means      = NaN(size(uniq_delta_ts));
    for i = 1:length(uniq_delta_ts)
        this_delta_t            = uniq_delta_ts(i);
        q_delta_t               = delta_t == this_delta_t;
        spike_jitter_means(i) 	= nanmean(delta_spike_jitter(q_delta_t & q_all'));
    end
    %%
    
    figure(5)
    subplot(1,length(merged_groups),a)
    plot(delta_t(q_all,:)',delta_spike_jitter(q_all,:)','Color',[0 0 0 .5])
    hold on
    
    plot(uniq_delta_ts,spike_jitter_means,'Color',[1 0 0],'LineWidth',4)
    
    title([group_folders{a} ', n = ' num2str(sum(q_all)) ' units'])
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('First spike jitter')
    xlabel('Delay')
    fixplot
    %
end

plot_depths     = [];
depth_groups  	= [];
plot_opto_resp 	= [];
opto_groups     = [];
for a = 1:length(group_data)
    plot_depths     = [plot_depths; group_data(a).unit_depths(:)];
    depth_groups    = [depth_groups; repmat(a,size(group_data(a).unit_depths(:)))];
    plot_opto_resp  = [plot_opto_resp; unique(group_data(a).select_opto(:))];
   	opto_groups     = [opto_groups; repmat(a,size(unique(group_data(a).select_opto(:))))];
    
end

figure
beeswarmplot(plot_depths, depth_groups, group_folders)
title('Unit depths by group')
ylabel('Unit depth')
fixplot

figure
beeswarmplot(plot_opto_resp, opto_groups, group_folders)
title('L5 opto response by group')
ylabel('Spike rate')
fixplot