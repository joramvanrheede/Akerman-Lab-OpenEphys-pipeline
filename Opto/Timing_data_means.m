
%% This script will load the preprocessed data

% To do - Include option to detect threshold value in opto_pulse_function
% To do - Automatically select 'nearest' power on a by experiment basis (if we don't have a power for all experiments)
% Overlay of psths

%% Changing any of these parameters means you need to 'do reload'
do_reload           = true;

% List data folders for the experiment types
data_folders        = { '/Users/Joram/Data/Preprocessed/POM/Timing'};

whisk_resp_win      = [0.006 0.025];
psth_bins           = [-0.25:0.001:.5];
artifact_win        = [-0.002 0.006];

% threshold_units     = 'SD'; % SD or Hz
% threshold_measure   = 'rate'; % 'peak' or 'rate'
%% Response thresholding

%% These parameters can be changed without reloading

% Nicknames for the groups
group_names         = {'GTRB'};

% Response measure to plot / analyse
resp_measure_name 	= 'delta_spike_rate'; % 'delta_spike_rate', 'raw_spike_rate', 'peak_spike_rate', 'corr_peak_spike_rate', or 'peak_time'
target_power     	= 100; % [a number, or 'threshold' to use a threshold]
resp_threshold      = 2;

target_delta_t      = -0.02;


% Channels
L23_chans           = 1:7;
L5_chans            = 17:28;

% PSTH visualisation
psth_smoothwin      = [7]; % nr of samples
psth_smooth_method  = 'gaussian';

% save_figs           = false;
% fig_save_dir       	= '/Users/Joram/Dropbox/Akerman Postdoc/Figures/2019_09_25_LabM';
% expt_name           = 'Laser_pulse_';

%% Code execution starts here
close all;
%% Function: load pulse_data (with multiple time win inputs)



if do_reload
    groups  = struct;
    for a = 1:length(data_folders)
        groups(a).timing_data         = get_opto_timing_data(data_folders{a}, whisk_resp_win, psth_bins);
    end
end


group_number                    = [];
L23_all                         = [];
L5_all                          = [];
all_L23_resps                   = [];
all_L5_resps                    = [];
all_delta_ts                    = [];
plot_L23_resp_measure           = [];
plot_L5_resp_measure            = [];
plot_L23_spike_density_rates    = [];
plot_L23_control_density_rates  = [];
plot_L23_diff_density_rates     = [];
plot_L5_spike_density_rates     = [];
plot_L5_control_density_rates   = [];
plot_L5_diff_density_rates      = [];
all_psths                       = [];
for a = 1:length(groups)
    timing_data     = groups(a).timing_data;
    
    L23_psths   = [];
    L5_psths    = [];
    for b=1:length(timing_data)
        
        %% Select the experiment file with the right opto power
        opto_powers                     = [timing_data(b).experiment(:).opto_power];
        power_diffs                     = target_power - opto_powers;
        abs_power_diffs                 = abs(power_diffs);
        min_power_diff                  = min(abs_power_diffs);
        q_min_power_diff                = abs_power_diffs == min_power_diff;
        closest_power_ind               = find(q_min_power_diff,1,'first');
        best_power                      = opto_powers(closest_power_ind);
        
        
        
        %% Select the appropriate time difference between opto stim and whisk stim
        experiment                      = timing_data(b).experiment(closest_power_ind);
        
        
        
        %%
        delta_ts                        = experiment.delta_t;
        delta_t_diffs                   = target_delta_t - delta_ts;
        abs_delta_t_diffs               = abs(delta_t_diffs);
        min_delta_t_diff                = min(abs_delta_t_diffs);
        q_min_delta_t_diff              = abs_delta_t_diffs == min_delta_t_diff;
        closest_delta_t_ind             = find(q_min_delta_t_diff,1,'first');
        best_delta_t                    = delta_ts(closest_delta_t_ind);
        
        
        
        
        %% Response measure for comparison plots
        
        control_spike_rate          	= experiment.control_spike_rate;
        spont_spike_std              	= experiment.spont_spike_std;
        
        experiment.data_folder
        
        L23_control_spike_rate       	= mean(control_spike_rate(L23_chans))
        L23_spont_spike_std             = mean(spont_spike_std(L23_chans))            
        L23_resp                        = L23_control_spike_rate > resp_threshold * L23_spont_spike_std
        
        L5_control_spike_rate       	= mean(control_spike_rate(L5_chans))
        L5_spont_spike_std              = mean(spont_spike_std(L5_chans))             
        L5_resp                         = L5_control_spike_rate > resp_threshold * L5_spont_spike_std
        
        resp_measure                    = experiment.(resp_measure_name);
        
        
        if ~L23_resp
            disp([experiment.data_folder ' OUT!'])
            continue
        end
        
        all_delta_ts                    = [all_delta_ts; delta_ts(:)];
        
        L23_resps                       = mean(resp_measure(L23_chans,:))
        all_L23_resps                   = [all_L23_resps; L23_resps(:)];
        
        L5_resps                        = mean(resp_measure(L5_chans,:))
        all_L5_resps                    = [all_L5_resps; L5_resps(:)];
        
        L23_resp_measure                = mean(resp_measure(L23_chans,closest_delta_t_ind));
        
        L5_resp_measure                 = mean(resp_measure(L5_chans,closest_delta_t_ind));
        
        plot_L23_resp_measure           = [plot_L23_resp_measure; L23_resp_measure];
        plot_L5_resp_measure            = [plot_L5_resp_measure; L5_resp_measure];
        
        group_number                    = [group_number; a];
        
        L23_spike_density_rates         = mean(experiment.spike_density_rates(L23_chans,:,closest_delta_t_ind));
        L23_control_density_rates       = mean(experiment.spike_density_rates(L23_chans,:,end));
        L23_diff_density_rates          = L23_spike_density_rates - L23_control_density_rates;
        
        L5_spike_density_rates          = mean(experiment.spike_density_rates(L5_chans,:,closest_delta_t_ind));
        L5_control_density_rates      	= mean(experiment.spike_density_rates(L5_chans,:,end));
        L5_diff_density_rates           = L5_spike_density_rates - L5_control_density_rates;
        
        light_onset                     = 0 + best_delta_t;
        q_artifact                      = psth_bins >= (artifact_win(1) + light_onset) & psth_bins <= (artifact_win(2) + light_onset);
        
        L23_spike_density_rates(q_artifact)     = NaN;
        L23_control_density_rates(q_artifact)   = NaN;
        L23_diff_density_rates(q_artifact)      = NaN;
        
        L5_spike_density_rates(q_artifact)      = NaN;
        L5_control_density_rates(q_artifact)    = NaN;
        L5_diff_density_rates(q_artifact)       = NaN;
        
        
        plot_L23_spike_density_rates    = [plot_L23_spike_density_rates L23_spike_density_rates(:)];
        plot_L23_control_density_rates  = [plot_L23_control_density_rates L23_control_density_rates(:)];
        plot_L23_diff_density_rates     = [plot_L23_diff_density_rates L23_diff_density_rates(:)] ;
        plot_L5_spike_density_rates     = [plot_L5_spike_density_rates L5_spike_density_rates(:)];
        plot_L5_control_density_rates   = [plot_L5_control_density_rates L5_control_density_rates(:)];
        plot_L5_diff_density_rates      = [plot_L5_diff_density_rates L5_diff_density_rates(:)];
        
    end
%     L23_group_psths(a)  = {L23_psths};
%     L5_group_psths(a)   = {L5_psths};


end

% Make resp measure name a string suitable for plots
resp_measure_name(ismember(resp_measure_name,'_')) = ' ';
resp_measure_name(1) = upper(resp_measure_name(1));


%% 
[uniq_delta_ts, ~, group_nrs]    = unique(all_delta_ts);

figure
point_plot(all_L23_resps, group_nrs, uniq_delta_ts)
title(['Layer 23 - opto delay vs ' resp_measure_name])
fixplot
xlabel('Opto power')
ylabel(resp_measure_name)
xlim([-0.12 0.02])
hold on
plot(xlim,[0 0 ],'r:','LineWidth',2)

L23_p = [];
for a = 1:length(uniq_delta_ts)
    this_delta_t        = uniq_delta_ts(a);
    q_delta_t           = all_delta_ts == this_delta_t;
    
    these_L23_resps     = all_L23_resps(q_delta_t);
    
    [h(a),L23_p(a)] = ttest(these_L23_resps); % ttest against mean of 0
    
end
    

%%
figure
point_plot(all_L5_resps, group_nrs, uniq_delta_ts)
title(['Layer 5 - opto delay vs ' resp_measure_name])
fixplot
xlabel('Opto power')
ylabel(resp_measure_name)
xlim([-0.12 0.02])
hold on
plot(xlim,[0 0 ],'r:','LineWidth',2)

L5_p = [];
for a = 1:length(uniq_delta_ts)
    this_delta_t        = uniq_delta_ts(a);
    q_delta_t           = all_delta_ts == this_delta_t;
    
    these_L5_resps     = all_L5_resps(q_delta_t);
    
    [h(a),L5_p(a)] = ttest(these_L5_resps); % ttest against mean of 0
    
end

L23_p
L5_p

%% Compare data points and means for L23
figure
beeswarmplot(plot_L23_resp_measure,group_number,group_names)
hold on
plot(xlim,[0 0],'k:','LineWidth',2)
title(['L23; Delta t = ' num2str(target_delta_t) ';  Laser power = ' num2str(target_power)])
ylabel(resp_measure_name)
fixplot

%% Compare data points and means for L5
figure
beeswarmplot(plot_L5_resp_measure,group_number,group_names)
hold on
plot(xlim,[0 0],'k:','LineWidth',2)
title(['L5; Delta t = ' num2str(target_delta_t) ';  Laser power = ' num2str(target_power)])
ylabel(resp_measure_name)
fixplot

%% Post stimulus time histograms

figure
set(gcf,'Units','Normalized','Position',[.1 .4 .8 .4])
for a = 1:length(groups)
    subplot(length(groups),2,a)
    q_group         = group_number == a;
    this_psth       = plot_L23_spike_density_rates(:,q_group); 
    
    smooth_psths     = smoothdata(this_psth,1,'gaussian',psth_smoothwin);
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.25,'BarWidth',1);
        hold on
    end
    xlim([-0.1 0.5])
    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',3)
    title([group_names{a} ' - L23 PSTH'])
    fixplot
    
    xlim([0.000 0.1])
    % shaded_region(whisk_resp_win)
    
    title([group_names{a} 'L23 PSTH'])
    fixplot
    
    subplot(length(groups),2,length(groups)+a)
    control_psths           = plot_L23_control_density_rates(:,q_group); 
    
    smooth_psths     = smoothdata(control_psths,1,'gaussian',psth_smoothwin);
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.25,'BarWidth',1);
        hold on
    end
    xlim([-0.1 0.5])
    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',3)
    title([group_names{a} ' - L23 PSTH'])
    fixplot
    
    xlim([0.000 0.1])
    % shaded_region(whisk_resp_win)
    
    title([group_names{a} 'L23 Control PSTH'])
    fixplot
end
subplot_equal_y


figure
set(gcf,'Units','Normalized','Position',[.1 .4 .8 .4])
for a = 1:length(groups)
    subplot(length(groups),2,a)
    q_group         = group_number == a;
    this_psth       = plot_L5_spike_density_rates(:,q_group); 
    
    smooth_psths     = smoothdata(this_psth,1,'gaussian',psth_smoothwin);
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.25,'BarWidth',1);
        hold on
    end
    xlim([-0.1 0.5])
    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',3)
    title([group_names{a} ' - L5 PSTH'])
    fixplot
    
    xlim([0.000 0.1])
    % shaded_region(whisk_resp_win)
    
    title([group_names{a} 'L5 PSTH'])
    fixplot
    
    subplot(length(groups),2,length(groups)+a)
    control_psths           = plot_L5_control_density_rates(:,q_group); 
    
    smooth_psths     = smoothdata(control_psths,1,'gaussian',psth_smoothwin);
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.25,'BarWidth',1);
        hold on
    end
    xlim([-0.1 0.5])
    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',3)
    title([group_names{a} ' - L5 PSTH'])
    fixplot
    
    xlim([0.000 0.1])
    % shaded_region(whisk_resp_win)
    
    title([group_names{a} 'L5 Control PSTH'])
    fixplot
end
subplot_equal_y

%% L23 Diff plot

figure
%set(gcf,'Units','Normalized','Position',[.1 .4 .8 .4])
for a = 1:length(groups)
    subplot(1,length(groups),a)
    q_group         = group_number == a;
    this_psth       = plot_L23_diff_density_rates(:,q_group); 
    
    smooth_psths    = smoothdata(this_psth,1,'gaussian',psth_smoothwin);
    
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.4,'BarWidth',1);
        hold on
        %plot(psth_bins(1:end-1),smooth_psths,'LineWidth',3)
    end
    xlim([0.000 0.1])

    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',2)
    % shaded_region(whisk_resp_win)
    
    title([group_names{a} 'L23 Diff PSTH'])
    fixplot
    
end
subplot_equal_y

%% L5 Diff plot

figure
%set(gcf,'Units','Normalized','Position',[.1 .4 .8 .4])
for a = 1:length(groups)
    subplot(1,length(groups),a)
    q_group         = group_number == a;
    this_psth       = plot_L5_diff_density_rates(:,q_group); 
    
    smooth_psths    = smoothdata(this_psth,1,'gaussian',psth_smoothwin);
    
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.4,'BarWidth',1);
        hold on
        %plot(psth_bins(1:end-1),smooth_psths,'LineWidth',3)
    end
    xlim([0.000 0.1])

    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',2)
    % shaded_region(whisk_resp_win)
    
    title([group_names{a} 'L5 Diff PSTH'])
    fixplot
    
end
subplot_equal_y
