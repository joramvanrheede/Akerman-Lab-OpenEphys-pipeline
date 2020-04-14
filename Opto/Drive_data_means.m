
%% This script will load the preprocessed data

% To do - Include option to detect threshold value in opto_pulse_function
% To do - Automatically select 'nearest' power on a by experiment basis (if we don't have a power for all experiments)
% Overlay of psths

%% Changing any of these parameters means you need to 'do reload'
do_reload           = false;

% List data folders for the experiment types
data_folders        = {'/Users/Joram/Data/Preprocessed/POM/Drive' '/Users/Joram/Data/Preprocessed/M1/Drive'};

whisk_resp_win      = [0.006 0.025];
psth_bins           = [-0.25:0.001:.2];
artifact_win        = [-0.002 0.006];


%% These parameters can be changed without reloading

% Nicknames for the groups
group_names         = {'POM' 'M1'};

% Response measure to plot / analyse
resp_measure_name 	= 'delta_spike_rates'; % 'delta_spike_rates', 'delta_peak_spike_rates'
target_power     	= 5; % [a number, or 'threshold' to use a threshold]
resp_threshold      = 1;

% Channels
L23_chans           = 1:7;
L4_chans            = 9:15;
L5_chans            = 17:30;

% PSTH visualisation
psth_smoothwin      = [10]; % nr of samples
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
        groups(a).drive_data        = get_opto_drive_data(data_folders{a}, whisk_resp_win, psth_bins);
    end
end


group_number                    = [];
all_psths                       = [];

plot_L23_resp_measure           = [];
plot_L23_spike_density_rates    = [];
plot_L23_control_density_rates  = [];
plot_L23_diff_density_rates     = [];

plot_L4_resp_measure           = [];
plot_L4_spike_density_rates    = [];
plot_L4_control_density_rates  = [];
plot_L4_diff_density_rates     = [];

plot_L5_resp_measure            = [];
plot_L5_spike_density_rates     = [];
plot_L5_control_density_rates   = [];
plot_L5_diff_density_rates      = [];

for a = 1:length(groups)
    drive_data     = groups(a).drive_data;
    
    L23_psths   = [];
    L4_psths    = [];
    L5_psths    = [];
    for b = 1:length(drive_data)
        
        %% Select the experiment file with the right opto power
        opto_powers                     = [drive_data(b).experiment(:).opto_power];
        power_diffs                     = target_power - opto_powers;
        abs_power_diffs                 = abs(power_diffs);
        min_power_diff                  = min(abs_power_diffs);
        q_min_power_diff                = abs_power_diffs == min_power_diff;
        closest_power_ind               = find(q_min_power_diff,1,'first');
        best_power                      = opto_powers(closest_power_ind);
        
        %% Select the appropriate time difference between opto stim and whisk stim
        experiment                      = drive_data(b).experiment(closest_power_ind);
        
        
        %% Response measure for comparison plots
        
%         control_spike_rate          	= experiment.control_whisk_spike_rates;
%         spont_spike_std              	= experiment.spont_spike_std;
%         
%         experiment.data_folder
%         
%         L23_control_spike_rate       	= mean(control_spike_rate(L23_chans))
%         L23_spont_spike_std             = mean(spont_spike_std(L23_chans))            
%         L23_resp                        = L23_control_spike_rate > resp_threshold * L23_spont_spike_std
%         
%         L5_control_spike_rate       	= mean(control_spike_rate(L5_chans))
%         L5_spont_spike_std              = mean(spont_spike_std(L5_chans))             
%         L5_resp                         = L5_control_spike_rate > resp_threshold * L5_spont_spike_std
%         
        resp_measure                    = experiment.(resp_measure_name);
        
        
        L23_resp_measure                = mean(resp_measure(L23_chans));
        L4_resp_measure                 = mean(resp_measure(L4_chans));
        L5_resp_measure                 = mean(resp_measure(L5_chans));
        
        plot_L23_resp_measure           = [plot_L23_resp_measure; L23_resp_measure];
        plot_L4_resp_measure            = [plot_L4_resp_measure; L4_resp_measure];
        plot_L5_resp_measure            = [plot_L5_resp_measure; L5_resp_measure];
        
        group_number                    = [group_number; a];
        
        
        %% 
        q_artifact                      = psth_bins > artifact_win(1) & psth_bins < artifact_win(2);
        
        %%
        L23_spike_density_rates         = mean(experiment.spike_density_rates(L23_chans,:,1));
        L23_control_density_rates       = mean(experiment.spike_density_rates(L23_chans,:,2));
        L23_spike_density_rates(q_artifact)     = 0;
        L23_control_density_rates(q_artifact)   = 0;
        L23_diff_density_rates          = L23_spike_density_rates - L23_control_density_rates;
        
        plot_L23_spike_density_rates    = [plot_L23_spike_density_rates L23_spike_density_rates(:)];
        plot_L23_control_density_rates  = [plot_L23_control_density_rates L23_control_density_rates(:)];
        plot_L23_diff_density_rates     = [plot_L23_diff_density_rates L23_diff_density_rates(:)];
        
        %% 
        
        L4_spike_density_rates         = mean(experiment.spike_density_rates(L4_chans,:,1));
        L4_control_density_rates       = mean(experiment.spike_density_rates(L4_chans,:,2));
        L4_spike_density_rates(q_artifact)     = 0;
        L4_control_density_rates(q_artifact)   = 0;
        L4_diff_density_rates          = L4_spike_density_rates - L4_control_density_rates;
        
        plot_L4_spike_density_rates    = [plot_L4_spike_density_rates L4_spike_density_rates(:)];
        plot_L4_control_density_rates  = [plot_L4_control_density_rates L4_control_density_rates(:)];
        plot_L4_diff_density_rates     = [plot_L4_diff_density_rates L4_diff_density_rates(:)];
        
        %%
        
        L5_spike_density_rates         = mean(experiment.spike_density_rates(L5_chans,:,1));
        L5_control_density_rates       = mean(experiment.spike_density_rates(L5_chans,:,2));
        L5_spike_density_rates(q_artifact)     = 0;
        L5_control_density_rates(q_artifact)   = 0;
        L5_diff_density_rates          = L5_spike_density_rates - L5_control_density_rates;
        
        plot_L5_spike_density_rates    = [plot_L5_spike_density_rates L5_spike_density_rates(:)];
        plot_L5_control_density_rates  = [plot_L5_control_density_rates L5_control_density_rates(:)];
        plot_L5_diff_density_rates     = [plot_L5_diff_density_rates L5_diff_density_rates(:)];

        
    end

end

% Make resp measure name a string suitable for plots
resp_measure_name(ismember(resp_measure_name,'_')) = ' ';
resp_measure_name(1) = upper(resp_measure_name(1));

x_lims      = [psth_bins(1) psth_bins(end)];

%% Compare data points and means for L23
figure
beeswarmplot(plot_L23_resp_measure,group_number,group_names)
hold on
plot(xlim,[0 0],'k:','LineWidth',2)
title(['L23'])
ylabel(resp_measure_name)
fixplot

%% Compare data points and means for L4
figure
beeswarmplot(plot_L4_resp_measure,group_number,group_names)
hold on
plot(xlim,[0 0],'k:','LineWidth',2)
title(['L4'])
ylabel(resp_measure_name)
fixplot

%% Compare data points and means for L5
figure
beeswarmplot(plot_L5_resp_measure,group_number,group_names)
hold on
plot(xlim,[0 0],'k:','LineWidth',2)
title(['L5'])
ylabel(resp_measure_name)
fixplot


%% OVERLAY PSTHS

psth_overlay_plot(plot_L23_control_density_rates,plot_L23_spike_density_rates,psth_bins,group_number,psth_smoothwin,'L23 PSTHs',group_names)
psth_overlay_plot(plot_L4_control_density_rates,plot_L4_spike_density_rates,psth_bins,group_number,psth_smoothwin,'L4 PSTHs',group_names)
psth_overlay_plot(plot_L5_control_density_rates,plot_L5_spike_density_rates,psth_bins,group_number,psth_smoothwin,'L5 PSTHs',group_names)


%% DIFF PSTHS

diff_psth(plot_L23_diff_density_rates, group_number, psth_smoothwin, psth_bins, group_names,'L23 Difference PSTH (Test - control)')
diff_psth(plot_L4_diff_density_rates, group_number, psth_smoothwin, psth_bins, group_names,'L4 Difference PSTH (Test - control)')
diff_psth(plot_L5_diff_density_rates, group_number, psth_smoothwin, psth_bins, group_names,'L5 Difference PSTH (Test - control)')



%% Functions from here


function psth_overlay_plot(control_psth_data,test_psth_data,psth_bins,group_number,psth_smoothwin,label,group_names)

figure
set(gcf,'Units','Normalized','Position',[.1 .4 .8 .4])
for a = 1:max(group_number)
    
    q_group         = group_number == a;
    
    control_psths   = control_psth_data(:,q_group);
    smooth_control  = smoothdata(control_psths,1,'gaussian',psth_smoothwin);
    
    test_psths     	= test_psth_data(:,q_group);
    smooth_test     = smoothdata(test_psths,1,'gaussian',psth_smoothwin);
    
    %%
    subplot(max(group_number),2,a)
    for b = 1:size(smooth_control,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_control(:,b),'FaceColor',[0 0 0],'FaceAlpha',.25,'BarWidth',1);
        hold on
    end
    xlim([-0.05 0.2])
    plot(psth_bins(1:end-1),mean(smooth_control,2),'k-','LineWidth',3)
    title([group_names{a} ' ' label ' - CONTROL'])
    fixplot
    
    %%
    subplot(max(group_number),2,max(group_number)+a)
    for b = 1:size(smooth_test,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_test(:,b),'FaceColor',[0 0 0],'FaceAlpha',.25,'BarWidth',1);
        hold on
    end
    xlim([-0.05 0.2])
    plot(psth_bins(1:end-1),mean(smooth_test,2),'k-','LineWidth',3)
    title([group_names{a} ' ' label ' - TEST'])
    fixplot
    
    
end
subplot_equal_y
end

%%
function diff_psth(psth_diff_data, group_number, psth_smoothwin, psth_bins, group_names, label)
n_groups = max(group_number);
figure
set(gcf,'Units','Normalized','Position',[.3 .4 .4 .4])
for a = 1:n_groups
    subplot(1,n_groups,a)
    q_group         = group_number == a;
    this_psth       = psth_diff_data(:,q_group); 
    
    smooth_psths    = smoothdata(this_psth,1,'gaussian',psth_smoothwin);
    
    for b = 1:size(smooth_psths,2)
        bar_h           = bar(psth_bins(1:end-1),smooth_psths(:,b),'FaceColor',[0 0 0],'FaceAlpha',.4,'BarWidth',1);
        hold on
    end
    xlim([0.000 0.2])

    plot(psth_bins(1:end-1),mean(smooth_psths,2),'k-','LineWidth',2)
    
    title([group_names{a} ' ' label])
    fixplot
    
end
subplot_equal_y
end
