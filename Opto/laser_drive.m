% Laser_Pulse analysis


data_dir            = ['/Volumes/Akermanlab/Joram/RJB/Drive/2019_01_24'];

% data_file       = ['2019_01_22-16-Laser_pulse'];

hist_binsize        = [0.005];
rate_window         = [0.020];

resp_win            = [0.050 0.100];
spont_win           = [0 .9];
win_margin          = [0.2];

smoothwin           = 7;

qreload             = 1;

target_channels     = [1:16];

%% Make a 'sum' histogram picture of the summed post-LED window and post whisk window

if qreload
    
    folder_files    = dir([data_dir filesep '*.mat']); % read .mat files from folder
    
    laser_power     = [];
    laser_psth      = [];
    control_psth    = [];
    for a = 1:length(folder_files)
        
        this_file           = folder_files(a).name;
        full_filenm         = fullfile(data_dir, this_file);
        
        disp(['Loading ' this_file]);
        load(full_filenm); % this should load variables 'ephys_data' and 'parameters'
        
        ephys_data          = ephys_data_psths(ephys_data, hist_binsize, rate_window);
        
        laser_power(a)      = ephys_data.parameters.LED_power;              % What power was the LASER/LED set to?
        
        control_laser_onset = max(ephys_data.condition_values(:,1));      	% First column of ephys_data.condition_values contains LED/LASER onset time; control trials are distinguished by their later LASER onset time
        test_laser_onset    = min(ephys_data.condition_values(:,1));      	% First column of ephys_data.condition_values contains LED/LASER onset time; control trials are distinguished by their later LASER onset time
        
        laser_duration      = mean(ephys_data.condition_values(:,3));       % All conditions in this protocol should have the same laser duration
        
        time_min            = control_laser_onset - win_margin;
        time_max            = time_min + laser_duration + win_margin;
        laser_control_win   = [time_min time_max];
        
        stim_win            = laser_control_win - (control_laser_onset - test_laser_onset);
        
        
        for b = 1:length(ephys_data.conditions)
            
            time_bins   = ephys_data.conditions(b).psth_bins;
            
            q_control_win       = time_bins > min(laser_control_win) & time_bins < max(laser_control_win);
            q_test_win          = time_bins > min(stim_win) & time_bins < max(stim_win);
            
            test_psth(a,b,:)    = smooth(mean(ephys_data.conditions(b).expt_psths(target_channels,:)),smoothwin);
            
            control_psth(a,b,:) = smooth(mean(ephys_data.conditions(b).expt_psths(target_channels,:)),smoothwin);
            
        end
    end
end

% 
% laser_plus_whisk
% laser_with_whisk
% 
% ephys_data.condition_values(:,1)

laser_psth = test_psth;


figure
plot(squeeze(laser_psth(1,[1],:))','k-','LineWidth',2)
hold on
plot(squeeze(laser_psth(1,[2],:))','r-','LineWidth',2)
hold off
% figure
% plot(squeeze(laser_psth(1,[2],:))','k-','LineWidth',2)
% hold on
% plot(squeeze(laser_psth(1,[4],:))','r-','LineWidth',2)
% hold off
clunk

for c = 1:size(laser_psth,2)
    whisk_onset     = ephys_data.conditions(c).LED_onset;
    pulse_duration  = ephys_data.conditions(c).LED_duration;
    psth_bins       = ephys_data.conditions(c).psth_bins - pulse_onset;
    
    q_pulse      	= psth_bins >= 0 & psth_bins <= pulse_duration;
    q_respwin       = psth_bins >= resp_win(1) & psth_bins <= resp_win(2);
    q_spontwin   	= psth_bins + pulse_onset >= spont_win(1) & psth_bins + pulse_onset <= spont_win(2);
    
    laser_psth(:,q_pulse)   = NaN;
    
    uniq_powers     = unique(laser_power);
    n_powers        = length(uniq_powers);
    
    mean_psth_by_power  = [];
    spike_sum_by_power  = [];
    spont_sum_by_power  = [];
    
    this_power                  = uniq_powers(b);
    q_power                     = laser_power == this_power;
    
    these_psths                 = laser_psth(q_power,:);
    
    mean_psth_by_power(b,:)     = smooth(mean(these_psths,1),smoothwin);
    
    spike_sum_by_power(b)       = mean(mean(these_psths(:,q_respwin),1));
    
    spont_sum_by_power(b)       = mean(mean(these_psths(:,q_spontwin),1));
    
end




% plot(psth_bins,mean_psth_by_power','LineWidth',2)
% hold on, plot(psth_bins,mean_psth_by_power(1,:),'k-','LineWidth',2)
% title('PSTH of spikes following 3ms LASER pulse')
% xlabel('Post-pulse Time (s)')
% ylabel('Spike count in time bin')
% set(gca,'FontName','Helvetica','FontSize',14)
% 
% %legend(mat2cell(uniq_powers));
% 
% figure
% plot(uniq_powers,spike_sum_by_power - spont_sum_by_power,'k.-','LineWidth',2,'MarkerSize',20)
% title('LASER Power setting vs binned spiking response','FontName','Helvetica','FontSize',20)
% xlabel('LASER Power setting','FontName','Helvetica','FontSize',16)
% ylabel('Spiking response','FontName','Helvetica','FontSize',16)
% set(gca,'FontName','Helvetica','FontSize',14)


