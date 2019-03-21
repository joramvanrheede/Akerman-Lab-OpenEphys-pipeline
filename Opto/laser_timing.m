% Laser_Pulse analysis


data_dir        = ['/Volumes/Akermanlab/Joram/RBSN/Timing/2019_02_15'];

% data_file       = ['2019_01_22-16-Laser_pulse'];

hist_binsize    = [0.001];
rate_window     = [0.020];

resp_win        = [0.050 0.100];
spont_win       = [0 .9];

smoothwin       = 7;

qreload         = 1;

if qreload
    
    folder_files    = dir([data_dir filesep '*.mat']); % read .mat files from folder
    
    laser_power     = [];
    laser_psth      = [];
    for a = 1:length(folder_files)
        
        this_file       = folder_files(a).name;
        full_filenm     = fullfile(data_dir, this_file);
        
        disp(['Loading ' this_file]);
        load(full_filenm); % this should load variables 'ephys_data' and 'parameters'
        
        ephys_data      = ephys_data_psths(ephys_data, hist_binsize, rate_window);
        
        laser_power(a)  = ephys_data.parameters.LED_power;
        
        for b = 1:length(ephys_data.conditions)
            
            laser_psth(a,b,:) = smooth(mean(ephys_data.conditions(b).expt_psths),smoothwin);
            
        end
    end
end

subplot(1,2,1)
plot(squeeze(laser_psth(1,[2:2:12],:))')
xlim([950 1250])
subplot(1,2,2)
plot(squeeze(laser_psth(1,[1:2:11],:))','LineWidth',2)
xlim([950 1250])

clunk

for c = 1:size(laser_psth,2)
    whisk_onset     = ephys_data.conditions(c).whisk_onset;
    pulse_onset     = ephys_data.conditions(c).LED_onset;
    
    pulse_duration  = ephys_data.conditions(c).LED_duration;
    
    whisk_bins     	= ephys_data.conditions(c).psth_bins - whisk_onset;
    laser_bins      = ephys_data.conditions(c).psth_bins - pulse_onset;
    
    q_pulse      	= laser_bins >= 0 & laser_bins <= pulse_duration;
    q_respwin       = whisk_bins >= resp_win(1) & whisk_bins <= resp_win(2);
    q_spontwin   	= whisk_bins + whisk_onset >= spont_win(1) & whisk_bins + whisk_onset <= spont_win(2);
    
    laser_psth(:,q_pulse)   = NaN;
    
    uniq_powers     = unique(laser_power);
    n_powers        = length(uniq_powers);
    
    mean_psth_by_power  = [];
    spike_sum_by_power  = [];
    spont_sum_by_power  = [];
    
    this_power                  = uniq_powers(c);
    q_power                     = laser_power == this_power;
    
    these_psths                 = laser_psth(q_power,:);
    
    mean_psth_by_power(c,:)     = smooth(mean(these_psths,1),smoothwin);
    
    spike_sum_by_power(c)       = mean(mean(these_psths(:,q_respwin),1));
    
    spont_sum_by_power(c)       = mean(mean(these_psths(:,q_spontwin),1));
    
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


