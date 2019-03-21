% Laser_Pulse analysis


data_dir        = '/Volumes/Akermanlab/Joram/RJB/Laser_Pulse/2019_01_24';

% data_file       = ['2019_01_22-16-Laser_pulse'];

% Do LFP

hist_binsize    = [0.001];
rate_window     = [0.020];

resp_win        = [0.003 0.500];
spont_win       = [0 .9];

smoothwin       = 5;

qreload         = 1;

if qreload
    
    folder_files    = dir([data_dir filesep '*.mat']); % read .mat files from folder
    
    laser_power     = [];
    laser_psth      = [];
    for a = 1%1:length(folder_files)
        
        this_file       = folder_files(a).name;
        full_filenm     = fullfile(data_dir, this_file);
        
        disp(['Loading ' this_file]);
        load(full_filenm); % this should load variables 'ephys_data' and 'parameters'
        
        ephys_data      = ephys_data_psths(ephys_data, hist_binsize, rate_window);
        
        laser_power(a)  = ephys_data.parameters.LED_power;
        
        
        laser_psth(a,:) = mean(ephys_data.conditions(1).expt_psths);
    end
    
end

pulse_onset     = ephys_data.conditions(1).LED_onset;
pulse_duration  = ephys_data.conditions(1).LED_duration;
psth_bins       = ephys_data.conditions(1).psth_bins - pulse_onset;

q_pulse      	= psth_bins >= 0 & psth_bins <= pulse_duration;
q_respwin       = psth_bins >= resp_win(1) & psth_bins <= resp_win(2);
q_spontwin   	= psth_bins + pulse_onset >= spont_win(1) & psth_bins + pulse_onset <= spont_win(2);

% laser_psth(:,q_pulse)   = NaN;

uniq_powers     = unique(laser_power);
n_powers        = length(uniq_powers);


%% Loop to make PSTH for each laser power
mean_psth_by_power  = [];
spike_sum_by_power  = [];
spont_sum_by_power  = [];
for b = 1:n_powers
    this_power                  = uniq_powers(b);
    q_power                     = laser_power == this_power; 
    
    these_psths                 = laser_psth(q_power,:);
    
    mean_psth_by_power(b,:)     = smooth(mean(these_psths,1),smoothwin);
    
    spike_sum_by_power(b)       = mean(mean(these_psths(:,q_respwin),1));
    
    spont_sum_by_power(b)       = mean(mean(these_psths(:,q_spontwin),1));
end

% plot all psths
plot(psth_bins,mean_psth_by_power','LineWidth',2)
title('PSTH of spikes following 3ms LASER pulse')
xlabel('Post-pulse Time (s)')
ylabel('Spike count in time bin')
set(gca,'FontName','Helvetica','FontSize',14)




%legend(mat2cell(uniq_powers));

figure
plot(uniq_powers,spike_sum_by_power - spont_sum_by_power,'k.-','LineWidth',2,'MarkerSize',20)
title('LASER Power setting vs binned spiking response','FontName','Helvetica','FontSize',20)
xlabel('LASER Power setting','FontName','Helvetica','FontSize',16)
ylabel('Spiking response','FontName','Helvetica','FontSize',16)
set(gca,'FontName','Helvetica','FontSize',14)


