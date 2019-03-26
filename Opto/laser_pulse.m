% Laser_Pulse analysis

% Where are the preprocessed data?
% AVK:
% POM
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Laser_pulse/2018_12_03
% 
% 
% S1
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV cS1/Laser_pulse/2019_02_07_P1 --> nothing
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV cS1/Laser_pulse/2019_02_07_P2 --> barely anything
% 
% 
% 


data_dir        = '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Laser_pulse/2018_12_03';

% data_file       = ['2019_01_22-16-Laser_pulse'];

% Do LFP

hist_binsize    = [0.001];
rate_window     = [0.020];

resp_win        = [0.0035 0.020]; % [0.001 0.003];
spont_win       = [0 .9];

plot_win        = [-0.1 0.2];

smoothwin       = 1;

qreload         = 1;

%% 
close all


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
figure
plot(psth_bins,mean_psth_by_power')
title('PSTH of spikes following 3ms LASER pulse')
xlabel('Post-pulse Time (s)')
ylabel('Spike count in time bin')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')
xlim(plot_win)


%legend(mat2cell(uniq_powers));

figure
plot(uniq_powers,spike_sum_by_power - spont_sum_by_power,'k.-','MarkerSize',20)
title('LASER Power setting vs binned spiking response')
xlabel('LASER Power setting')
ylabel('Spiking response')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')


