% opto_power_single_file

% Where are the preprocessed data?
% AVK:
% POM OLD 2018_09_10
% '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/LED power/2018_09_10/2018_09_10-1-LED power.mat' High pwr --> interesting disconnect between early and late response window
% '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/LED power/2018_09_10/2018_09_10-2-LED power.mat'

% 2019_09_12
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/LED power/2018_09_12/2018_09_12-3-LED power.mat

data_file        = '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/LED power/2018_09_10/2018_09_10-1-LED power.mat';


psth_bin_size  	= [0.001];

resp_win        = [0 0.04]% [0.05 0.500]% [0 0.04] % [0 0.04] % [0.05 0.500]; 
spont_win       = [-1 -0.1];

plot_win        = [-0.5 1.5];

psth_win        = [-1 2];

channels        = [1:16];

smoothing       = 9;

qreload         = 1;

%%

close all

if qreload
    load(data_file)
end

opto_powers     = [ephys_data.conditions(:).LED_power]';

opto_powers(isnan(opto_powers)) = 2.5;

n_powers        = length(opto_powers);

%% Loop to make PSTH for each laser power
psth_by_power       = [];
spike_sum_by_power  = [];
spont_sum_by_power  = [];
for b = 1:n_powers
    this_power                  = opto_powers(b)
    q_power                     = opto_powers == this_power; 
    
    opto_onset              	= ephys_data.conditions(q_power).LED_onset;
    
    spikes                   	= ephys_data.conditions(q_power).spikes - opto_onset;
    
    n_trials                    = ephys_data.conditions(q_power).n_trials
    
    figure
    [plot_handle, psth_counts, psth_bins]   = psth(spikes, psth_bin_size, psth_win);
    close(gcf)
    
    q_resp_win                  = psth_bins > resp_win(1) & psth_bins <= resp_win(2);
    q_spont_win                 = psth_bins > spont_win(1) & psth_bins <= spont_win(2);
    
    psth_by_power(b,:)          = psth_counts / n_trials / length(channels);
    
    spike_sum_by_power(b)       = mean(psth_by_power(b,q_resp_win),2);
    
    spont_sum_by_power(b)       = mean(psth_by_power(b,q_spont_win),2);
    
    psth_by_power(b,:)          = smooth(psth_by_power(b,:),smoothing); % smooth for plotting only
end

% plot all psths
figure(1)
plot(psth_bins(1:end-1),psth_by_power')
title('PSTH of spikes following 3ms LASER pulse')
xlabel('Post-pulse Time (s)')
ylabel('Spike count in time bin')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')
xlim(plot_win)
yzero
scf

%legend(mat2cell(uniq_powers));

[opto_powers, sortinds] = sort(opto_powers);

figure(2)
plot(opto_powers,spike_sum_by_power(sortinds) - spont_sum_by_power(sortinds),'k.-','MarkerSize',20)
title('LASER Power setting vs binned spiking response')
xlabel('LASER Power setting')
ylabel('Spiking response')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','TickDir','out','box','off')
yzero
scf
