% opto_pulse_multi

% List of experiments to be included?
data_folder         = '/Volumes/Akermanlab/Joram/Preprocessed data/AVK Cortex/Laser_pulse';

do_reload           = true;

psth_smoothing      = 11;

% Channel groups / work out a way to do this automatically from data, for each individual experiment. In metadata file?
L23_chans           = 1:8;
L4_chans            = 9:14;
L5a_chans           = 16:23;
L5b_chans           = 16:30;
all_chans           = 1:32;

opto_early_win      = [0.007 0.200];
opto_late_win       = [0.050 0.150];

psth_bins           = [-0.25:0.002:.5];
artifact_win        = [-0.002 0.007];

psth_power          = 5;

save_figs           = false;
fig_save_dir       	= '/Users/Joram/Dropbox/Akerman Postdoc/Figures/2019_09_25_LabM';
expt_name           = 'Laser_pulse_';

%% Code execution starts here

opto_wins           = [opto_early_win; opto_late_win];

%% Function: load pulse_data (with multiple time win inputs)
if do_reload
    pulse_data          = get_opto_pulse_data(data_folder, opto_wins, psth_bins, artifact_win);
end


for a = 1:2
    fig_h(a) = figure(a);
    set(fig_h(a),'Units','Normalized','Position',[.05 .3 .9 .6])
end

q_early_win     = psth_bins > opto_early_win(1) & psth_bins <= opto_early_win(2);
q_late_win      = psth_bins > opto_late_win(1) & psth_bins <= opto_late_win(2);

for a = 1:length(pulse_data)
    
    max_power_spike_rate        = pulse_data(a).delta_spike_rate(:,end);
    channel_spike_profile(:,a) 	= smooth(max_power_spike_rate,3);
    [max_resp(a), max_chan(a)]  = max(channel_spike_profile(:,a));
    
    
    early_win_size          = opto_early_win(2) - opto_early_win(1);
    late_win_size           = opto_late_win(2) - opto_late_win(1);
    
    L23_early_spikes        = mean(pulse_data(a).delta_spike_rate(L23_chans,:,1),1);
    L23_late_spikes         = mean(pulse_data(a).delta_spike_rate(L23_chans,:,2),1);

    L4_early_spikes         = mean(pulse_data(a).delta_spike_rate(L4_chans,:,1),1);
    L4_late_spikes          = mean(pulse_data(a).delta_spike_rate(L4_chans,:,2),1);

    L5a_early_spikes        = mean(pulse_data(a).delta_spike_rate(L5a_chans,:,1),1);
    L5a_late_spikes         = mean(pulse_data(a).delta_spike_rate(L5a_chans,:,2),1);

    L5b_early_spikes        = mean(pulse_data(a).delta_spike_rate(L5b_chans,:,1),1);
    L5b_late_spikes         = mean(pulse_data(a).delta_spike_rate(L5b_chans,:,2),1);

    %% 
    figure(1)
    
    subplot(1,2,1)
    plot(pulse_data(a).opto_power, L23_early_spikes,'LineWidth',2)
    title('L2/3 early spike rate (- spont)')
    ylabel('Delta spike rate (Hz)')
    xlabel('Laser power')
    axis tight
    xlim([4.5 10.5])
    hold on
    fixplot
    
%     subplot(1,4,2)
%     plot(pulse_data(a).opto_power, L4_early_spikes,'LineWidth',2)
%     title('L4 early spikes')
%     xlim([4.5 10.5])
%     hold on
%     fixplot
%     
%     subplot(1,4,3)
%     plot(pulse_data(a).opto_power, L5a_early_spikes,'LineWidth',2)
%     title('L5a early spikes')
%     xlim([4.5 10.5])
%     hold on
%     fixplot
    
    subplot(1,2,2)
    plot(pulse_data(a).opto_power, L5b_early_spikes,'LineWidth',2)
    title('L5 early spike rate (- spont)')
    ylabel('Delta spike rate (Hz)')
    xlabel('Laser power')
    axis tight
    xlim([4.5 10.5])
    hold on
    fixplot
    
    %% 
    figure(2)
    
    subplot(1,2,1)
    plot(pulse_data(a).opto_power, L23_late_spikes,'LineWidth',2)
    title('L2/3 late spike rate (- spont)')
    ylabel('Delta spike rate (Hz)')
    xlabel('Laser power')
    axis tight
    xlim([4.5 10.5])
    hold on
    fixplot
    
%     subplot(1,4,2)
%     plot(pulse_data(a).opto_power, L4_late_spikes,'LineWidth',2)
%     title('L4 late spikes')
%     xlim([4.5 10.5])
%     hold on
%     fixplot
%     
%     subplot(1,4,3)
%     plot(pulse_data(a).opto_power, L5a_late_spikes,'LineWidth',2)
%     title('L5a late spikes')
%     xlim([4.5 10.5])
%     hold on
%     fixplot
    
    subplot(1,2,2)
    plot(pulse_data(a).opto_power, L5b_late_spikes,'LineWidth',2)
    title('L5 late spike rate (- spont)')
    ylabel('Delta spike rate (Hz)')
    xlabel('Laser power')
    axis tight
    xlim([4.5 10.5])
    hold on
    fixplot
    
    
end

opto_powers         = [pulse_data.opto_power];
early_spike_rates   = [pulse_data.delta_spike_rate];
late_spike_rates    = [pulse_data.delta_spike_rate];

% Remove more finely spaced data
uniq_opto_powers    = unique(ceil(opto_powers));
uniq_opto_powers(uniq_opto_powers < 5) = [];


mean_early_L23_rate = [];
mean_early_L4_rate  = [];
mean_early_L5a_rate = [];
mean_early_L5b_rate = [];

mean_late_L23_rate  = [];
mean_late_L4_rate   = [];
mean_late_L5a_rate  = [];
mean_late_L5b_rate  = [];

serr_early_L23_rate	= [];
serr_early_L4_rate  = [];
serr_early_L5a_rate	= [];
serr_early_L5b_rate	= [];

serr_late_L23_rate  = [];
serr_late_L4_rate   = [];
serr_late_L5a_rate  = [];
serr_late_L5b_rate  = [];

for a = 1:length(uniq_opto_powers)
    this_power  = uniq_opto_powers(a);
    q_power     = opto_powers == this_power;
    
    these_early_L23_rates       = mean(early_spike_rates(L23_chans,q_power),1);
    these_early_L4_rates        = mean(early_spike_rates(L4_chans,q_power),1);
    these_early_L5a_rates       = mean(early_spike_rates(L5a_chans,q_power),1);
    these_early_L5b_rates       = mean(early_spike_rates(L5b_chans,q_power),1);
    
    these_late_L23_rates        = mean(late_spike_rates(L23_chans,q_power),1);
    these_late_L4_rates         = mean(late_spike_rates(L4_chans,q_power),1);
    these_late_L5a_rates        = mean(late_spike_rates(L5a_chans,q_power),1);
    these_late_L5b_rates        = mean(late_spike_rates(L5b_chans,q_power),1);
    
    % Early
    mean_early_L23_rate(a)      = mean(these_early_L23_rates);
    serr_early_L23_rate(a)      = serr(these_early_L23_rates);
    
    mean_early_L4_rate(a)       = mean(these_early_L4_rates);
    serr_early_L4_rate(a)       = serr(these_early_L4_rates);
    
    mean_early_L5a_rate(a)      = mean(these_early_L5a_rates);
    serr_early_L5a_rate(a)      = serr(these_early_L5a_rates);
    
    mean_early_L5b_rate(a)      = mean(these_early_L5b_rates);
    serr_early_L5b_rate(a)      = serr(these_early_L5b_rates);
    
    % Late
    mean_late_L23_rate(a)       = mean(these_late_L23_rates);
    serr_late_L23_rate(a)       = serr(these_late_L23_rates);
    
    mean_late_L4_rate(a)        = mean(these_late_L4_rates);
    serr_late_L4_rate(a)        = serr(these_late_L4_rates);
    
    mean_late_L5a_rate(a)       = mean(these_late_L5a_rates);
    serr_late_L5a_rate(a)       = serr(these_late_L5a_rates);
    
    mean_late_L5b_rate(a)       = mean(these_late_L5b_rates);
    serr_late_L5b_rate(a)       = serr(these_late_L5b_rates);

end

% Early
figure(1)
subplot(1,2,1)
errorbar(uniq_opto_powers, mean_early_L23_rate, serr_early_L23_rate,'Color',[0 0 0],'LineWidth',3)
% subplot(1,4,2)
% errorbar(uniq_opto_powers, mean_early_L4_rate, serr_early_L4_rate,'Color',[0 0 0],'LineWidth',3)
% subplot(1,4,3)
% errorbar(uniq_opto_powers, mean_early_L5a_rate, serr_early_L5a_rate,'Color',[0 0 0],'LineWidth',3)
subplot(1,2,2)
errorbar(uniq_opto_powers, mean_early_L5b_rate, serr_early_L5b_rate,'Color',[0 0 0],'LineWidth',3)



% Late
figure(2)
subplot(1,2,1)
errorbar(uniq_opto_powers, mean_late_L23_rate, serr_late_L23_rate,'Color',[0 0 0],'LineWidth',3)
% subplot(1,4,2)
% errorbar(uniq_opto_powers, mean_late_L4_rate, serr_late_L4_rate,'Color',[0 0 0],'LineWidth',3)
% subplot(1,4,3)
% errorbar(uniq_opto_powers, mean_late_L5a_rate, serr_late_L5a_rate,'Color',[0 0 0],'LineWidth',3)
subplot(1,2,2)
errorbar(uniq_opto_powers, mean_late_L5b_rate, serr_late_L5b_rate,'Color',[0 0 0],'LineWidth',3)


norm_channel_spike_profile   = bsxfun(@rdivide, channel_spike_profile, max(channel_spike_profile));

[[1:length(max_chan)]' max_chan(:) max_resp(:)]


L23_hist_stack  = [];
L4_hist_stack   = [];
L5a_hist_stack  = [];
L5b_hist_stack  = [];
for a = 1:length(pulse_data)
    these_opto_powers	= pulse_data(a).opto_power;
    spike_counts        = pulse_data(a).density_counts;
    
    q_power             = these_opto_powers == psth_power;
    
    L23_spike_counts    = mean(spike_counts(L23_chans,:,q_power));
    L4_spike_counts     = mean(spike_counts(L4_chans,:,q_power));
    L5a_spike_counts    = mean(spike_counts(L5a_chans,:,q_power));
    L5b_spike_counts    = mean(spike_counts(L5b_chans,:,q_power));
    
%     L23_hist_smooth         = smooth(L23_spike_counts,psth_smoothing)';
%     L23_hist_smooth_norm    = 
    
    L23_hist_stack      = [L23_hist_stack; smooth(L23_spike_counts,psth_smoothing)'];
    L4_hist_stack       = [L4_hist_stack; smooth(L4_spike_counts,psth_smoothing)'];
    L5a_hist_stack      = [L5a_hist_stack; smooth(L5a_spike_counts,psth_smoothing)'];
    L5b_hist_stack      = [L5b_hist_stack; smooth(L5b_spike_counts,psth_smoothing)'];
end
figure
set(gcf,'Units','Normalized','Position',[.05 .3 .9 .6])
psth_x  = psth_bins(1:end-1);
subplot(1,2,1)
plot(psth_x,L23_hist_stack,'LineWidth',2)
fixplot
xlabel('Time (s)')
ylabel('Spike count')
title('L2/3 PSTH')

% subplot(1,4,2)
% plot(psth_x,L4_hist_stack,'LineWidth',2)
% fixplot
% subplot(1,4,3)
% plot(psth_x,L5a_hist_stack,'LineWidth',2)
% fixplot
subplot(1,2,2)
plot(psth_x,L5b_hist_stack,'LineWidth',2)
fixplot
xlabel('Time (s)')
ylabel('Spike count')
title('L5 PSTH')


%%
if save_figs
    figure(1)
    print(gcf,'-r300','-dpng',[fig_save_dir filesep expt_name '_RateXPower_early'])
    figure(2)
    print(gcf,'-r300','-dpng',[fig_save_dir filesep expt_name '_RateXPower_late'])
    figure(3)
    print(gcf,'-r300','-dpng',[fig_save_dir filesep expt_name '_PSTHs'])
end





