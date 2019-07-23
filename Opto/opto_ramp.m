% opto_plus_whisk

data_file           = '/Volumes/Akermanlab/Joram/Preprocessed data/AVK POM/LED_ramp/2019_07_19/2019_07_19-36-LED_ramp.mat';

channels            = 1:16;
trials              = 'all';

plot_lims           = [0 18]; % Relative to whisk onset
bin_size            = [0.1];

cond_names          = {'Ramp only'};

spike_resp_win      = [0.006 0.030];
opto_resp_win       = [0.150 0.300];
spont_resp_win    	= [0 0.7];

LFP_P2P_win         = [0.005 0.15];
LFP_freq_bins       = [0:1:150];
LFP_freq            = 1000;

opto_LFP_P2P_win    = [0.005 0.2];

clim_perc           = 0.5;

expt_name           = 'MAX CHRONOS-OFF';
fig_save_dir        = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/2019_05_29';

save_figs           = true;

%%
close all

load(data_file)

if ischar(trials) && strcmpi(trials,'all')
    trials              = 1:size(ephys_data.conditions(1).spikes,2);
end

cond_data   = ephys_data.conditions;


mean_LFP_peaks          = [];
serr_LFP_peaks          = [];

mean_spike_counts       = [];
serr_spike_counts       = [];

all_spike_counts        = [];
all_LFP_peaks           = [];

all_opto_spike_counts   = [];
all_opto_LFP_peaks      = [];

all_spont_spike_counts  = [];

for a = 1:length(cond_data)
    spikes                  = cond_data(a).spikes(channels,trials,:) - cond_data(a).whisk_onset;
    spike_bins              = plot_lims(1):bin_size:plot_lims(2);
    
    LFPs                    = cond_data(a).LFP_trace(channels,trials,:);
    LFP_timestamps          = [1:length(LFPs)]/1000 - cond_data(a).whisk_onset;
    
    q_LFP_resp              = LFP_timestamps > LFP_P2P_win(1) & LFP_timestamps < LFP_P2P_win(2);
    
    % Quantify LFP
    LFP_neg_peaks           = squeeze(min(LFPs(:,:,q_LFP_resp),[],3));
    LFP_pos_peaks           = max(LFPs(:,:,q_LFP_resp),[],3);
    LFP_p2p                 = LFP_pos_peaks - LFP_neg_peaks;
    
    all_LFP_peaks(:,a)      = mean(LFP_p2p,1);
    
    % Quantify spikes
    spike_counts            = histc(spikes,spike_resp_win,3);
    spike_counts            = squeeze(spike_counts(:,:,1));
    
    all_spike_counts(:,a)   = mean(spike_counts,1);
    
    %% 
    opto_spikes             = cond_data(a).spikes(channels,trials,:) - cond_data(a).LED_onset;
    
    % opto LFP
    opto_LFP_timestamps     = [1:length(LFPs)]/1000 - cond_data(a).LED_onset;
    q_opto_LFP_resp       	= opto_LFP_timestamps > opto_LFP_P2P_win(1) & opto_LFP_timestamps < opto_LFP_P2P_win(2);
    
    % Quantify opto LFP
    opto_LFP_neg_peaks    	= squeeze(min(LFPs(:,:,q_opto_LFP_resp),[],3));
    opto_LFP_pos_peaks   	= max(LFPs(:,:,q_opto_LFP_resp),[],3);
    opto_LFP_p2p         	= opto_LFP_pos_peaks - opto_LFP_neg_peaks;
    
    all_opto_LFP_peaks(:,a) = mean(opto_LFP_p2p,1);
    
    
    % Quantify opto spikes
    opto_spike_counts    	= histc(opto_spikes,opto_resp_win,3);
    opto_spike_counts    	= squeeze(opto_spike_counts(:,:,1));
    
    opto_win_size           = opto_resp_win(2)-opto_resp_win(1);
    whisk_win_size          = spike_resp_win(2)-spike_resp_win(1);
    
    opto_win_ratio          = whisk_win_size / opto_win_size;
    
    all_opto_spike_counts(:,a)   = mean(opto_spike_counts,1) * opto_win_ratio;
    
    spont_spikes            = cond_data(a).spikes(channels,trials,:);
    spont_spike_counts    	= histc(spont_spikes,spont_resp_win,3);
    spont_spike_counts    	= squeeze(spont_spike_counts(:,:,1));
    
    spont_win_size        	= spont_resp_win(2)-spont_resp_win(1);
    spont_win_ratio         = whisk_win_size / spont_win_size;
    
    all_spont_spike_counts(:,a)     = mean(spont_spike_counts,1) * spont_win_ratio;
    
    % Make plots:
    figure(1)
    set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
    subplot(1,length(cond_data),a)
    psth(spikes,bin_size,plot_lims)
    title(cond_names{a})
    fixplot
    
    figure(2)
    set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
    subplot(1,length(cond_data),a)
    [image_handle, count_data] = spike_density_plot(spikes,1,spike_bins);
    ylabel('Channel')
    xlabel('Time (s)')
    color_lims = [0 robust_max(count_data(:),clim_perc)];
    set(gca,'CLim',color_lims)
    colorbar
    title([cond_names{a} ' spike count'])
    
    
    
    figure(3)
    set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
    subplot(1,length(cond_data),a)
    raster_plot(spikes,2)
    xlim(plot_lims)
    ylabel('Trial number')
    xlabel('Time (s)')
    fixplot
    title(cond_names{a})
    
    figure(4)
    set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
    subplot(1,length(cond_data),a)
    plot_LFP_traces(LFPs, 1, LFP_timestamps)
    xlim(plot_lims)
    fixplot
    title(cond_names{a})
    
    figure(5)
    set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
    subplot(1,length(cond_data),a)
    plot_LFP_traces(LFPs, 1, LFP_timestamps, .05)
    xlim(plot_lims)
    fixplot
    title(cond_names{a})
    
    figure(6)
    set(gcf,'Units','normalized','Position',[.1 .2 .8 .4],'PaperPositionMode','auto')
    subplot(1,length(cond_data),a)
    [power_spectra, mean_power_spectrum] = LFP_power_spectra(LFPs,LFP_freq_bins,LFP_freq,100,0);
    mean_rows   = mean(mean_power_spectrum,2);
    mean_rows   = repmat(mean_rows,1,size(mean_power_spectrum,2));
    mean_power_spectrum = mean_power_spectrum ./ mean_rows;
    imagesc(mean_power_spectrum)
    axis xy
    fixplot
    xlabel('Time bin nr (x100ms)')
    ylabel('LFP frequency')
    
    title(['Normalised LFP power change - ' cond_names{a}])

end

figure(1)
subplot_equal_y
figure(2)
subplot_equal_clims
figure(4)
subplot_equal_y
figure(5)
subplot_equal_y




if save_figs
    print(1,fullfile(fig_save_dir,['Drive ' expt_name ' PSTHs']),'-dpng')
end

