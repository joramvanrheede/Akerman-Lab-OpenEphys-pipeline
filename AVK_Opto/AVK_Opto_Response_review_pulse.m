%%AVK_Opto_Response_review_pulse
clear all; 
close all;
str_Directory = ['D:\Multi_unit Coalated\2_Cortical\LED Powers\2020_01_22_A2']; %['D:\Multi_unit Coalated\1_POM\LED Powers\2019_11_26_A1'];
str_Experiment = ['2020_01_22-20-LED_powers'];%%['2019_11_28-11-LED_powers']%
str_Output = ['D:\Multi_unit Coalated\2_Cortical\Opto_response_Test\'];
load([str_Directory '\' str_Experiment '.mat']);


%% Some hardcoded variables; consider whether they should be inputs?

smoothing               = 1;    % smoothing on superimposed PSTHs
color_max_percentile    = .5;   % reduces colour range on spike density plots such that the highest colour_max_percentile doesn't influence the colour range
lowest_power            = 1;    % default intensity value for any conditions where LED signal was too low to register 

%% Check for inputs and set defaults

    channels        = ':';
    psth_bins       = [-0.05:0.001:0.2];
    artifact_win    = [-0.002 0.002];
    Layer_5_chans = 16:25;
plot_win        = [psth_bins(1) psth_bins(end)];
psth_bin_size   = mean(diff(psth_bins));

%% loads variable

LED_Powers = [ephys_data.conditions(:).LED_power];
Max_Power = max(LED_Powers);
iMax_Power = find(LED_Powers == Max_Power,1,'first');
Max_Power = Max_Power(1);
Max_Cond = ephys_data.conditions(iMax_Power);
opto_onset              = nanmedian([Max_Cond.LED_onset]);
early_resp_win = [artifact_win(2) artifact_win(2)+Max_Cond.LED_duration+0.005];
resp_win        = [artifact_win(2)+Max_Cond.LED_duration artifact_win(2)+Max_Cond.LED_duration+0.02];
    
% Pre-make figures
raster_fig          = figure;
set(gcf,'Units','normalized','Position',[.3 .1 .4 .9],'PaperPositionMode','auto')
spike_density_fig   = figure;
set(gcf,'Units','normalized','Position',[.3 .1 .4 .9],'PaperPositionMode','auto')
Beeswarm_fig   = figure;
set(gcf,'Units','normalized','Position',[.2 .1 .6 .9],'PaperPositionMode','auto')

%% gets spikes
spikes                      = [];
    
these_spikes                = Max_Cond.spikes(channels,:,:);
spikes(1:size(these_spikes,1),size(spikes,2)+(1:size(these_spikes,2)),1:size(these_spikes,3)) = these_spikes;

spikes(spikes == 0)     = NaN; % Default var gets padded with zeros; change padding zeros to NaN
    
    spont_spikes            = spikes-0.005;
    spikes                  = spikes - opto_onset;
    
    n_trials                = size(spikes,2);
    %% removes artifact
        q_artifact            	= spikes > artifact_win(1) & spikes < artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > artifact_win(1) & spont_spikes < artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
        q_artifact            	= spikes > Max_Cond.LED_duration + artifact_win(1) & spikes < Max_Cond.LED_duration +artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > Max_Cond.LED_duration + artifact_win(1) & spont_spikes < Max_Cond.LED_duration +artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
        
    %% plot Figures
figure(raster_fig)
    subplot(2,1,1);
    raster_plot(spikes,1);
    xlim([psth_bins(1) psth_bins(end)])
    fixplot
    set(gca,'FontSize',12)
    xlabel('Time (s)')
    ylabel('Chan no');
    title(['Raster plot for Max Power = ' num2str(Max_Power)])
    
    subplot(2,1,2);
    raster_plot(spikes,2);
    xlim([psth_bins(1) psth_bins(end)])
    fixplot
    set(gca,'FontSize',12)
    xlabel('Time (s)')
    ylabel('Trial No');
    title(['Raster plot for Max Power = ' num2str(Max_Power)])
    
% Spike density counts, channels x time
    figure(spike_density_fig)
    subplot(2,1,1)
    [~, spike_density_counts(:,:)]                    = spike_density_plot(spikes,1,psth_bins);
    set(gca,'FontSize',12)
    title(['Spike density for Max power = ' num2str(Max_Power)])
    xlabel('Time (s)');
    ylabel('Channel No');
    
    % PSTH_counts
    subplot(2,1,2);
    [~, psth_by_power(:)]                           	= psth(spikes, psth_bins);
    fixplot
    set(gca,'FontSize',12)
    
   
    
    %% quantify Response in Layer 5
        % Binned Early spike rate
    early_spike_rate_Max_Power(1:n_trials)                               = spike_rate_by_trial(spikes(Layer_5_chans,:,:), early_resp_win);
    
    % Spontaneous Early spike rate
    early_spont_rate_Max_Power(1:n_trials)                               = spike_rate_by_trial(spont_spikes(Layer_5_chans,:,:), early_resp_win);
    
        % Binned spike rate
    spike_rate_Max_Power(1:n_trials)                               = spike_rate_by_trial(spikes(Layer_5_chans,:,:), resp_win);
    
    % Spontaneous spike rate
    spont_rate_Max_Power(1:n_trials)                               = spike_rate_by_trial(spont_spikes(Layer_5_chans,:,:), resp_win);
    
    % PSTH_counts
    figure(Beeswarm_fig)
    subplot(2,2,1);
      [~, psth_by_power(:)]                           	= psth(spikes(Layer_5_chans,:,:), psth_bins);
       title('Opto Response L5');
     
    fixplot;
    set(gca,'FontSize',12)
    
      subplot(2,2,3);
    [~, psth_by_power(:)]                           	= psth(spont_spikes(Layer_5_chans,:,:), psth_bins);
      title('Spont Activity L5'); 
    fixplot;
    subplot_equal_y();
    set(gca,'FontSize',12)
    
    %% Some stats on the firing rate differences compared to control
    norm_data(1) =  kstest(early_spike_rate_Max_Power);
    norm_data(2) =  kstest(early_spont_rate_Max_Power);
    norm_data(3) =  kstest(spike_rate_Max_Power);
    norm_data(4) =  kstest(spont_rate_Max_Power);
    norm_data = ~(sum(norm_data)); %% KS test returns 1 if data not normally distributed so if any data points are 1 result will be 1. 
    % not() this data so that it returns true only if all data sets are
    % normal.
    
    if norm_data
        disp('Data Normally distributed, performing two sample t-test');
     [h_early,p_early] = ttest2(early_spike_rate_Max_Power,early_spont_rate_Max_Power,'tail','right'); 
     [h_rate,p_rate] = ttest2(spike_rate_Max_Power,spont_rate_Max_Power,'tail','right'); 
    else 
        disp('Data Not normally distributed performing Wilcoxon Ranksum test');
     [p_early,h_early] = ranksum(early_spike_rate_Max_Power,early_spont_rate_Max_Power,'tail','right'); 
    [p_rate,h_rate] = ranksum(spike_rate_Max_Power,spont_rate_Max_Power,'tail','right'); 
    end; 

% Bonferroni correction for multiple comparisons; this should be on the conservative side

%% Beeswarmplot peak spike rate
figure(Beeswarm_fig)
group_labels = [ones(1,numel(spike_rate_Max_Power)) 2*ones(1,numel(spont_rate_Max_Power))];
groups = [spike_rate_Max_Power spont_rate_Max_Power];
Labels = {'Opto';'Spont'};
subplot(2,2,2)
beeswarmplot(groups,group_labels,Labels);
fixplot
ylabel('Binned Response rate (Hz)')
title (['Difference significance p = ' num2str(p_rate)]);

group_labels = [ones(1,numel(early_spike_rate_Max_Power)) 2*ones(1,numel(early_spont_rate_Max_Power))];
groups = [early_spike_rate_Max_Power early_spont_rate_Max_Power];
Labels = {'Opto';'Spont'};
subplot(2,2,4)
beeswarmplot(groups,group_labels,Labels);
fixplot
title (['Difference significance p = ' num2str(p_early)]);
ylabel('Early Binned Response rate (Hz)')


%% save data
disp('Saving...');
disp(['Saving as...' str_Output str_Experiment '_Opto_Test.mat']);
save ([str_Output str_Experiment '_Opto_Test.mat'],'early_spike_rate_Max_Power','early_spont_rate_Max_Power','spike_rate_Max_Power','spont_rate_Max_Power','Max_Power','p_early','p_rate','h_early','h_rate');

if ~exist([str_Output str_Experiment], 'dir')
       mkdir([str_Output str_Experiment])
end
    
disp(['Saving figures to...' str_Output str_Experiment '\']);
disp('Beeswarm..');
saveas(Beeswarm_fig,[str_Output str_Experiment '\' str_Experiment '_Beeswarm.fig']);
saveas(Beeswarm_fig,[str_Output str_Experiment '\' str_Experiment '_Beeswarm.png']);
disp('Raster...');
saveas(raster_fig,[str_Output str_Experiment '\' str_Experiment '_raster.fig']);
saveas(raster_fig,[str_Output str_Experiment '\' str_Experiment '_raster.png']);
disp('Spike_Density');
saveas(spike_density_fig,[str_Output str_Experiment '\' str_Experiment '_spike_density.fig']);
saveas(spike_density_fig,[str_Output str_Experiment '\' str_Experiment '_spike_density.png']);
disp('Saved');
%clear all;
%close all;
