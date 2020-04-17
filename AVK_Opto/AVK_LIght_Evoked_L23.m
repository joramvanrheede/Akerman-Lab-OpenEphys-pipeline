%%AVK_Opto_Response_review_pulse
clear all; 
close all;
str_Directory = ['D:\Multi_unit Coalated\1_POM\LED Powers\']; %['D:\Multi_unit Coalated\1_POM\LED Powers\2019_11_26_A1'];
str_L5_Directory = ['D:\Multi_unit Coalated\1_POM\Layer_5_pulse_response'];
str_Experiment = ['2019_11_26_A1'];%%['2019_11_28-11-LED_powers']%
str_Output = ['D:\Multi_unit Coalated\2_Cortical\L2_3_evoked\'];


Target_Response = 90; %25 150 
L2_3_Evoked = struct
L2_3_Evoked.Target_Response = Target_Response;
%% Some hardcoded variables; consider whether they should be inputs?

smoothing               = 1;    % smoothing on superimposed PSTHs
color_max_percentile    = .5;   % reduces colour range on spike density plots such that the highest colour_max_percentile doesn't influence the colour range
lowest_power            = 1;    % default intensity value for any conditions where LED signal was too low to register 

%% Check for inputs and set defaults

    channels        = ':';
    psth_bins       = [-0.05:0.001:0.2];
    artifact_win    = [-0.002 0.002];
    Layer_23_chans = 1:9;
plot_win        = [psth_bins(1) psth_bins(end)];
psth_bin_size   = mean(diff(psth_bins));

%% load layer 5 and find target Power
load([str_L5_Directory '\' str_Experiment '_L5_Evoked.mat']);
L2_3_Evoked.Layer5_Response = File_results.Mean_Light_Delta;

%% loads variable
matfiles = dir(fullfile(str_Directory,str_Experiment, '*.mat'));
nfiles = length(matfiles);
for i = 1 : nfiles
    disp(['Loading...' matfiles(i).name]);
    load(fullfile(matfiles(i).folder, matfiles(i).name));
    cond_data = ephys_data.conditions
    Layer_5_resp = Layer_response(cond_data,16:25);
    [a(i) ,power_index(i)] = min(abs(Layer_5_resp-Target_Response))
   
end;
   if numel(power_index) >1 
       [~, a] = min(a);
       power_index = power_index(a);
       load(fullfile(matfiles(a).folder, matfiles(a).name));
   end;
cond_data = ephys_data.conditions
    Layer_5_resp = Layer_response(cond_data,16:25);
  
    opto_onset              = nanmedian([cond_data(:).LED_onset]);
    opto_duration           = nanmedian([cond_data(:).LED_duration]);         
       opto_powers     = [cond_data(:).LED_power]';
    opto_powers(isnan(opto_powers)) = lowest_power;
    
    % if any powers were NaN and replaced with a minimum nr, powers and conditions struct
    % will not be in ascending order; fix:
    [opto_powers]           = unique(opto_powers);
    cond_data_powers      	= [cond_data.LED_power];
    cond_data_powers(isnan(cond_data_powers)) = lowest_power;
  
    L2_3_Evoked.Layer5_Response_Calculated = Layer_5_resp(power_index);
    
    this_power                  = opto_powers(power_index);
    q_power                     = cond_data_powers == this_power;
    
%% gets spikes
    
    power_inds                  = find(q_power);
    Target_Power = cond_data_powers(power_inds);
    spikes                      = [];
        Target_Cond = cond_data(power_inds)
        these_spikes                = Target_Cond.spikes(Layer_23_chans,:,:);

        spikes(1:size(these_spikes,1),size(spikes,2)+(1:size(these_spikes,2)),1:size(these_spikes,3)) = these_spikes;
        
    spikes(spikes == 0)     = NaN; % Default var gets padded with zeros; change padding zeros to NaN
    
    spont_spikes            = spikes - 0.005;
    spikes                  = spikes - opto_onset;
    
    n_trials                = size(spikes,2);
  
early_resp_win = [artifact_win(2) artifact_win(2)+Target_Cond.LED_duration+0.005];
resp_win        = [artifact_win(2)+Target_Cond.LED_duration artifact_win(2)+Target_Cond.LED_duration+0.02];% 20ms following offset
mid_resp        = [Target_Cond.LED_duration+0.02 Target_Cond.LED_duration+0.05];% 20 -50ms following offset
late_resp        = [Target_Cond.LED_duration+0.05 Target_Cond.LED_duration+0.1];% 50 -100ms following offset

    
% Pre-make figures
raster_fig          = figure;
set(gcf,'Units','normalized','Position',[.3 .1 .4 .9],'PaperPositionMode','auto')
spike_density_fig   = figure;
set(gcf,'Units','normalized','Position',[.3 .1 .4 .9],'PaperPositionMode','auto')
Beeswarm_fig   = figure;
set(gcf,'Units','normalized','Position',[.2 .1 .6 .9],'PaperPositionMode','auto')

    %% removes artifact
        q_artifact            	= spikes > artifact_win(1) & spikes < artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > artifact_win(1) & spont_spikes < artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
        q_artifact            	= spikes > Target_Cond.LED_duration + artifact_win(1) & spikes < Target_Cond.LED_duration +artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > Target_Cond.LED_duration + artifact_win(1) & spont_spikes < Target_Cond.LED_duration +artifact_win(2);
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
    title(['Raster plot for Max Power = ' num2str(Target_Power)])
    
    subplot(2,1,2);
    raster_plot(spikes,2);
    xlim([psth_bins(1) psth_bins(end)])
    fixplot
    set(gca,'FontSize',12)
    xlabel('Time (s)')
    ylabel('Trial No');
    title(['Raster plot for Max Power = ' num2str(Target_Power)])
    
% Spike density counts, channels x time
    figure(spike_density_fig)
    
    % PSTH_counts
    subplot(3,1,2);
    [~, psth_L2_3(:)]                           	= psth(spikes, psth_bins);
    fixplot
    set(gca,'FontSize',12)
    L2_3_Evoked.psth_L2_3 = psth_L2_3;
    ylim([0 50]);
    
    % PSTH_counts
      subplot(3,1,3);
    [~, psth_Spont(:)]                           	= psth(spont_spikes(:,:,:), psth_bins);
     title('Spont Activity L2_3'); 
    fixplot;
    subplot_equal_y();
    set(gca,'FontSize',12)
    L2_3_Evoked.psth_Spont = psth_Spont;
    ylim([0 50]);
    
   subplot(3,1,1)
    [~, spike_density_counts(:,:)]                    = spike_density_plot(spikes,1,psth_bins);
    set(gca,'FontSize',12)
    title(['Spike density for Max power = ' num2str(Target_Power)])
    xlabel('Time (s)');
    ylabel('Channel No');
    
   
    
    %% quantify Response in Layer 2_3
      %   Binned Early spike rate
    early_spike_rate(1:n_trials)                               = spike_rate_by_trial(spikes(:,:,:), early_resp_win);
    early_spont_rate(1:n_trials)                               = spike_rate_by_trial(spont_spikes(:,:,:), early_resp_win);
    L2_3_Evoked.early_spike_rate = early_spike_rate;
    L2_3_Evoked.early_spont_rate = early_spont_rate;
    
    
        % Binned spike rate
    spike_rate(1:n_trials)                               = spike_rate_by_trial(spikes(:,:,:), resp_win);
    spont_rate(1:n_trials)                               = spike_rate_by_trial(spont_spikes(:,:,:), resp_win);
    L2_3_Evoked.spike_rate = spike_rate;
    L2_3_Evoked.spont_rate = spont_rate;
    
     % Binned Mid Spike Rate
    mid_spike_rate(1:n_trials)                               = spike_rate_by_trial(spikes(:,:,:), mid_resp);
    mid_spont_rate(1:n_trials)                               = spike_rate_by_trial(spont_spikes(:,:,:), mid_resp);
     L2_3_Evoked.mid_spike_rate = mid_spike_rate;
    L2_3_Evoked.mid_spont_rate = mid_spont_rate;
   
    
    %Binned Late Spike Rate
    late_spike_rate(1:n_trials)                               = spike_rate_by_trial(spikes(:,:,:), late_resp);
    late_spont_rate(1:n_trials)                               = spike_rate_by_trial(spont_spikes(:,:,:), late_resp);
    L2_3_Evoked.late_spike_rate = late_spike_rate;
    L2_3_Evoked.late_spont_rate = late_spont_rate;
   
    
    %% Some stats on the firing rate differences compared to control
    norm_data(1) =  kstest(early_spike_rate);
    norm_data(2) =  kstest(early_spont_rate);
    norm_data(3) =  kstest(spike_rate);
    norm_data(4) =  kstest(spont_rate);
    norm_data(5) =  kstest(mid_spike_rate);
    norm_data(6) =  kstest(mid_spont_rate);
    norm_data(7) =  kstest(late_spike_rate);
    norm_data(8) =  kstest(late_spike_rate);    
    norm_data = ~(sum(norm_data)); %% KS test returns 1 if data not normally distributed so if any data points are 1 result will be 1. 
    % not() this data so that it returns true only if all data sets are
    % normal.
    
    if norm_data
        disp('Data Normally distributed, performing two sample t-test');
     [h_early,p_early] = ttest2(early_spike_rate,early_spont_rate); 
     [h_rate,p_rate] = ttest2(spike_rate,spont_rate);
     [h_mid,p_mid] = ttest2(mid_spike_rate,mid_spont_rate);
     [h_late,p_late] = ttest2(late_spike_rate,late_spont_rate);
     
    else 
        disp('Data Not normally distributed performing Wilcoxon Ranksum test');
     [p_early,h_early] = ranksum(early_spike_rate,early_spont_rate); 
     [p_rate,h_rate] = ranksum(spike_rate,spont_rate); 
     [p_mid,h_mid] = ranksum(mid_spike_rate,mid_spont_rate);
     [p_late,h_late] = ranksum(late_spike_rate,late_spont_rate);
     
    end; 

    L2_3_Evoked.p_early = p_early;
    L2_3_Evoked.p_rate = p_early;
    L2_3_Evoked.p_mid = p_mid;
    L2_3_Evoked.p_late = p_late;
    
    L2_3_Evoked.h_early = h_early;
    L2_3_Evoked.h_rate = h_early;
    L2_3_Evoked.h_mid = h_mid;
    L2_3_Evoked.h_late = h_late;
    

%% Beeswarmplot peak spike rate
figure(Beeswarm_fig)
group_labels = [ones(1,numel(spike_rate)) 2*ones(1,numel(spont_rate))];
groups = [spike_rate spont_rate];
Labels = {'Opto';'Spont'};
subplot(2,2,1)
beeswarmplot(groups,group_labels,Labels);
fixplot
ylabel('Binned Response rate (Hz)')
title (['Difference significance p = ' num2str(p_rate)]);

group_labels = [ones(1,numel(early_spike_rate)) 2*ones(1,numel(early_spont_rate))];
groups = [early_spike_rate early_spont_rate];
Labels = {'Opto';'Spont'};
subplot(2,2,2)
beeswarmplot(groups,group_labels,Labels);
fixplot
title (['Difference significance p = ' num2str(p_early)]);
ylabel('Early Binned Response rate (Hz)')

group_labels = [ones(1,numel(mid_spike_rate)) 2*ones(1,numel(mid_spont_rate))];
groups = [mid_spike_rate mid_spont_rate];
Labels = {'Opto';'Spont'};
subplot(2,2,3)
beeswarmplot(groups,group_labels,Labels);
fixplot
title (['Difference significance p = ' num2str(p_mid)]);
ylabel('mid Binned Response rate (Hz)')

group_labels = [ones(1,numel(late_spike_rate)) 2*ones(1,numel(late_spont_rate))];
groups = [late_spike_rate late_spont_rate];
Labels = {'Opto';'Spont'};
subplot(2,2,4)
beeswarmplot(groups,group_labels,Labels);
fixplot
title (['Difference significance p = ' num2str(p_late)]);
ylabel('late Binned Response rate (Hz)')


%% save data
disp('Saving...');
disp(['Saving as...' str_Output str_Experiment '_' num2str(Target_Response) '_Opto_Test.mat']);
save ([str_Output str_Experiment '_' num2str(Target_Response) '_L23_Resp.mat'],'L2_3_Evoked');

if ~exist([str_Output str_Experiment], 'dir')
       mkdir([str_Output str_Experiment])
end
    
disp(['Saving figures to...' str_Output str_Experiment '\']);
disp('Beeswarm..');
saveas(Beeswarm_fig,[str_Output str_Experiment '\' str_Experiment '_' num2str(Target_Response) '_Beeswarm.fig']);
saveas(Beeswarm_fig,[str_Output str_Experiment '\' str_Experiment '_' num2str(Target_Response) '_Beeswarm.png']);
disp('Raster...');
saveas(raster_fig,[str_Output str_Experiment '\' str_Experiment '_' num2str(Target_Response) '_raster.fig']);
saveas(raster_fig,[str_Output str_Experiment '\' str_Experiment '_' num2str(Target_Response) '_raster.png']);
disp('Spike_Density');
saveas(spike_density_fig,[str_Output str_Experiment '\' str_Experiment '_' num2str(Target_Response) '_spike_density.fig']);
saveas(spike_density_fig,[str_Output str_Experiment '\' str_Experiment '_' num2str(Target_Response) '_spike_density.png']);
disp('Saved');
%clear all;
%close all;
    %}