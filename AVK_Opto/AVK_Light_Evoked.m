%%AVK_Opto_Response_review_pulse
clear all; 
close all;
str_Directory = ['D:\Multi_unit Coalated\2_Cortical\Laser_pulse']; %['D:\Multi_unit Coalated\1_POM\LED Powers\2019_11_26_A1'];
str_Experiment = ['2019_05_21'];
str_Output = ['D:\Multi_unit Coalated\2_Cortical\Layer_5_pulse_response\'];

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

matfiles = dir(fullfile(str_Directory,str_Experiment, '*.mat'));
nfiles = length(matfiles);
File_results = struct();
for i =1: nfiles
    disp(['Loading...' matfiles(i).name]);
    load(fullfile(matfiles(i).folder, matfiles(i).name));
    
    cond_data       = ephys_data.conditions;
    opto_powers     = [cond_data(:).LED_power]';
    opto_powers(isnan(opto_powers)) = lowest_power;
     
    % if any powers were NaN and replaced with a minimum nr, powers and conditions struct
    % will not be in ascending order; fix:
    [opto_powers]           = unique(opto_powers);
    cond_data_powers      	= [cond_data.LED_power];
    cond_data_powers(isnan(cond_data_powers)) = lowest_power;

    
    no_Powers = numel(opto_powers);
    opto_onset              = nanmedian([cond_data(:).LED_onset]);
    opto_duration           = nanmedian([cond_data(:).LED_duration]);         
    early_resp_win = [artifact_win(2) artifact_win(2)+opto_duration+0.005];
    resp_win        = [artifact_win(2)+opto_duration artifact_win(2)+opto_duration+0.02];

    File_results(i).Light_evoked = NaN*ones(no_Powers,200);
    File_results(i).Spont = NaN*ones(no_Powers,200);
    File_results(i).Light_Delta =NaN*ones(no_Powers,200);

    for k = 1 : no_Powers
    this_power                  = opto_powers(k);
    File_results(i).Powers(k) = this_power;
    
    q_power                     = cond_data_powers == this_power;
    
    
    power_inds                  = find(q_power);
    spikes                      = [];
    for c = 1:length(power_inds)
        
        these_spikes                = cond_data(power_inds(c)).spikes(channels,:,:);

        spikes(1:size(these_spikes,1),size(spikes,2)+(1:size(these_spikes,2)),1:size(these_spikes,3)) = these_spikes;
        
    end
    spikes(spikes == 0)     = NaN; % Default var gets padded with zeros; change padding zeros to NaN
    
    spont_spikes            = spikes;
    spikes                  = spikes - opto_onset;
    
    n_trials                = size(spikes,2);
       
        %% removes artifact
        q_artifact            	= spikes > artifact_win(1) & spikes < artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > artifact_win(1) & spont_spikes < artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
        q_artifact            	= spikes > opto_duration + artifact_win(1) & spikes < opto_duration +artifact_win(2);
        spikes(q_artifact)    	= NaN;
        
        q_artifact            	= spont_spikes > opto_duration + artifact_win(1) & spont_spikes < opto_duration +artifact_win(2);
        spont_spikes(q_artifact) =NaN;
        
            % Binned spike rate
            spike_rate(1:n_trials)                               = spike_rate_by_trial(spikes(Layer_5_chans,:,:), resp_win);
    
            % Spontaneous spike rate
            spont_rate(1:n_trials)                               = spike_rate_by_trial(spont_spikes(Layer_5_chans,:,:), resp_win);
     
      File_results(i).Light_evoked(k,1:numel(spike_rate)) = spike_rate;
       File_results(i).Spont(k,1:numel(spont_rate)) = spont_rate;
      File_results(i).Light_Delta(k,1:numel(spont_rate)) = spike_rate-spont_rate;
      
    File_results(i).Median_Light_evoked(k) = nanmedian(spike_rate);
    File_results(i).Mean_Light_evoked(k) = nanmean(spont_rate);
    
     File_results(i).Median_Light_Delta(k) = nanmedian(File_results(i).Light_Delta(k,:));
    File_results(i).Mean_Light_Delta(k) = nanmean(File_results(i).Light_Delta(k,:));
   
      %% Some stats on the firing rate differences compared to control
    norm_data(1) =  kstest(spike_rate);
    norm_data(2) =  kstest(spont_rate);
    norm_data = ~(sum(norm_data)); %% KS test returns 1 if data not normally distributed so if any data points are 1 result will be 1. 
    % not() this data so that it returns true only if all data sets are
    % normal.
    
    if norm_data
        disp('Data Normally distributed, performing two sample t-test');
     [h_rate,p_rate] = ttest2(spike_rate,spont_rate,'tail','right'); 
    else 
        disp('Data Not normally distributed performing Wilcoxon Ranksum test');
     [p_rate,h_rate] = ranksum(spike_rate,spont_rate,'tail','right'); 
    end; 
     File_results(i).Diff_From_Spont(k) = h_rate;
    File_results(i).Diff_From_Spont_p(k) = p_rate;  
    end;
 end;

%%Plotting
All_Powers = [];
All_Light_Delta =[];
All_Light_err = [];
All_Diff = [];
for k = 1: numel(File_results)
    All_Powers = [All_Powers File_results(k).Powers];
    All_Light_Delta = [All_Light_Delta File_results(k).Mean_Light_Delta];
    All_Light_err = [All_Light_err;(nanstd(File_results(k).Light_Delta,0,2))./sqrt(size(File_results(k).Light_Delta,2))];
    All_Diff = [All_Diff File_results(k).Diff_From_Spont];
end;

file_x = 1 : numel(All_Powers);
Star_index = NaN*ones(numel(All_Powers),1);

Star_index(find(All_Diff==1)) = max(All_Light_Delta)+max(All_Light_err)+5;
Star_index(Star_index ==0) =NaN;

ha = figure();
set(gcf,'Units','normalized','Position',[.1 .3 .8 .4],'PaperPositionMode','auto')
hold on;
if numel(unique(All_Powers)) >1
errorbar(All_Powers,All_Light_Delta,All_Light_err','k','LineWidth',2,'MarkerSize',30)
plot(All_Powers,Star_index,'k*','LineWidth',2,'MarkerSize',10);
title('Opto power vs Layer 5 Spike Rate Increase')
xlabel('Opto power')
else
errorbar(file_x,All_Light_Delta,All_Light_err,'k','LineWidth',2,'MarkerSize',30)
plot(file_x,Star_index,'k*','LineWidth',2,'MarkerSize',10);
title('Laser Power vs Layer 5 Spike Rate Increase')
xlabel('Recording number')
end;
    
hline = refline([0 0]);
hline.Color = 'r';
ylabel('Delta Spike Rate (hz)')
fixplot

%% Save_data
disp(['Saving_data']);
disp(['Saving as...' str_Output str_Experiment '_L5_Evoked.mat']);
save ([str_Output str_Experiment '_L5_Evoked.mat'],'File_results');
disp(['Saved']);

if ~exist([str_Output str_Experiment], 'dir')
       mkdir([str_Output str_Experiment])
end
disp(['Saving Figures']);
saveas(ha,[str_Output str_Experiment '\' str_Experiment '_Response.fig']);
saveas(ha,[str_Output str_Experiment '\' str_Experiment '_Response.png']);
disp(['Saved']);



    
