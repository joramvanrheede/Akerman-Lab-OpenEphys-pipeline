% single experiment plots for frequency experiment

close all

experiment          = sdata(11).expt(1); % which experiment

%% User parameters
summarise_channels  = [1:16]; % include these channels

% for raster plots:
trialrange          = [1 10]; % [min max] - don't exceed max nr of trials; else errors result.
x_ax_lims           = [0 12]; % limits for x-axes
condition_name    	= 'Frequency';
condition_units   	= 'Hz';

% figure saving:
save_figs        	= true; % if false, no figures will be saved
save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output'; % figures will be saved in this folder
figure_dpi          = 300; % dots per inch for saved figures (Standard = 150, HQ = 300, XHQ = 600)

respwinsize         = .04;

%% Fixed for experiment type 'Frequency'
split_conditions    = [1 4 6];  % split by these conditions, summarise over others
split_plots         = [4 6];    % [4 6] works

% for rasterplots:
split_raster_plots  = [4 1];    % split plots by these conditions
split_figures       = [6];      % split 

%% End of user input, code execution starts here

% Work out how to separate conditions

% retrieve condition matrix
condition_mat       = experiment.condition_mat;

% separate out the column with frequencies (relevant for this experiment type)
frequencies         = condition_mat(:,4);

% find unique values for separating axes
split_cond_mat      = condition_mat(:,split_conditions);
[split_cond_rows, indxa, cond_inds] = unique(split_cond_mat,'rows');

% find unique values for separating plot lines
split_plot_mat      = condition_mat(:,split_plots);
[split_plot_rows, indxa, cond_plot_inds] = unique(split_plot_mat,'rows');

[a,b,cond_plot_inds] = unique(cond_plot_inds);

%% Make PSTH figures
figure(1)
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .6])
figure(2)
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .6])
figure(3)
set(gcf,'Units','normalized')
set(gcf,'Position',[.3 0 .4 1])
profile_plots   = [];
for a = 1:size(split_cond_rows,1)
    
    these_conds         = split_cond_rows(a,:);
    sum_inds            = cond_inds == a;
    
    this_whisk_psth  	= mean(experiment.whisk_win_rel(sum_inds,:,:),1);
    
    figure(1)
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(experiment.whiskwinedges(1:end-1),smooth(mean(squeeze(this_whisk_psth(:,summarise_channels,:)),1),5),'LineWidth',2)
    hold on
    this_frequency      = mean(frequencies(sum_inds));
    
    stim_times          = [0:1/round(this_frequency):5.5];
    
    plot(stim_times,zeros(size(stim_times))+0.2,'k^','MarkerSize',5,'LineWidth',2)
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    
    xlim([-1 4.5])

    profile_plots(a).stim_times     = stim_times;
    for b = 1:length(stim_times)
        
        qstimwin            = experiment.whiskwinedges > stim_times(b) & experiment.whiskwinedges < (stim_times(b) + respwinsize);
        profile_segment     = squeeze(mean(mean(experiment.whisk_profile(sum_inds,summarise_channels,qstimwin),1),2));
        
        if isempty(profile_segment)
            profile_segment = NaN;
        end
        
        profile_plots(a).resp_peaks(b)  = max(profile_segment);
    end
    
    axis tight
    ylims = ylim;
    ylim([0 ylims(2)*1.1])
    
    figure(2)
    xlim([0 4])
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(profile_plots(a).stim_times,profile_plots(a).resp_peaks,'LineWidth',2,'MarkerSize',25);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
    axis tight
    ylims = ylim;
    ylim([0 ylims(2)*1.1])
    
    figure(3)
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    xlim([0 6])
    profile_plots(a).resp_peaks = profile_plots(a).resp_peaks(~isnan(profile_plots(a).resp_peaks));
    
    % hacked this to get response at ~3 seconds even if stim_times has a
    % value not quite equal to 3...
    plot(1:5,[profile_plots(a).resp_peaks([1:4]) profile_plots(a).resp_peaks(find(round(stim_times*100) == 300))],'LineWidth',2,'MarkerSize',25);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
    axis tight
    ylims = ylim;
    ylim([0 ylims(2)*1.1])
end

%% Aesthetic stuff for figure 1 + save option

% select figure 1
figure(1) 

% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);

% Background colour
set(gcf,'Color',[1 1 1]) 
% Title for column 1
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
% Title for column 2
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')

% save figure
figure(1)
if save_figs
    fileseps = regexp(experiment.filename,'/.');
    experiment_plot_folder      = [save_folder filesep experiment.filename((fileseps(end)+1):end-4)];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
    print(gcf,[experiment_plot_folder filesep ' PSTH plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
end

%% Aesthetic stuff for figure 2 + save option

% select figure
figure(2) 

% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);

% Background colour
set(gcf,'Color',[1 1 1])
% Title for column 1
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
% Title for column 2
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')

%% Aesthetic stuff for figure 3 + save option

% select figure
figure(3) 

% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);


%% rasterplots

% use channels_to_rasterplots figure to generate 
fighandles = channels_to_rasterplots(experiment.filename,summarise_channels,split_raster_plots,split_figures,trialrange,x_ax_lims,condition_name, condition_units);

% save figure
if save_figs
    for a = 1:length(fighandles)
        figure(fighandles(a))
        experiment_plot_folder      = [save_folder filesep experiment.filename(1:end-4)];
        if ~isdir(experiment_plot_folder)
            mkdir(experiment_plot_folder)
        end
        print(gcf,[experiment_plot_folder filesep 'Rasterplot stimulator ' num2str(a) ' chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    end
end

%% Look at peak instantaneous firing rate between LED and no LED

contrast_conds      = [1 6];
contrast_cond_mat   = condition_mat(:,contrast_conds);
[contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
contrast_rates      = [];
contrast_times      = [];
for a  = unique(cond_inds)'
    
    contrast_rates  = [contrast_rates; mean(mean(experiment.whisk_peak_rate(cond_inds == a,summarise_channels)))];
    contrast_times  = [contrast_times; mean(mean(experiment.whisk_peak_time(cond_inds == a,summarise_channels)))];
end

if contrast_rates(3) > contrast_rates(4)
    P_whisk_stim = 1;
    A_whisk_stim = 2;
else
    P_whisk_stim = 2;
    A_whisk_stim = 1;
end

P_rate_LED_on       = contrast_rates(P_whisk_stim);     % 
P_rate_LED_off      = contrast_rates(P_whisk_stim+2);   % 
A_rate_LED_on       = contrast_rates(A_whisk_stim);     % 
A_rate_LED_off      = contrast_rates(A_whisk_stim+2);   % 

PA_ratio_LED_on     = P_rate_LED_on / A_rate_LED_on;
PA_ratio_LED_off    = P_rate_LED_off / A_rate_LED_off;

P_LED_onoff_ratio   = P_rate_LED_on / P_rate_LED_off;
A_LED_onoff_ratio 	= A_rate_LED_on / A_rate_LED_off;

% 

figure
bar_handle  = bar([P_rate_LED_on P_rate_LED_off; A_rate_LED_on A_rate_LED_off],'LineWidth',2);
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','bold')
ylim([0 1000])
title('Response size, P vs. A whisker, LED ON vs. OFF')
ylabel('Peak stimulus-evoked firing rate (Hz)')
xlabel('Whisker (1 = principal, 2 = adjacent)')
legend({'LED ON' 'LED OFF'})

if save_figs
    print(gcf,[experiment_plot_folder filesep ' Bar graph P vs A - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
end

