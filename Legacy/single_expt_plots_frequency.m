% single experiment plots for frequency experiment

close all

experiment          = sdata(1).expt(1); % which experiment

%% User parameters
summarise_channels  = [13:32]; % include these channels

% for raster plots:
trialrange          = [1 10]; % [min max] - don't exceed max nr of trials; else errors result.
x_ax_lims           = [0 12]; % limits for x-axes
condition_name    	= 'Frequency';
condition_units   	= 'Hz';

% figure saving:
save_figs        	= false; % if false, no figures will be saved
save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output'; % figures will be saved in this folder
figure_dpi          = 300; % dots per inch for saved figures (Standard = 150, HQ = 300, XHQ = 600)

respwinsize         = .04; % post-stimulus window for assessing response to successive stimuli
prewinsize        	= .03; % pre-stimulus window for assessing response size minus spont

qchannelim          = false;
qpsth               = true;
qtrackingplots      = false;
qraster             = false;

first_vs_nr         = 4;
timed_vs_nr         = 4;

%% Fixed for experiment type 'Frequency'
split_conditions    = [1 4 6];  % split by these conditions, summarise over others
split_plots         = [4 6];    % [4 6] works

% for rasterplots:
split_raster_plots  = [4 1];    % split plots by these conditions
split_figures       = [6];      % split

%%
if save_figs
    fileseps = regexp(experiment.filename,'/.');
    experiment_plot_folder      = [save_folder filesep experiment.filename((fileseps(end)+1):end-4)];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
end

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
    profile_plots(a).frequency      = these_conds(2);
    profile_plots(a).LEDtime        = these_conds(1);
    profile_plots(a).stimulator     = these_conds(3);
    for b = 1:length(stim_times)
        
        qstimwin            = experiment.whiskwinedges > stim_times(b) & experiment.whiskwinedges < (stim_times(b) + respwinsize);
        qprestimwin         = experiment.whiskwinedges > (stim_times(b) - prewinsize) & experiment.whiskwinedges < (stim_times(b));
        
        profile_segment     = squeeze(mean(mean(experiment.whisk_profile(sum_inds,summarise_channels,:),1),2));
        pre_profile_segment = squeeze(mean(mean(experiment.whisk_profile(sum_inds,summarise_channels,:),1),2));
        
        if isempty(profile_segment)
            profile_segment = NaN;
        end
        
        if isempty(pre_profile_segment)
            pre_profile_segment = NaN;
        end
        
        profile_plots(a).resp_peaks(b)      = max(profile_segment);
        profile_plots(a).pre_resps(b)       = mean(pre_profile_segment);
        profile_plots(a).adjusted_resps(b)  = profile_plots(a).resp_peaks(b) - profile_plots(a).pre_resps(b);
    end
    
    axis tight
    ylims = ylim;
    ylim([0 ylims(2)*1.1])
    
    figure(2)
    xlim([0 4])
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    plot(profile_plots(a).stim_times,profile_plots(a).adjusted_resps,'LineWidth',2,'MarkerSize',25);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
    axis tight
    ylims = ylim;
    ylim([0 ylims(2)*1.1])
    
    figure(3)
    subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),cond_plot_inds(a))
    xlim([0 6])
    profile_plots(a).adjusted_resps = profile_plots(a).adjusted_resps(~isnan(profile_plots(a).adjusted_resps));
    
    
    % find stimulus times at t=0 and at at approsimately t=1s, 2s, etc...
    timed_inds  = [1 find(round(stim_times*100) == 100) find(round(stim_times*100) == 200) find(round(stim_times*100) == 300)];
    
    plot(1:4,[profile_plots(a).adjusted_resps(timed_inds)],'LineWidth',2,'MarkerSize',25);
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',16)
    hold on
    axis tight
    ylims = ylim;
    if ylims(2) == 0
        ylims(2) = .1;
    end
    ylim([0 ylims(2)*1.1])
    
    profile_plots(a).first_resps    = profile_plots(a).adjusted_resps(1:4);
    profile_plots(a).timed_resps    = profile_plots(a).adjusted_resps(timed_inds);
    
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

% Title for column 1
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),1)
title(gca,'Stimulator 1')
% Title for column 2
subplot(length(unique(condition_mat(:,split_plots(1)))),length(unique(condition_mat(:,split_plots(2)))),2)
title(gca,'Stimulator 2')

% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);

% save figure
if save_figs
    print(1,[experiment_plot_folder filesep ' PSTH plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    print(2,[experiment_plot_folder filesep ' Peakresp plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
end


%%
first_pks   = cell2mat({profile_plots.first_resps}');
timed_pks   = cell2mat({profile_plots.timed_resps}');

first_ratio = first_pks(:,first_vs_nr)./first_pks(:,1);
timed_ratio = first_pks(:,timed_vs_nr)./timed_pks(:,1);
frequencies = [profile_plots.frequency]';
LEDtime     = [profile_plots.LEDtime]';
LEDtime     = LEDtime < 2; % if LED is on early, we're in the LED on condition, so 1 = on, 0 = off
stimulator  = [profile_plots.stimulator]';

% if principal whisker is on stimulator 2, invert stimulator code so that 1
% is  principal and 2 is adjacent
if mean(mean(experiment.whisk_peak_rate(stimulator == 1,:))) < mean(mean(experiment.whisk_peak_rate(stimulator == 2,:)))
    stimulator = 3 - stimulator;
end

uniqfreqs   = unique(frequencies);
uniqLEDs    = unique(LEDtime);
uniqstims   = unique(stimulator);
nLEDconds   = length(uniqLEDs);
nstimulators = length(uniqstims);

figure
for a = 1:nLEDconds
    for b = 1:nstimulators
        qLED    = LEDtime == uniqLEDs(a);
        qstim   = stimulator == uniqstims(b);
        
        % set up variables for loop over uniqfreqs
        first_resp_plotmeans    = [];
        first_resp_plotserrs    = [];
        timed_resp_plotmeans    = [];
        timed_resp_plotserrs    = [];
        plotfreqs               = [];
        
        % loop over uniqfreqs
        for c = 1:length(uniqfreqs)
            qfreq   = frequencies == uniqfreqs(c);
            qall    = qfreq & qLED & qstim;
            
            % if no data meet these criteria, go to next iteration
            if sum(qall) == 0
                continue
            end
            
            % get ratios that meet criteria
            these_first_ratios      = first_ratio(qall);
            these_timed_ratios      = timed_ratio(qall);
            
            % get means and standard errors for plotting
            first_resp_plotmeans    = [first_resp_plotmeans; mean(these_first_ratios)];
            first_resp_plotserrs    = [first_resp_plotserrs; serr(these_first_ratios)];
            timed_resp_plotmeans    = [timed_resp_plotmeans; mean(these_timed_ratios)];
            timed_resp_plotserrs    = [timed_resp_plotserrs; serr(these_timed_ratios)];
            plotfreqs               = [plotfreqs, uniqfreqs(c)];
        end
        
        figure(21)
        subplot(1,2,b)
        if ~uniqLEDs(a)
            errorbar(plotfreqs,first_resp_plotmeans,first_resp_plotserrs,'k-','LineWidth',2)
        else
            errorbar(plotfreqs,first_resp_plotmeans,first_resp_plotserrs,'r-','LineWidth',2)
        end
        xlim([0 16])
        set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
        if b == 1
            title('Principal whisker, first vs n-th response')
        elseif b == 2
            title('Adjacent whisker, first vs n-th response')
        end
        xlabel('Stimulus trigger frequency')
        ylabel('Response ratio (relative to 1st)')
        hold on
        
        figure(22)
        subplot(1,2,b)
        if ~uniqLEDs(a)
            errorbar(plotfreqs,timed_resp_plotmeans,timed_resp_plotserrs,'k-','LineWidth',2)
        else
            errorbar(plotfreqs,timed_resp_plotmeans,timed_resp_plotserrs,'r-','LineWidth',2)
        end
        xlim([0 16])
        set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
        
        if b == 1
            title('Principal Whisker - First vs n seconds response')
        elseif b ==2 
            title('Adjacent Whisker - First vs n seconds response')
        end
        
        xlabel('Stimulus trigger frequency')
        ylabel('Response ratio (relative to 1st)')
        hold on
    end
end
figure(21)
set(gcf,'Units','Normalized')
set(gcf,'Position',[.2 .4 .6 .4])
% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);
legend({'LED OFF' 'LED ON'})

figure(22)
set(gcf,'Units','Normalized')
set(gcf,'Position',[.2 .4 .6 .4])
% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);
legend({'LED OFF' 'LED ON'})

if save_figs 
    print(21,[experiment_plot_folder filesep 'First Tracking plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
    print(22,[experiment_plot_folder filesep 'Timed Tracking plot - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
end



%% rasterplots

if qraster
    % use channels_to_rasterplots figure to generate
    fighandles = channels_to_rasterplots(experiment.filename,summarise_channels,split_raster_plots,split_figures,trialrange,x_ax_lims,condition_name, condition_units);

    % save figure
    if save_figs
        for a = 1:length(fighandles)
            figure(fighandles(a))
            print(gcf,[experiment_plot_folder filesep 'Rasterplot stimulator ' num2str(a) ' chans ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng',['-r' num2str(figure_dpi)])
        end
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

% visual representation of responses across channels
if qchannelim
    whisk_resps     = experiment.whisk_peak_rate;
    P_whisk_resps   = whisk_resps(P_whisk_stim:2:end,:);
    A_whisk_resps   = whisk_resps(A_whisk_stim:2:end,:);
    
    P_whisk_resps_LED   = whisk_resps(stimulator == 1 & LEDtime,:);
    P_whisk_resps_nLED  = whisk_resps(stimulator == 1 & ~LEDtime,:);
    
    A_whisk_resps_LED   = whisk_resps(stimulator == 2 & LEDtime,:);
    A_whisk_resps_nLED  = whisk_resps(stimulator == 2 & ~LEDtime,:);
    
    LED_resps           = experiment.LED_rate;
    
    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'Units','normalized')
    set(gcf,'Position',[.2 .2 .6 .6])
    
    subplot(1,7,1)
    imagesc(P_whisk_resps_LED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Principal whisker + LED')
    ylabel('Channel number')
    colorbar
    % axis off
    subplot(1,7,2)
    imagesc(P_whisk_resps_nLED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Principal whisker NO LED')
    ylabel('Channel number')
    colorbar
    % axis off
    
    subplot(1,7,3)
    imagesc(A_whisk_resps_LED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Adjacent whisker + LED')
    colorbar
    
    % axis off
    subplot(1,7,4)
    imagesc(A_whisk_resps_nLED')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('Adjacent whisker NO LED')
    colorbar
    % axis off
    
    subplot(1,7,5)
    imagesc(log(P_whisk_resps_LED' ./ A_whisk_resps_LED'))
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('P/A ratio WITH LED')
    colorbar
    % axis off
    
    subplot(1,7,6)
    imagesc(log(P_whisk_resps_nLED' ./ A_whisk_resps_nLED'))
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('P/A ratio NO LED')
    colorbar
    % axis off
    
    subplot(1,7,7)
    imagesc(LED_resps')
    colormap gray
    set(gca,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
    title('LED resp size')
    colorbar
    % axis off
    
    [r_P_whisk_LED, p_P_whisk_LED] = corr(mean(LED_resps)',mean(P_whisk_resps)')
    
    
    [r_A_whisk_LED, p_A_whisk_LED] = corr(mean(LED_resps)',mean(A_whisk_resps)')
    
end

