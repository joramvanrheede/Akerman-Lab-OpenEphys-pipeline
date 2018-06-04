% Multiple experiment plots for drive experiment

%% Add LED response size
% LED ON response peak
% LED ON sustained response (in control condition!)
% LED OFF response
% some responsiveness criterion implementation

%% Add examination of tracking repeat stimuli
% Plot response with error bars to first 4-5 stimuli
% Plot response to first stimulus vs response to stimulus 3 seconds later
% (by frequency)
% Make some metrics: first/fourth stim ratio, first vs 3 secs ratio
%

save_folder         = '/Users/Joram/Dropbox/Akerman Postdoc/Figures/Matlab output';
save_expt_name      = 'Frequency multi 20180502';
save_figs        	= false;

summarise_channels  = [1:16]; % include these channels

% apply threshold for LED responsiveness?
q_check_LED_resp    = true; % 
LEDresp_threshold   = 3;  % threshold relative to spontaneous

% response window for assessing whisk response to stimuli for purpose of
% comparison of 1st, 2nd, 3rd, last...
respwinsize         = 0.04;

% Zoomed in values for visualising instantaneous LED response
z_LED_tmin          = -.01;
z_LED_tmax          = .2;

%% Set for this type of experiment
split_conditions    = [1 6]; % split by these conditions, summarise over others
split_plots         = [6]; % [4 6] works

%%
close all

P_rate_LED_on       = [];
P_rate_LED_off      = [];
A_rate_LED_on       = [];
A_rate_LED_off      = [];
PA_ratio_LED_on     = [];
PA_ratio_LED_off    = [];
P_LED_onoff_ratio   = [];
A_LED_onoff_ratio   = [];
P_pktime_LED_on     = [];
P_pktime_LED_off	= [];
A_pktime_LED_on     = [];
A_pktime_LED_off	= [];
LED_rel_mean        = [];
LED_sust_mean       = [];
LED_OFF_mean       	= [];
LED_rate_traces   	= [];
LED_win_edges       = [];
counter             = 0;
LEDresp           	= [];
freq_resp_measures  = [];
for i = 1:length(sdata)
    experiment              = sdata(i).expt;
    
    for j = 1:length(experiment)

        counter             = counter+1;
        
        condition_mat       = experiment(j).condition_mat
        
        first_resp_peaks    = [];
        timed_resp_peaks    = [];
        for k = 1:size(condition_mat,1)
            
            these_conds         = condition_mat(k,:);
            
            this_LED_delay      = these_conds(1);
            this_frequency      = these_conds(4);
            this_stim           = these_conds(6);
            
            stim_times          = [0:1/round(this_frequency):5.5]; % recreate stimulus times from frequency
            
            resp_peaks          = [];
            for l = 1:length(stim_times)
                
                % determine response assessment window for this stimulus time
                qstimwin            = experiment(j).whiskwinedges > stim_times(l) & experiment(j).whiskwinedges < (stim_times(l) + respwinsize);
                
                % get segment of whisk_profile to determine max instantaneous rate of fire in this window
                profile_segment     = squeeze(mean(mean(experiment(j).whisk_profile(k,summarise_channels,qstimwin),1),2)); 
                
                % if we've requested data beyond whisker profile, use NaN
                % rather than empty values
                if isempty(profile_segment)
                    profile_segment = NaN;
                end
                
                % find max instantaneous ROF
                resp_peaks(l)  = max(profile_segment);
            end
            
            first_stim_times    = [stim_times(1:4)];
            first_resp_peaks    = [first_resp_peaks; resp_peaks(1:4)];
            
            timed_stim_times    = [0:1:3];
            timed_resp_peaks    = [timed_resp_peaks; resp_peaks(1) resp_peaks(find(round(stim_times*100) == 100)) resp_peaks(find(round(stim_times*100) == 200)) resp_peaks(find(round(stim_times*100) == 300))];

        end
        
        freq_resp_measures(counter).first_resp_peaks    = first_resp_peaks;
        freq_resp_measures(counter).timed_resp_peaks    = timed_resp_peaks;
        freq_resp_measures(counter).firstratio          = first_resp_peaks(:,4)./first_resp_peaks(:,1);
        freq_resp_measures(counter).timedratio          = timed_resp_peaks(:,4)./timed_resp_peaks(:,1);
        freq_resp_measures(counter).frequencies         = condition_mat(:,4);
        freq_resp_measures(counter).LEDtime             = condition_mat(:,1);
        freq_resp_measures(counter).stimulator          = condition_mat(:,6);

        
        %% Look at peak instantaneous firing rate of first response between LED and no LED
        contrast_conds      = [1 6];
        contrast_cond_mat   = condition_mat(:,contrast_conds);
        [contrast_cond_rows, indxa, cond_inds] = unique(contrast_cond_mat,'rows');
        contrast_rates      = [];
        contrast_times      = [];
        for a  = unique(cond_inds)'
            contrast_rates      = [contrast_rates; mean(mean(experiment(j).whisk_peak_rate(cond_inds == a,summarise_channels)))];
            contrast_times      = [contrast_times; mean(mean(experiment(j).whisk_peak_time(cond_inds == a,summarise_channels)))];
        end
        
        if contrast_rates(3) > contrast_rates(4)
            P_whisk_stim = 1;
            A_whisk_stim = 2;
        else
            P_whisk_stim = 2;
            A_whisk_stim = 1;
            % make sure that principal stimulator is stimulator 1 in freq_resp_measures 
            freq_resp_measures(counter).stimulator = 3 -freq_resp_measures(counter).stimulator;
        end
        
        P_rate_LED_on(counter)   	= contrast_rates(P_whisk_stim);
        P_rate_LED_off(counter)  	= contrast_rates(P_whisk_stim+2);
        A_rate_LED_on(counter)  	= contrast_rates(A_whisk_stim);
        A_rate_LED_off(counter)   	= contrast_rates(A_whisk_stim+2);
        
        PA_ratio_LED_on(counter)  	= P_rate_LED_on(counter) / A_rate_LED_on(counter);
        PA_ratio_LED_off(counter)  	= P_rate_LED_off(counter) / A_rate_LED_off(counter);
        
        P_LED_onoff_ratio(counter) 	= P_rate_LED_on(counter) / P_rate_LED_off(counter);
        A_LED_onoff_ratio(counter) 	= A_rate_LED_on(counter) / A_rate_LED_off(counter);
        
        
        P_pktime_LED_on(counter)  	= contrast_times(P_whisk_stim);
        P_pktime_LED_off(counter)  	= contrast_times(P_whisk_stim+2);
        A_pktime_LED_on(counter)  	= contrast_times(A_whisk_stim);
        A_pktime_LED_off(counter)  	= contrast_times(A_whisk_stim+2);
        
        
        %% LED resp values
        
        LEDresp(counter)         	= mean(experiment(j).LED_rel(:));
        
        LED_delays                  = condition_mat(:,1);
        [uniq_LED_delays, indxa, LED_inds] = unique(LED_delays);
        LED_control_inds            = LED_inds == LED_inds(end);
        
        LED_rate_mat             	= experiment(j).LED_rate(LED_control_inds,summarise_channels);
        LED_mean(counter)           = mean(LED_rate_mat(:));
        
        LED_sust_mat                = experiment(j).LED_sust_rates(LED_control_inds,summarise_channels);
        LED_sust_mean(counter)      = mean(LED_sust_mat(:));
        
        LED_OFF_mat                 = experiment(j).LED_OFF_rates(LED_control_inds,summarise_channels);
        LED_OFF_mean(counter)       = mean(LED_OFF_mat(:));
        
        %% LED PSTHs
        LED_rate_PSTHs            	= experiment(j).LED_win_rates(LED_control_inds,summarise_channels,:);
        mean_LED_rate_PSTHs      	= squeeze(mean(LED_rate_PSTHs,2));
        LED_rate_traces(:,counter) 	= mean(mean_LED_rate_PSTHs);
        LED_win_edges(:,counter)   	= experiment(j).LEDwinedges(1:end-1);
        
    end
end

qLEDresp            = LEDresp > LEDresp_threshold;

if q_check_LED_resp
    P_rate_LED_on       = P_rate_LED_on(qLEDresp);
    P_rate_LED_off      = P_rate_LED_off(qLEDresp);
    A_rate_LED_on       = A_rate_LED_on(qLEDresp);
    A_rate_LED_off      = A_rate_LED_off(qLEDresp);
    
    PA_ratio_LED_on     = PA_ratio_LED_on(qLEDresp);
    PA_ratio_LED_off	= PA_ratio_LED_off(qLEDresp);
    
    P_LED_onoff_ratio   = P_LED_onoff_ratio(qLEDresp);
    A_LED_onoff_ratio	= A_LED_onoff_ratio(qLEDresp);
    
    
    P_pktime_LED_on     = P_pktime_LED_on(qLEDresp);
    P_pktime_LED_off  	= P_pktime_LED_off(qLEDresp);
    A_pktime_LED_on     = A_pktime_LED_on(qLEDresp);
    A_pktime_LED_off	= A_pktime_LED_off(qLEDresp);
    
    freq_resp_measures  = freq_resp_measures(qLEDresp);
    
end

%%
first_pks   = cell2mat({freq_resp_measures.first_resp_peaks}');
timed_pks   = cell2mat({freq_resp_measures.timed_resp_peaks}');
first_ratio = cell2mat({freq_resp_measures.firstratio}');
timed_ratio = cell2mat({freq_resp_measures.timedratio}');
frequencies = round(cell2mat({freq_resp_measures.frequencies}'));
LEDtime     = cell2mat({freq_resp_measures.LEDtime}');
LEDtime     = LEDtime < 2;
stimulator  = cell2mat({freq_resp_measures.stimulator}');

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

figure(22)
set(gcf,'Units','Normalized')
set(gcf,'Position',[.2 .4 .6 .4])
% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);



%% Response size plotting

figure
set(gcf,'Units','normalized')
set(gcf,'Position',[0 .4 1 .4])
set(gcf,'Color',[1 1 1])

subplot(1,3,1)
pairedlineplot(P_rate_LED_off,P_rate_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal Whisker ROF')
title('Principal whisker response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

subplot(1,3,2)
pairedlineplot(A_rate_LED_off,A_rate_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Adjacent Whisker ROF')
title('Adjacent whisker response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylim([0 ylimits(2)*1.2])

subplot(1,3,3)
pairedlineplot(PA_ratio_LED_off,PA_ratio_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal/Adjacent ratio')
title('Principal/Adjacent ratio')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

if save_figs
    experiment_plot_folder      = [save_folder filesep save_expt_name];
    if ~isdir(experiment_plot_folder)
        mkdir(experiment_plot_folder)
    end
    print(gcf,[experiment_plot_folder filesep 'Principal vs Adjacent resp paired plots - chan ' num2str(summarise_channels(1)) '-' num2str(summarise_channels(end))],'-dpng','-r300')
end

figure
set(gcf,'Color',[1 1 1])
pairedlineplot(P_LED_onoff_ratio,A_LED_onoff_ratio,{'Principal' 'Adjacent'},'Whisker','LED effect on whisker')
title('Effect of LED on P vs A whisker')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

%% Peak time plotting

figure
set(gcf,'Units','normalized')
set(gcf,'Position',[.167 .4 .67 .4])
set(gcf,'Color',[1 1 1])

subplot(1,2,2)
pairedlineplot(A_pktime_LED_off,A_pktime_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Adjacent Whisker peak time')
title('Adjacent whisker peak time')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

subplot(1,2,1)
pairedlineplot(P_pktime_LED_off,P_pktime_LED_on,{'LED OFF' 'LED ON'},'Opto condition','Principal Whisker peak time')
title('Principal whisker peak time')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
ylim([0 ylimits(2)*1.2])


%% LED plotting
figure
set(gcf,'Units','normalized')
set(gcf,'Position',[.167 .4 .67 .4])
set(gcf,'Color',[1 1 1])

subplot(1,2,1)
plot_handle = plot(LED_win_edges(:,~qLEDresp),LED_rate_traces(:,~qLEDresp),'r-','LineWidth',2);
hold on
plot_handle = plot(LED_win_edges(:,qLEDresp),LED_rate_traces(:,qLEDresp),'k-','LineWidth',2);
title('Average LED response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
xlabel('Time (s)')
ylabel('Spike rate relative to spontaneous')
axis tight
ylimits = ylim;
ylim([0 ylimits(2)*1.2])



qzoom           = LED_win_edges(:,1) >= z_LED_tmin & LED_win_edges(:,1) <= z_LED_tmax;
z_LED_tvals     = LED_win_edges(qzoom,1);
z_LED_traces    = LED_rate_traces(qzoom,:);

subplot(1,2,2)
plot(z_LED_tvals,z_LED_traces(:,~qLEDresp),'r-','LineWidth',2)
hold on
plot(z_LED_tvals,z_LED_traces(:,qLEDresp),'k-','LineWidth',2)
title('Average LED response')
set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20,'FontWeight','Bold')
xlabel('Time (s)')
ylabel('Spike rate relative to spontaneous')
axis tight
ylimits = ylim;
ylim([0 ylimits(2)*1.2])

% subplot(1,2,2)
% pairedlineplot(LED_rel_mean,LED_sust_rel_mean,{'LED ON' 'LED SUSTAINED'},'Opto dynamics','Spike rate relative to mean')
% title('Adjacent whisker peak time')
% set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
% ylimits = ylim;
% ylim([0 ylimits(2)*1.2])


