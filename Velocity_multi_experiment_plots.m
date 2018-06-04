% Velocity_multi_experiment_plots
%

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
save_expt_name      = 'Velocity multi 08052018';
save_figs        	= false;

summarise_channels  = [6]; % include these channels

% apply threshold for LED responsiveness?
q_check_LED_resp    = true; % 
LEDresp_threshold   = 3;  % threshold relative to spontaneous

%% Set for this type of experiment
split_conditions    = [1 5]; % split by these conditions, summarise over others
split_plots         = [6]; % [4 6] works

%%
close all

LED_rel_mean        = [];
LED_sust_mean       = [];
LED_OFF_mean       	= [];
LED_rate_traces   	= [];
LED_win_edges       = [];
counter             = 0;
LEDresp           	= [];
vel_resp_measures   = [];
for i = 1:length(sdata)
    experiment              = sdata(i).expt;
    
    for j = 1:length(experiment)
        counter             = counter+1; % keep track of total nr of experiments
        
        condition_mat       = experiment(j).condition_mat;
        
        peak_rates  = [];
        velocities  = [];
        stimulator  = [];
        LEDtimes    = [];
        for k = 1:size(condition_mat,1)
            
            these_conds         = condition_mat(k,:);
            
            LEDtimes(k)         = these_conds(1);
            velocities(k)       = these_conds(5);
            stimulator(k)     	= these_conds(6);
            
            peak_rates(k)       = mean(experiment(j).whisk_peak_rate(k,summarise_channels));
            
        end
        
        stim_1_resp     = mean(peak_rates(stimulator == 1));
        stim_2_resp     = mean(peak_rates(stimulator == 2));
        
        if stim_1_resp > stim_2_resp
            P_whisk_stim = 1;
            A_whisk_stim = 2;
            vel_resp_measures(counter).stimulator   = stimulator';
        else
            P_whisk_stim = 2;
            A_whisk_stim = 1;
            vel_resp_measures(counter).stimulator   = 3-stimulator';
        end
        
        vel_resp_measures(counter).peak_rates   = peak_rates';
        vel_resp_measures(counter).velocities   = velocities';
        vel_resp_measures(counter).LEDtime      = LEDtimes';
        vel_resp_measures(counter).PA_ratios    = [peak_rates(stimulator == P_whisk_stim) ./ peak_rates(stimulator == A_whisk_stim)]';

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

qLEDresp    = LEDresp > LEDresp_threshold;

if q_check_LED_resp
 
    vel_resp_measures   = vel_resp_measures(qLEDresp);
    
end

peak_rates  = cell2mat({vel_resp_measures.peak_rates}');
velocities  = round(cell2mat({vel_resp_measures.velocities}'));
PA_ratios   = cell2mat({vel_resp_measures.PA_ratios}');
LEDtime     = cell2mat({vel_resp_measures.LEDtime}');
LEDtime     = LEDtime < 2;
stimulator  = cell2mat({vel_resp_measures.stimulator}');

uniqvels        = unique(velocities);
uniqLEDs        = unique(LEDtime);
uniqstims       = unique(stimulator);
nLEDconds       = length(uniqLEDs);
nstimulators    = length(uniqstims);

figure
for a = 1:nLEDconds
    for b = 1:nstimulators
        qLED    = LEDtime == uniqLEDs(a);
        qstim   = stimulator == uniqstims(b);
        
        % set up variables for loop over uniqfreqs
        plot_rate_means         = [];
        plot_rate_serrs         = [];
        plot_vels               = [];
        
        % loop over uniqvels
        for c = 1:length(uniqvels)
            qvel   = velocities == uniqvels(c);
            qall    = qvel & qLED & qstim;
            
            % if no data meet these criteria, go to next iteration
            if sum(qall) == 0
                continue
            end
            
            % get ratios that meet criteria
            these_peak_rates        = peak_rates(qall);
            
            % get means and standard errors for plotting
            plot_rate_means    = [plot_rate_means; mean(these_peak_rates)];
            plot_rate_serrs    = [plot_rate_serrs; serr(these_peak_rates)];

            plot_vels        	= [plot_vels, uniqvels(c)];
        end
        
        figure(21)
        subplot(1,3,b)
        if ~uniqLEDs(a)
            errorbar(log(plot_vels),plot_rate_means,plot_rate_serrs,'k-','LineWidth',2)
        else
            errorbar(log(plot_vels),plot_rate_means,plot_rate_serrs,'r-','LineWidth',2)
        end
        xlim([0 6])
        set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
        if b == 1
            title('Principal whisker, Peak response rate')
        elseif b == 2
            title('Adjacent whisker, Peak response rate')
        end
        xlabel('Log Velocity')
        ylabel('Peak response rate')
        hold on
        
    end
    
    plot_PAR_means  = [];
    plot_PAR_serrs  = [];
    plot_vels       = [];
    for c = 1:length(uniqvels)
        qvel    = velocities == uniqvels(c);
        qall    = qvel & qLED;
        qall    = qall(stimulator == 1); 
        
        % if no data meet these criteria, go to next iteration
        if sum(qall) == 0
            continue
        end
        
        % get ratios that meet criteria
        these_PA_ratios  	= PA_ratios(qall);
        
        % get means and standard errors for plotting
        plot_PAR_means      = [plot_PAR_means; mean(these_PA_ratios)];
        plot_PAR_serrs      = [plot_PAR_serrs; serr(these_PA_ratios)];
        
        plot_vels        	= [plot_vels, uniqvels(c)];
    end
    
    subplot(1,3,3)
    if ~uniqLEDs(a)
        errorbar(log(plot_vels),plot_PAR_means,plot_PAR_serrs,'k-','LineWidth',2)
    else
        errorbar(log(plot_vels),plot_PAR_means,plot_PAR_serrs,'r-','LineWidth',2)
    end
    xlim([0 6])
    set(gca,'LineWidth',2,'FontName','Garamond','FontSize',20)
    xlabel('Log Velocity')
    ylabel('P-A Ratio')
    title('P-A Ratio')
    hold on
    
end
figure(21)
set(gcf,'Units','Normalized')
set(gcf,'Position',[.1 .4 .8 .4])
% Set all y axes to the same range (based on the largest range)
plotaxes    = get(gcf,'Children');
maxy        = cellfun(@max,get(plotaxes,'Ylim'));
set(plotaxes,'Ylim',[0 max(maxy)]);
ylim([0 3]); % adjust final axes

