function[Stim_response] = Stim_Responsive(spikes,resp_win,control_win,n_trials,delta_t, plot_figs,stim_type,fn)
% Stim_Responsive    checks the response of a sorted unit to a stimulus response
% Inputs
% spikes : N (units) x M (trials) x C (spike times) Matrix of spikes
% Response window: window in Ms to evaluate spike response following
% stimulus at time point 0
% Control window : control time point to assess spike response against,
% preceeds stimulus onset
%
%plot_figs : should figures be plotted. Default true.
% stim_type : simply for figure output, figures will be labelled according
%             to stim_type and saved with file name accordingly. 
% fn :  directory where output figures should be stored. 
%
% Alexander von klemperer 2020

% Default to all units
    if nargin < 6 || isempty(plot_figs)
        plot_figs        = true;
    end

    psth_bins       = [-0.3:0.001:0.3]; % sets bin size in seconds to be displayed in raster plots


% Binned spike rate
    [spike_rates]        = spike_rates_individual(spikes, resp_win); 
    [spike_probs]        = spike_prob_by_channel(spikes, resp_win);
    
    
    [spike_rates_control]        = spike_rates_individual(spikes, control_win);
    [spike_probs_control]        = spike_prob_by_channel(spikes, control_win);
    
    
    for k = 1 : numel(spike_probs) % takes spike probablities from all units
        [p(k),h(k)] = Chi2_test(spike_probs(k).*n_trials,n_trials,spike_probs_control(k).*n_trials,n_trials,0.05); % Statistically tests whether the probalitiy of spiking in post stimulus period is different than control period
        if (h(k) == 1) && spike_probs(k)>spike_probs_control(k) % if null hypothesis rejected and spiking probability greater than control
            Good_resp_prob(k) = 1; 
        else
            Good_resp_prob(k) = 0;
         end;
     
         if (kstest(spike_rates(k,:))==0) && (kstest(spike_rates_control(k,:))==0); % tests whether spike rates across all trials are normally distributed
           [Good_resp_rate(k),pt(k)] = ttest2(spike_rates(k,:),spike_rates_control(k,:),'Tail','right'); % tests whether spike rates are larger following stimulus than control
         else
         [pt(k),Good_resp_rate(k)] = signrank(spike_rates(k,:),spike_rates_control(k,:),'Tail','right');
        end;
         
     
    end;
        %sets output structure
        Stim_response.Good_resp_prob = Good_resp_prob;
        Stim_response.Good_resp_rate = Good_resp_rate;
        Stim_response.Spike_rates = spike_rates;
        Stim_response.Spike_probs = spike_probs;
        Stim_response.Spike_rates_control = spike_rates_control;
        Stim_response.Spike_probs_control = spike_probs_control;
        
   if (plot_figs == true)    
    % Raster plots
    Response_plot_h         	= figure('Name',stim_type,'NumberTitle','off');
    set(gcf,'Units','Normalized','Position',[.3 0 .4 1],'PaperPositionMode','auto')
  
    figure(Response_plot_h)
    subplot(2,2,1)
    raster_plot(spikes,1);
    xlim([min(psth_bins) max(psth_bins)])
    title(['Opto-whisk delay = ' num2str(-delta_t * 1000) 'ms'])
    ylabel('Unit number')
    set(gca,'FontName','Helvetica','FontWeight','Bold','box','off')
    
    subplot(2,2,2)
    raster_plot(spikes(1,:,:),2);
    xlim([min(psth_bins) max(psth_bins)])
    title(['Opto-whisk delay = ' num2str(-delta_t * 1000) 'ms'])
    ylabel('Trial number')
    set(gca,'FontName','Helvetica','FontWeight','Bold','box','off')
    
    
       % Set sensible axis and labels
    title(['Opto-whisk delay = ' num2str(-delta_t * 1000) 'ms'])
    ylabel('Spike count')
    
    set(gca,'FontName','Helvetica','FontWeight','Bold','box','off')
    
    subplot(2,2,3);
    hold on;
    plot(nanmean(spike_rates,2),'g^');
    plot(nanmean(spike_rates_control,2),'go');
    plot(Good_resp_rate*1.1*(max(nanmean(spike_rates,2))), 'g*');
      title('Rates');
    
    subplot(2,2,4);
    hold on;
    plot(spike_probs,'b^');
    plot(spike_probs_control,'bo');
    plot(Good_resp_prob,'b*')
    title('Probability');
    
    if ~exist(fn,'dir')
        mkdir(fn)
    end;
    saveas(Response_plot_h,[fn '\' stim_type '_response_plots.fig']);
    saveas(Response_plot_h,[fn '\' stim_type '_response_plots.png']);
    
   end;