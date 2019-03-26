%% AVK POM

% 2018_09_10
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Timing pre/2018_09_10/2018_09_10-4-Timing pre.mat

% 2018_12_03
% '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Timing_pre/2018_12_03/2018_12_03-17-Timing_pre.mat' -- > beauty!
% '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Timing_pre/2018_12_03/2018_12_03-18-Timing_pre.mat'
% '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Timing_pre/2018_12_03/2018_12_03-19-Timing_pre.mat'

% 15/02/2019
% '/Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Timing_pre/2019_02_15/2019_02_15-18-Timing_pre.mat'
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV PoM/Timing/2019_02_15/2019_02_15-22-Timing.mat


%% AVK Cx M1
% 2019_01_22
% /Volumes/Akermanlab/Joram/AVK RBSN rAAV M1/Timing_pre/2019_01_22/2019_01_22-20-Timing_pre.mat


%% GTRB
% 2018_10_02
% '/Volumes/Akermanlab/Joram/AVK GTRB/Timing_LED_Power/2018_10_02/2018_10_02-5-Timing_LED_Power.mat' -->
% '/Volumes/Akermanlab/Joram/AVK GTRB/Timing_LED_Power/2018_10_02/2018_10_02-6-Timing_LED_Power.mat'

% /Volumes/Akermanlab/Joram/AVK GTRB/Timing_LED_Power/2018_10_16/2018_10_16-3-Timing_LED_Power.mat --> also good, 10ms whisk stim
% /Volumes/Akermanlab/Joram/AVK GTRB/Timing_LED_Power/2018_10_16/2018_10_16-5-Timing_LED_Power.mat--> 25ms whisk stim



%% Input variables 
data_file       = ['/Volumes/Akermanlab/Joram/AVK GTRB/Timing_LED_Power/2018_10_02/2018_10_02-5-Timing_LED_Power.mat'];

q_reload        = 1;

%% PSTH settings
early_resp_win  = [0.006 0.03]; % response window for assessing 'early' response (spike count is from this time window)
late_resp_win   = [0.03 0.1];   % response window for assessing 'late' response (spike count is from this time window)
opto_resp_win   = [0 0.03];     % window for assessing opto response

psth_win        = [-0.2 .5];   % window for plotting PSTH
psth_bin_size   = [0.001];  % bin size for PSTH

do_max_power_only = true;

do_one_whisker_only = true;
target_whisker_nr   = 1;

channels        = [1:16];    % which recording sites / channels to include


%% Code execution starts here

close all

n_channels      = length(channels);


if q_reload
    load(data_file)
end

early_spike_count   = [];
late_spike_count    = [];
opto_spike_count    = [];
delta_t             = [];
whisker_nr          = [];
opto_power          = [];
n_trials            = [];
all_psth_counts    	= [];
for a = 1:length(ephys_data.conditions)
    
        this_cond               = ephys_data.conditions(a);

        this_t_whisk            = this_cond.whisk_onset;
        this_t_opto             = this_cond.LED_onset;
        this_whisker_nr       	= this_cond.whisk_stimulator;
        this_n_trials           = this_cond.n_trials;
        
        spikes                  = this_cond.spikes(channels, :, :) - this_t_whisk;
        
        delta_t(a)              = this_t_opto - this_t_whisk;
        whisker_nr(a)           = this_whisker_nr;
        opto_power(a)           = this_cond.LED_power;
        
        spikes_in_early_win     = spikes > early_resp_win(1) & spikes <= early_resp_win(2);
        spikes_in_late_win      = spikes > late_resp_win(1) & spikes <= late_resp_win(2);
        
        early_spike_count(a)    = sum(spikes_in_early_win(:)) / n_channels / this_n_trials;
        late_spike_count(a)     = sum(spikes_in_late_win(:)) / n_channels / this_n_trials;
        
        n_trials(a)           	= this_n_trials;
        
        % get full post_stimulus_time_histogram
        figure
        [plot_handle, psth_counts, psth_bins]  = psth(spikes, psth_bin_size, psth_win);
        close(gcf)
        
        all_psth_counts         = [all_psth_counts; psth_counts(:)'];
        
        %% opto resp
        
        spikes                  = (spikes + this_t_whisk) - this_t_opto;
        spikes_in_opto_win      = spikes > opto_resp_win(1) & spikes <= opto_resp_win(2);
        opto_spike_count(a)     = sum(spikes_in_opto_win(:)) / n_channels / this_n_trials;
end
% 
% opto_power(isnan(opto_power))   = 2.5;
% delta_t(isnan(delta_t))         = 2.5;


uniq_delta_ts       = unique(delta_t);
uniq_whisker_nrs    = unique(whisker_nr);
uniq_opto_powers 	= unique(opto_power);

if do_max_power_only
    uniq_opto_powers    = max(uniq_opto_powers);
end

if do_one_whisker_only
    uniq_whisker_nrs    = target_whisker_nr;
end

q_control         	= uniq_delta_ts > 1;
uniq_delta_ts       = uniq_delta_ts(uniq_delta_ts <= 0 & ~q_control);

uniq_delta_ts    	= uniq_delta_ts(~isnan(uniq_delta_ts));
uniq_whisker_nrs  	= uniq_whisker_nrs(~isnan(uniq_whisker_nrs));
uniq_opto_powers   	= uniq_opto_powers(~isnan(uniq_opto_powers));


n_delta_ts          = length(uniq_delta_ts);
n_whiskers          = length(uniq_whisker_nrs);
n_opto_powers       = length(uniq_opto_powers);

for i = 1:n_whiskers
    for j = 1:n_opto_powers
        
        q_whisker_nr  	= whisker_nr == uniq_whisker_nrs(i);
        q_opto_power    = opto_power == uniq_opto_powers(j);
        
        cond_early_spike_count      = early_spike_count(q_whisker_nr & q_opto_power);
        cond_late_spike_count       = late_spike_count(q_whisker_nr & q_opto_power);
        cond_delta_t                = delta_t(q_whisker_nr & q_opto_power);
        
        q_below_0                   = cond_delta_t <= 0;
        
        subplot_ind     = ((i - 1) * n_opto_powers) + j;
        
        figure(1)
        subplot(n_whiskers,n_opto_powers,subplot_ind)
        plot(cond_delta_t(q_below_0),cond_early_spike_count(q_below_0),'k.-','LineWidth',2,'MarkerSize',20)
        xlimits     = xlim;
        line([xlimits],[cond_early_spike_count(end) cond_early_spike_count(end)],'Color',[1 0 0],'LineWidth',2,'LineStyle',':')
        title('Early response window')
        ylabel('Spike count')
        xlabel('Opto-whisk time delay')
        fixplot
        yzero
        
        figure(2)
        subplot(n_whiskers,n_opto_powers,subplot_ind)
        plot(cond_delta_t(q_below_0),cond_late_spike_count(q_below_0),'k.-','LineWidth',2,'MarkerSize',20)
        xlimits     = xlim; 
        line([xlimits],[cond_late_spike_count(end) cond_late_spike_count(end)],'Color',[1 0 0],'LineWidth',2,'LineStyle',':')
        title('Late response window')
        ylabel('Spike count')
        xlabel('Opto-whisk time delay')
        fixplot
        yzero
        
        for k = 1:n_delta_ts
            
            this_delta_t    = uniq_delta_ts(k);
            q_delta_t       = delta_t == this_delta_t;
            
            subplot_ind2    = ((k - 1) * n_opto_powers) + j;
            
            figure(2+i)
            subplot(n_delta_ts, n_opto_powers, subplot_ind2)
            
            this_psth       = mean(all_psth_counts(q_delta_t & q_whisker_nr & q_opto_power,:),1);
            
            this_n_trials   = n_trials(q_delta_t & q_whisker_nr & q_opto_power);
            
            plot(psth_bins(1:end-1),this_psth / n_channels / this_n_trials / psth_bin_size,'k-','LineWidth',2)
            title(['Opto-whisk delay = ' num2str(-this_delta_t * 1000) 'ms'])
            ylabel('Spike rate (Hz)')
            
            % fixplot
            ylim([0 1000])
            
        end
        
        
    end
end
figure(1)
subplot_equal_y

figure(2)
subplot_equal_y

figure(3)
subplot_equal_y
set(gcf,'Units','normalized','Position',[.3 0 .4 1])

figure(4)
subplot_equal_y
set(gcf,'Units','normalized','Position',[.3 0 .4 1])




% Opto_power
q_control = delta_t > 1;
counts_by_opto_power = [];
for l = 1:length(uniq_opto_powers)
    q_opto_power = opto_power == uniq_opto_powers(l);
    
    counts_by_opto_power(l)     = mean(opto_spike_count(q_opto_power & q_control));
    

end

figure(5)
plot(uniq_opto_powers,counts_by_opto_power,'k.-','LineWidth',2,'MarkerSize',20)
fixplot
yzero

