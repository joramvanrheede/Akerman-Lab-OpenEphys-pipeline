

group_folders   = {'F:\Sorted_Cortical\Timing'};
group_names     = {'Cortical'};

target_power    = 50;
target_delay    = [-0.02];

min_depth       = 0;
max_depth       = 200;

q_reload        = 1;

tic % time script


%% function group_data = get_unit_data(group_folders) ?
if q_reload
    
    groups = load_sorted_experiments(group_folders);
    
    clear group_data
    
    for a = 1:length(groups)
        group_data(a).unit_depths               = [];
        group_data(a).delta_t                   = [];
        
        group_data(a).spike_rate                = [];
        group_data(a).spike_rate_p              = [];
        group_data(a).spike_rate_ctrl           = [];
        
        group_data(a).spike_prob                = [];
        group_data(a).spike_prob_p              = [];
        group_data(a).spike_prob_ctrl           = [];
        
        group_data(a).first_spike_time          = [];
        group_data(a).first_spike_p             = [];
        group_data(a).first_spike_time_ctrl     = [];
        
        group_data(a).first_spike_jitter        = [];
        group_data(a).first_spike_jitter_ctrl   = [];
        
        group_data(a).peak_spike_rate           = [];
        group_data(a).peak_spike_rate_ctrl    	= [];
        
        group_units(a).spike_rate               = [];
        group_units(a).depths                   = [];
        
        group_data(a).select_power              = [];
        
        cum_n_units                             = 0;
        for b = 1:length(groups(a).prep)
            
            opto_power = [];
            for c = 1:length(groups(a).prep(b).expt_data)
                
                ephys_data      = groups(a).prep(b).expt_data(c).ephys_data;
                
                opto_power(c)   = mode([ephys_data.conditions.LED_power]);
                
            end
            
            target_power_diff   = abs(target_power - opto_power);
            target_power_ind    = find(target_power_diff == min(target_power_diff),1);
            select_power        = opto_power(target_power_ind);
            
            
            ephys_data          = groups(a).prep(b).expt_data(target_power_ind).ephys_data;
            
            delta_t             = [ephys_data.conditions.LED_onset] - [ephys_data.conditions.whisk_onset];
            delta_t             = round(delta_t,3); % round to nearest ms
            target_delta_t_ind 	= find(delta_t == target_delay,1);
            

            timing_unit_data                    = timing_by_unit(ephys_data);
            close all
            n_units                             = size(ephys_data.conditions(1).spikes,1);
            next_unit                           = n_units + 1;

            unit_inds                           = cum_n_units + [1:n_units];
            cum_n_units                         = cum_n_units + n_units;
            
            n_dts                               = size(timing_unit_data.spike_rate,2);
            max_dts                             = max(size(timing_unit_data.spike_rate,2),size(group_data(a).spike_rate,2));
            
            %% pre-allocate to max with NaNs
            group_data(a).spike_rate(unit_inds,1:max_dts)           = NaN;
            group_data(a).spike_rate_p(unit_inds,1:max_dts)         = NaN;
            
            group_data(a).spike_prob(unit_inds,1:max_dts)           = NaN;
            group_data(a).spike_prob_p(unit_inds,1:max_dts)         = NaN;
            
            group_data(a).first_spike_time(unit_inds,1:max_dts) 	= NaN;
            group_data(a).first_spike_p(unit_inds,1:max_dts)    	= NaN;
            
            group_data(a).first_spike_jitter(unit_inds,1:max_dts)	= NaN;
            
            group_data(a).peak_spike_rate(unit_inds,1:max_dts)  	= NaN;
            
            group_data(a).unit_depths(unit_inds)                    = NaN;
            group_data(a).delta_t(unit_inds,1:max_dts)          	= NaN;
            group_data(a).select_power(unit_inds)                	= NaN;
            
            %%
            group_data(a).spike_rate(unit_inds,1:n_dts)      	= timing_unit_data.spike_rate;
            group_data(a).spike_rate_p(unit_inds,1:n_dts)    	= timing_unit_data.spike_rate_p;
            group_data(a).spike_rate_ctrl(unit_inds)            = timing_unit_data.spike_rate(:,end);
            
            group_data(a).spike_prob(unit_inds,1:n_dts)      	= timing_unit_data.spike_probabilities;
            group_data(a).spike_prob_p(unit_inds,1:n_dts)       = timing_unit_data.spike_prob_p;
            group_data(a).spike_prob_ctrl(unit_inds)            = timing_unit_data.spike_probabilities(:,end);
            
            group_data(a).first_spike_time(unit_inds,1:n_dts) 	= timing_unit_data.first_spike_times;
            group_data(a).first_spike_p(unit_inds,1:n_dts)    	= timing_unit_data.first_spike_p;
            group_data(a).first_spike_time_ctrl(unit_inds)    	= timing_unit_data.first_spike_times(:,end);
            
            group_data(a).first_spike_jitter(unit_inds,1:n_dts)	= timing_unit_data.first_spike_jitter;
            group_data(a).first_spike_jitter_ctrl(unit_inds)    = timing_unit_data.first_spike_jitter(:,end);
            
            group_data(a).peak_spike_rate(unit_inds,1:n_dts)  	= timing_unit_data.peak_spike_rates;
            group_data(a).peak_spike_rate_ctrl(unit_inds)       = timing_unit_data.peak_spike_rates(:,end);
            
            group_data(a).unit_depths(unit_inds)              	= ephys_data.unit_depths;
            group_data(a).delta_t(unit_inds,1:n_dts)          	= repmat(delta_t,n_units,1);
            group_data(a).select_power(unit_inds)             	= select_power;
            
        end
        

        [group_data(a).unit_depths, depth_sort_inds]    = sort(group_data(a).unit_depths);
        
        group_data(a).spike_rate                        = group_data(a).spike_rate(depth_sort_inds,:);
        group_data(a).spike_rate_p                      = group_data(a).spike_rate_p(depth_sort_inds,:);
        group_data(a).spike_rate_ctrl                 	= group_data(a).spike_rate_ctrl(depth_sort_inds);
        
        group_data(a).spike_prob                        = group_data(a).spike_prob(depth_sort_inds,:);
        group_data(a).spike_prob_p                     	= group_data(a).spike_prob_p(depth_sort_inds,:);
        group_data(a).spike_prob_ctrl                 	= group_data(a).spike_prob_ctrl(depth_sort_inds);
        
        group_data(a).first_spike_time              	= group_data(a).first_spike_time(depth_sort_inds,:);
        group_data(a).first_spike_p                     = group_data(a).first_spike_p(depth_sort_inds,:);
        group_data(a).first_spike_time_ctrl            	= group_data(a).first_spike_time_ctrl(depth_sort_inds);
        
        group_data(a).first_spike_jitter                = group_data(a).first_spike_jitter(depth_sort_inds,:);
        group_data(a).first_spike_jitter_ctrl           = group_data(a).first_spike_jitter_ctrl(depth_sort_inds);
        
        group_data(a).peak_spike_rate                   = group_data(a).peak_spike_rate(depth_sort_inds,:);
        group_data(a).peak_spike_rate_ctrl              = group_data(a).peak_spike_rate_ctrl(depth_sort_inds);
        
        group_data(a).delta_t                           = group_data(a).delta_t(depth_sort_inds,:);
    end
end

figure(1)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(2)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(3)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(4)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(5)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
for a = 1:length(groups)
    unit_depths         = group_data(a).unit_depths;
    q_depth             = unit_depths <= max_depth & unit_depths >= min_depth;
    
    delta_t             = group_data(a).delta_t;
    
    %% Spike probability
    
    spike_prob          = group_data(a).spike_prob;
    spike_prob_ctrl     = group_data(a).spike_prob_ctrl;
    delta_spike_prob    = spike_prob - spike_prob_ctrl(:);
    
    spike_prob_p        = group_data(a).spike_prob_p;
    
    q_prob_p                        = spike_prob_p <= 0.05;
    sig_delta_spike_prob            = delta_spike_prob;
    sig_delta_spike_prob(~q_prob_p) = NaN;
    
    spike_prob_up       = delta_spike_prob >= 0;
    
    %%
    figure(1)
    subplot(1,length(groups),a)
    plot(delta_t(q_depth,:)',delta_spike_prob(q_depth,:)','Color',[0 0 0 .2])
    hold on
    title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
    plot(delta_t(q_depth,:),sig_delta_spike_prob(q_depth,:),'k.','MarkerSize',25)
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('Delta spike probability')
    xlabel('Delay')
    fixplot
    
    
    %% Spike rate
    spike_rate          = group_data(a).spike_rate;
    spike_rate_ctrl     = group_data(a).spike_rate_ctrl;
    delta_spike_rate    = spike_rate - spike_rate_ctrl(:);
    
    spike_rate_p        = group_data(a).spike_rate_p;
    
    q_rate_p                        = spike_rate_p <= 0.05;
    sig_delta_spike_rate            = delta_spike_rate;
    sig_delta_spike_rate(~q_rate_p) = NaN;
    
    spike_rate_up       = delta_spike_rate >= 0;
    
    %%
    figure(2)
    subplot(1,length(groups),a)
    plot(delta_t(q_depth,:)',delta_spike_rate(q_depth,:)','Color',[0 0 0 .2])
    hold on
    title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
    plot(delta_t(q_depth,:),sig_delta_spike_rate(q_depth,:),'k.','MarkerSize',25)
    xlim([-0.120 0.01])
    plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
    subplot_equal_y
    ylabel('Delta spike rateability')
    xlabel('Delay')
    fixplot
    
    
    %% 
%     first_spike_time        = group_data(a).first_spike_time;
%     delta_first_spike       = first_spike_time - first_spike_time(:,end);
%     first_spike_p           = group_data(a).first_spike_p;
%     
%     q_first_spike_p         = first_spike_p <= 0.05;
%     sig_delta_first_spike   = delta_first_spike;
%     sig_delta_first_spike(~q_first_spike_p) = NaN;
%     
%     first_spike_up        	= delta_first_spike >= 0;
%     
%     figure(3)
%     subplot(1,length(groups),a)
%     plot(delta_t,delta_first_spike(q_depth,:)','Color',[0 0 0 .2])
%     hold on
%     title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
%     plot(delta_t,sig_delta_first_spike(q_depth,:),'k.','MarkerSize',25)
%     xlim([-0.120 0.01])
%     plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
%     subplot_equal_y
%     ylabel('Delta first spike')
%     xlabel('Delay')
%     fixplot
%     
%     %% 
%     peak_spike_rate         = group_data(a).peak_spike_rate;
%     delta_peak_rate         = peak_spike_rate - peak_spike_rate(:,end);
%     
%     figure(4)
%     subplot(1,length(groups),a)
%     plot(delta_t,delta_peak_rate(q_depth,:)','Color',[0 0 0 .5])
%     hold on
%     title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
%     xlim([-0.120 0.01])
%     plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
%     subplot_equal_y
%     ylabel('Peak spike rate')
%     xlabel('Delay')
%     fixplot
%     
%     %% 
%     first_spike_jitter    	= group_data(a).first_spike_jitter;
%     delta_spike_jitter   	= first_spike_jitter - first_spike_jitter(:,end);
%     
%     figure(5)
%     subplot(1,length(groups),a)
%     plot(delta_t,delta_spike_jitter(q_depth,:)','Color',[0 0 0 .5])
%     hold on
%     title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
%     xlim([-0.120 0.01])
%     plot(xlim,[0 0],':','Color',[1 0 0],'LineWidth',2)
%     subplot_equal_y
%     ylabel('First spike jitter')
%     xlabel('Delay')
%     fixplot
%     
end



