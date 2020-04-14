
group_folders   = {'F:\Sorted_Cortical\Drive\Test'};
group_names     = {'Cortical'};

target_power    = 20;

min_depth       = 0;
max_depth       = 200;

q_reload        = 1;

tic % time script


%% function group_data = get_unit_data(group_folders) ?
if q_reload
    clear group_data
    
    groups = load_sorted_experiments(group_folders);
    
    group_units     = [];
    for a = 1:length(groups)
        group_data(a).unit_depths   = [];
        group_data(a).delta_t       = [];
        
        group_data(a).spike_rate    = [];
        group_data(a).spike_rate_p  = [];
        
        group_data(a).spike_prob    = [];
        group_data(a).spike_prob_p  = [];
        
        group_data(a).first_spike_time 	= [];
        group_data(a).first_spike_p     = [];
        
        group_data(a).first_spike_jitter	= [];
        group_data(a).peak_spike_rate       = [];
        
        group_units(a).spike_rate   = [];
        group_units(a).depths       = [];
        group_data(a).select_power  = [];
        
        for b = 1:length(groups(a).prep)
            
            opto_power = [];
            for c = 1:length(groups(a).prep(b).expt_data)
                
                ephys_data      = groups(a).prep(b).expt_data(c).ephys_data;
                
                
                
                opto_power(c)   = mode([ephys_data.conditions.LED_power])
                
            end
            
            target_power_diff   = abs(target_power - opto_power);
            target_power_ind    = find(target_power_diff == min(target_power_diff),1);
            select_power        = opto_power(target_power_ind)
            
            ephys_data          = groups(a).prep(b).expt_data(target_power_ind).ephys_data;
            
            
            drive_unit_data                     = drive_by_unit(ephys_data);
            close all
            group_data(a).spike_rate            = [group_data(a).spike_rate; drive_unit_data.mean_whisk_rates];
            group_data(a).spike_rate_p          = [group_data(a).spike_rate_p; drive_unit_data.spike_rate_p(:)];
            
            group_data(a).spike_prob            = [group_data(a).spike_prob; drive_unit_data.spike_probs];
            group_data(a).spike_prob_p          = [group_data(a).spike_prob_p; drive_unit_data.spike_prob_p(:)];
            
            group_data(a).first_spike_time      = [group_data(a).first_spike_time; drive_unit_data.first_spike_times];
            group_data(a).first_spike_p         = [group_data(a).first_spike_p; drive_unit_data.first_spike_p(:)];
%             
%             group_data(a).first_spike_jitter    = [group_data(a).first_spike_jitter; drive_unit_data.first_spike_jitter];
%             
%             group_data(a).peak_spike_rate       = [group_data(a).peak_spike_rate; drive_unit_data.peak_spike_rates];
            
            group_data(a).unit_depths           = [group_data(a).unit_depths; ephys_data.unit_depths];
            group_data(a).select_power          = [group_data(a).select_power; select_power];
        end
        
        %     [group_units(a).depths, depth_sort_inds]        = sort(group_units(a).depths);
        %     group_units(a).spikes                           = group_units(a).spikes(depth_sort_inds,:,:);
        
        [group_data(a).unit_depths, depth_sort_inds]    = sort(group_data(a).unit_depths);
        
        group_data(a).spike_rate                        = group_data(a).spike_rate(depth_sort_inds,:);
        group_data(a).spike_rate_p                      = group_data(a).spike_rate_p(depth_sort_inds,:);
        
        group_data(a).spike_prob                        = group_data(a).spike_prob(depth_sort_inds,:);
        group_data(a).spike_prob_p                     	= group_data(a).spike_prob_p(depth_sort_inds,:);
        
        group_data(a).first_spike_time              	= group_data(a).first_spike_time(depth_sort_inds,:);
        group_data(a).first_spike_p                     = group_data(a).first_spike_p(depth_sort_inds,:);
%         
%         group_data(a).first_spike_jitter                = [group_data(a).first_spike_jitter(depth_sort_inds,:)];
%         
%         group_data(a).peak_spike_rate                   = [group_data(a).peak_spike_rate(depth_sort_inds,:)];
%         
        
    end
end

figure(1)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(2)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])
figure(3)
set(gcf,'Units','normalized','Position',[.2 .4 .6 .4])

for a = 1:length(groups)
    unit_depths         = group_data(a).unit_depths;
    q_depth             = unit_depths <= max_depth & unit_depths >= min_depth;

    %% Spike probability
    
    spike_prob          = group_data(a).spike_prob;
    delta_spike_prob    = spike_prob(:,1) - spike_prob(:,2);
    spike_prob_p        = group_data(a).spike_prob_p;
    
    q_prob_p                        = spike_prob_p <= 0.05;
    sig_delta_spike_prob            = delta_spike_prob;
    sig_delta_spike_prob(~q_prob_p) = NaN;
    
    spike_prob_up       = delta_spike_prob >= 0;
    
    %%
    figure(1)
    subplot(1,length(groups),a)
    plot(ones(size(delta_spike_prob(q_depth))),delta_spike_prob(q_depth),'.','Color',[0 0 0 .5],'MarkerSize',25)
    hold on
    title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
    plot(ones(size(sig_delta_spike_prob(q_depth))),sig_delta_spike_prob(q_depth),'.','Color',[0 1 0 .5],'MarkerSize',25)
    subplot_equal_y
    ylabel('Delta spike probability')
    xlabel('Delay')
    fixplot
    
    
    %% Spike rate
    spike_rate          = group_data(a).spike_rate;
    delta_spike_rate    = spike_rate(:,1) - spike_rate(:,2);
    spike_rate_p        = group_data(a).spike_rate_p;
    
    q_rate_p                        = spike_rate_p <= 0.05;
    sig_delta_spike_rate            = delta_spike_rate;
    sig_delta_spike_rate(~q_rate_p) = NaN;
    
    spike_rate_up       = delta_spike_rate >= 0;
    
    %%
    figure(2)
    subplot(1,length(groups),a)
    plot(ones(size(delta_spike_rate(q_depth))),delta_spike_rate(q_depth),'.','Color',[0 0 0 .5],'MarkerSize',25)
    hold on
    title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])
    plot(ones(size(sig_delta_spike_rate(q_depth))),sig_delta_spike_rate(q_depth),'.','Color',[0 1 0 .5],'MarkerSize',25)
    subplot_equal_y
    ylabel('Delta spike rateability')
    xlabel('Delay')
    fixplot
    
    %% 
    first_spike_time        = group_data(a).first_spike_time;
    delta_first_spike       = first_spike_time(:,1) - first_spike_time(:,2);
    first_spike_p           = group_data(a).first_spike_p;
    
    q_first_spike_p         = first_spike_p <= 0.05;
    sig_delta_first_spike   = delta_first_spike;
    sig_delta_first_spike(~q_first_spike_p) = NaN;
    
    first_spike_up        	= delta_first_spike >= 0;
    
    figure(3)
    subplot(1,length(groups),a)
    plot(ones(size(delta_first_spike(q_depth))),delta_first_spike(q_depth),'.','Color',[0 0 0 .5],'MarkerSize',25)
    hold on
    plot(ones(size(sig_delta_first_spike(q_depth))),sig_delta_first_spike(q_depth),'.','Color',[0 1 0 .5],'MarkerSize',25)
    title([group_names{a} ', n = ' num2str(sum(q_depth)) ' units'])

    subplot_equal_y
    ylabel('Delta first spike')
    xlabel('Delay')
    fixplot
    
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



