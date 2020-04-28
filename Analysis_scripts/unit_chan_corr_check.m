function[unit_depths,channel_depths,channel_unit_correlation] = unit_chan_corr_check(groups);
%% Checks channel correlations with unit activity. 
chan_corr = [];
unit_depths = [];
channel_depths = [];
bin_size = 0.005;
resp_win = [0 0.1];

    for z = 1 : numel(groups);
        for i = 1 : numel(groups(z).ephys_data);

            unit_trial_spikes = [];
            MUA_trial_spikes = [];

            unit_depths = [unit_depths;groups(z).ephys_data(i).unit_depths];

                for j = 1 : numel(groups(z).ephys_data(i).conditions);
    
                    unit_spikes = groups(z).ephys_data(i).conditions(j).unit_spikes; 
                    MUA_spikes = groups(z).ephys_data(i).conditions(j).multiunit_spikes;
                    binned_unit_spikes = squeeze(sum(unit_trial_bins(unit_spikes,bin_size,resp_win),3));
                    binned_MUA_spikes = squeeze(sum(unit_trial_bins(MUA_spikes,bin_size,resp_win),3));
                    unit_trial_spikes = [unit_trial_spikes binned_unit_spikes];
                    MUA_trial_spikes = [MUA_trial_spikes binned_MUA_spikes];
                end;

            clear rho pval u_d m_d;
            if ~ isempty(unit_trial_spikes);
                [rho,pval] = corr(unit_trial_spikes',MUA_trial_spikes','type','Spearman');
                [~,cc] = max(rho,[],2);
                m_d = cc*25;
            else
                rho = NaN;
                m_d = [];
            end;
           channel_depths = [channel_depths;m_d];
        end;
   end;

figure();
unit_depths = double(unit_depths);
[channel_unit_correlation(1),channel_unit_correlation(2)] = corr(unit_depths,channel_depths);
plot(unit_depths,channel_depths,'bo');
title('Channel Depth in which activity best correlates to unit activity')
xlabel('unit depths(um)');
ylabel('channel depths(um)');
end