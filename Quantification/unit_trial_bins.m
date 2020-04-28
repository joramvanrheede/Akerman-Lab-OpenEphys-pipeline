function [unit_trial_bins,bin_edges] = unit_trial_bins(spikes,bin_size, resp_win);

bin_edges = resp_win(1) : bin_size : resp_win(2);
clear a;
    for j = 1 : size(spikes,1);
        for k  = 1 : size(spikes,2);
            trial_spikes = squeeze(spikes(j,k,:));
            a(j,k,:) = histcounts(trial_spikes,bin_edges);
        end;
    end;
    unit_trial_bins = a; 
end