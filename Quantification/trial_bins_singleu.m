function [unit_trial_bins,bin_edges] = trial_bins_singleu(spikes,bin_size, resp_win);

bin_edges = resp_win(1) : bin_size : resp_win(2);
clear a;
    for j = 1 : size(spikes,1);
            trial_spikes = squeeze(spikes(j,:));
            a(j,:) = histcounts(trial_spikes,bin_edges);
    end;
    unit_trial_bins = a; 
end