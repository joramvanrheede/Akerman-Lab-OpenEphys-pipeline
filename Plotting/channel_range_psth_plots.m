% channel_range_LFP_plots

%% Parameters
cond_nr         = 5; % which condition to plot this for?
expt_nr         = [2];
day_nr          = [2];

select_samples  = [500:1500];
kill_samples    = [];%[501 511 512];%[501:511]; % to remove artefacts

select_channels = [1:32];
kill_channels   = [];%[23 27]%[7 23 26 27];

chan_offset     = [.2];

% psth_smoothing  = 9;

%% execution starts here


[kill_chans,kill_inds,temp_inds]    = intersect(select_channels, kill_channels);
select_channels(kill_inds)  = [];

[kill_samps,kill_inds,temp_inds]    = intersect(select_samples, kill_samples);
% select_samples(kill_inds)   = NaN;


whisk_win_rates         = squeeze(sdata(day_nr).expt(expt_nr).whisk_profile(cond_nr,select_channels,select_samples));
whisk_win_rates(:,kill_inds) = NaN;

trace_offsets           = repmat([0:size(whisk_win_rates,1)-1]',1,size(whisk_win_rates,2));
whisk_win_rates         = whisk_win_rates + flipud(trace_offsets) * chan_offset;

figure
plot(whisk_win_rates','LineWidth',3)

% Aesthetics
axis tight
axis off
ylims = ylim;
set(gca,'YLim',[ylims(1)-chan_offset ylims(2)+chan_offset])

%% Give informative title
condition_mat   = sdata(day_nr).expt(expt_nr).condition_mat;
cond_row        = condition_mat(cond_nr,:);

cond_row_str    = 'Cond';
for a = 1:length(cond_row)
    cond_row_str = [cond_row_str ' ' num2str(cond_row(a))];
end

title(cond_row_str)
scf

