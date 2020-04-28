% check_whisk_resp_baseline

%% File I/O settings
expt_folder         = '/Volumes/Akermanlab/Joram/In_vivo_mouse_data/2019_04_02/RBSN_2019-04-02_15-25-26_14';
q_reload            = 0; % Set to 0 once the data are loaded to save reloading time

%% Visualisation settings
bin_size            = [0.001]; % For PSTH

stim_type           = 'whisk'; % 'whisk' or 'opto'

%% Responsiveness settings
stim_win            = [0 4]; % Window for assessing frequency tracking (counting from stimulus onset)
LFP_signal_freq     = 1000; % LFP signal sample rate in Hz
smooth_win          = 1; % Smoothing window for smoothing over channels

n_layer4_chans      = 6; % Set to a sensible (conservative) nr of channels likely to be in channel 4

min_freq_for_L4     = 10; % use only frequencies above or equal to this minimum to estimate L4 bounds

%% Data loading settings

% for synchronisation data:
events_chans        = [1 1 3 2]; % order of events channels: [trial, whisk, opto, whisk_stim_nr]

% for spike detection:
spike_thresh        = 4; % threshold for spike detection

% channel order for current probe:
channel_map         = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];


%% Load and preprocess data
tic
if q_reload
    ephys_data   	= get_ephys_data(expt_folder,events_chans,channel_map, spike_thresh);
end

%% Unpack data
conditions                      = ephys_data.conditions;

cond_tracking_by_channel        = NaN(length(channel_map),length(conditions));
peak_tracking                   = NaN(length(conditions),1);
upper_L4_bound                  = NaN(length(conditions),1);
lower_L4_bound                  = NaN(length(conditions),1);
L4_center                       = NaN(length(conditions),1);
top_channels                    = NaN(length(conditions),n_layer4_chans);
stim_freqs                      = NaN(length(conditions),1);

close all
figure(1)
set(gcf,'Units','Normalized','Position',[.1 .1 .8 .8])
figure(2)
set(gcf,'Units','Normalized','Position',[.1 .1 .8 .8])
figure(3)
set(gcf,'Units','Normalized','Position',[.1 .1 .8 .8])
for a = 1:length(conditions)
    
    spikes              = conditions(a).spikes;
    LFPs                = conditions(a).LFP_trace;

    switch stim_type
        case 'whisk'
            stim_onset      = conditions(a).whisk_start;
            stim_freq       = conditions(a).whisk_freq;
        case 'opto'
            stim_onset      = conditions(a).opto_onset;
            stim_freq       = conditions(a).opto_freq;
    end
    
    stim_freqs(a)   = stim_freq;
    
    spikes      	= spikes - stim_onset; % Set whisk onset time to 0
    LFP_times       = ([1:size(LFPs,3)]/LFP_signal_freq) - stim_onset; % Generate timestamps for LFP signal
    
    % Determine the time window of the stimulus to set limits of where to report tracking
    q_stim_win   	= LFP_times > stim_win(1) * LFP_times <= stim_win(2);
    
    stim_LFPs       = LFPs(:,:,q_stim_win);
    stim_LFP_times  = LFP_times(q_stim_win);
    
    % LFP_power_spectra   
    
    % Get the LFP band power in this band for each channel and each trial
    tracking_band_power     = LFP_band_power(stim_LFPs, stim_freq, LFP_signal_freq);
    
    % Average over trials to get mean for each channel
    tracking_by_channel     = mean(tracking_band_power,2);
    
    % Smooth
    tracking_by_channel     = smooth(tracking_by_channel,smooth_win);
    
    % Normalise
    tracking_by_channel     = tracking_by_channel / max(tracking_by_channel);
    
    % Find the best tracking channel
    peak_tracking(a)        = find(tracking_by_channel == max(tracking_by_channel));
    
    % Find the k best tracking channels
    [top_k, k_inds]         = maxk(tracking_by_channel,6);
    
    % Save the top channels in order
    top_channels(a,:)       = sort(k_inds);
    
    % Estimate 
    chans_above_half_max    = tracking_by_channel > 0.5;
    
    % find the longest stretch of channels where the tracking is above 50%
    % of the maximum tracking, as an estimate of which channels belong to
    % Layer 4 --> it gives an upper and lower bound as well as a centre
    is_above        = 0;
    longest_stretch = 0;
    longest_stretch_start   = 0;
    longest_stretch_end     = 0;
    for b = 1:length(chans_above_half_max)
        this_val = chans_above_half_max(b);
        
        if is_above == 0 && this_val == 1
            % we were not above half max, but now we are --> start counting
            start_ind   = b;
            is_above    = 1;
        end
        
        if is_above == 1 && this_val == 0
            % we were above half max, but not anymore --> calculate length
            % of this stretch and compare it to previous longest stretch
            % to see if this is now the longest. Set is_above back to 0.
            is_above = 0;
            end_ind     = b-1;
            this_stretch  = (end_ind - start_ind) + 1;
            if this_stretch > longest_stretch
                longest_stretch         = this_stretch;
                longest_stretch_start   = start_ind;
                longest_stretch_end     = end_ind;
            end
        end
        
        if b == length(chans_above_half_max) && this_val == 1
            % We have reached the end of the channels and we are still
            % tracking at above half max. Set end to this channel and see 
            % if this was the longest stretch.
            end_ind     = b;
            this_stretch  = (end_ind - start_ind) + 1;
            if this_stretch > longest_stretch
                longest_stretch         = this_stretch;
                longest_stretch_start   = start_ind;
                longest_stretch_end     = end_ind;
            end
        end
    end
    
    % Save the longest stretch boundaries and the center
    upper_L4_bound(a)   = longest_stretch_start;
    lower_L4_bound(a)   = longest_stretch_end;
    L4_center(a)        = mean([longest_stretch_start longest_stretch_end]);
    
    % Add to a variable that stores tracking by channel for all conditions
    cond_tracking_by_channel(:,a)   = tracking_by_channel;
    
    
    % Plot LFP traces using plot_LFP_traces function
    figure(1)
    subplot(length(conditions),1,a)
    plot_LFP_traces(LFPs,1,(1:size(LFPs,3))/1000 - stim_onset);
    x_limits = xlim;
    fixplot
    subplot_equal_y
    title(['LFP trace ' num2str(round(stim_freq)) 'Hz'])
    
    % Plot spiking response using psth function
    figure(2)
    subplot(length(conditions),1,a)
    [psth_handle, psth_counts, psth_bins] = psth(spikes,bin_size);
    fixplot
    subplot_equal_y
    title(['Post-stimulus time histogram ' num2str(round(stim_freq)) 'Hz'])
    
    % Plot spiking response using psth function
    figure(3)
    subplot(length(conditions),1,a)
    spike_density_plot(spikes,1,psth_bins);
    title(['Spike density plot ' num2str(round(stim_freq)) 'Hz'])
    
end


%% Determine putative L4 boundaries from the relevant tracking data

% create boolean to determine which frequencies to include in layer boundary
% analysis; and use to select only data for those frequencies from all 
% relevant variables:
q_use_freq      = stim_freqs >= min_freq_for_L4;
stim_freqs    	= stim_freqs(q_use_freq);

upper_L4_bound  = upper_L4_bound(q_use_freq);
lower_L4_bound  = lower_L4_bound(q_use_freq);
L4_center       = L4_center(q_use_freq);

top_channels    = top_channels(q_use_freq,:);
peak_tracking   = peak_tracking(q_use_freq);

cond_tracking_by_channel = cond_tracking_by_channel(:,q_use_freq);

% Determine the estimated upper and lower bounds
max_upper_bound = max(upper_L4_bound);
min_upper_bound = min(upper_L4_bound);

max_lower_bound = max(lower_L4_bound);
min_lower_bound = min(lower_L4_bound);

tick_channels   = [1 32 max_upper_bound min_upper_bound max_lower_bound min_lower_bound];
tick_channels   = unique(tick_channels);

%% Plot frequency response across channels and indicate layer boundary estimates

figure
set(gcf,'Units','Normalized','Position',[.3 .2 .4 .6])
n_channels  = length(channel_map);
plot(cond_tracking_by_channel,1:n_channels,'LineWidth',4)
hold on
axis ij

legendlabels    = {};
for a = 1:length(stim_freqs)
    legendlabels{a} = [num2str(round(stim_freqs(a))) 'Hz'];
end

title('LFP stim tracking by channel')
xlabel('Normalised power at stim frequency band (+/- 10%)')
ylabel('Channel number')
x_limits    = xlim;

line([0.5 0.5],ylim,'LineStyle',':','LineWidth',4,'Color',[.6 .6 .6])
xlim(x_limits)
shaded_y_region([min_upper_bound max_lower_bound],[1 1 0],0.25)
xlim(x_limits)
shaded_y_region([max_upper_bound min_lower_bound],[1 0 0],0.25)
xlim(x_limits)
ylim([0 n_channels+1])
legendlabels{length(stim_freqs)+1} = 'Half max (L4 cutoff)';
legendlabels{length(stim_freqs)+2} = 'Inclusive L4 bounds';
legendlabels{length(stim_freqs)+3} = 'Restrictive L4 bounds';
legend_handle = legend(legendlabels,'Location','SouthEast','AutoUpdate','off');

set(gca,'YTick',tick_channels)
fixplot

% Draw the recording electrode, for fun
electrode_shape = polyshape([0.03 0.01 0.01 0.05 0.05 0.03],[n_channels+1 n_channels 0 0 n_channels n_channels+1]);
plot(electrode_shape,'FaceColor',[.6 .6 .6],'EdgeAlpha',0,'FaceAlpha',1)
plot([0.03],[1:n_channels],'k.','MarkerSize',20)

hold off

%% Report back to user

disp(' ')
disp(['Stimulation frequencies used for boundary determination: ' num2str(round(stim_freqs')) ])
disp(['Respective peak tracking channels: ' num2str(peak_tracking')])
disp(' ')
disp(['Top ' num2str(n_layer4_chans) ' channels, by frequency (1 row per frequency):' ])
disp(top_channels)
disp(' ')
disp('Estimated Layer 4 boundaries:')
disp([' '])
disp('Most inclusive thresholded criterion for L4 (bounds in any freq):')
disp(['Upper channel: ' num2str(min_upper_bound)])
disp(['Lower channel: ' num2str(max_lower_bound)])
disp(['Estimated L4 centre location (channel): ' num2str(mean([min_upper_bound max_lower_bound]))])
disp([' '])
disp(['Most restrictive thresholded criterion for L4 (bounds in common between all freqs):'])
disp(['Upper channel: ' num2str(max_upper_bound)])
disp(['Lower channel: ' num2str(min_lower_bound)])
disp(['Estimated L4 centre location (channel): ' num2str(mean([max_upper_bound min_lower_bound]))])
disp([' '])

toc
