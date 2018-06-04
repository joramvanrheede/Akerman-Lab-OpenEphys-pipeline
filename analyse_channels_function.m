function channels = analyse_channels_function(channels,analysisparams)
% function channels = analyse_drive(channels,analysisparams)
% prepare spike count analysis of 'channels' structure using the parameters
% in 'analysisparams' struct
% 'analysisparams' has the folllowing fields:
% analysisparams.whiskwinedges
% analysisparams.LEDwin
% analysisparams.whiskwin
% analysisparams.LED_sust_win
% analysisparams.LEDwinedges
% analysisparams.samplerate - sample rate of data; is used to set the bin size to determine whisker spike rate for maximal resolution
% analysisparams.profile_smoothing - smoothing window for spike profile in milliseconds
% 
% returns 'channels' structure with additional fields 
% additions:
% PSTHs
% Spike counts in target windows
% Spike rates in target windows
% Spike rates relative to spontaneous activity on that channel
% Profile of whisker response over time
% This also needs to be done for LED - in particular, detection of peak
% response in different channels, to see temporal order of activation
% between channels if any
% --> should this just be a complete profile over time?

%% unpack parameters struct
whiskwinedges       = analysisparams.whiskwinedges;
LEDwin              = analysisparams.LEDwin;
whiskwin            = analysisparams.whiskwin;
LED_sust_win        = analysisparams.LED_sust_win;
LEDwinedges         = analysisparams.LEDwinedges;
samplerate          = analysisparams.samplerate; % sample rate of data; is used to set the bin size to determine whisker spike rate for maximal resolution
profile_smoothing   = analysisparams.profile_smoothing; % smoothing window for spike profile in milliseconds

%% Running code starts here

cond_nr             = length(channels(1).conditions);

conditions  = [];
for a = 1:length(channels(1).conditions)
    conditions(a,:) = channels(1).conditions(a).timings;
end

LED_conds       = conditions(:,1);
LED_wins        = repmat(LED_conds,1,2) + repmat(LEDwin,length(LED_conds),1);

[channels(1:end).LED_wins]   = deal(LED_wins);

for a = 1:length(channels)
    
    for b = 1:cond_nr
        
        n_episodes          = length(channels(a).conditions(b).episodes);
        
        condspikes          = [];
        rasterspikes        = [];
        for c = 1:n_episodes
            thesespikes         = channels(a).conditions(b).episodes(c).spikes;
            condspikes          = [condspikes; thesespikes];
            
            rasterspikes(1:length(thesespikes),c) = thesespikes;
            
        end
        
        rasterspikes(rasterspikes == 0) = NaN;
        
        whisk_spikes        = condspikes - conditions(b,2);
        
        channels(a).conditions(b).allspikes     = condspikes;
        channels(a).conditions(b).rasterspikes  = rasterspikes;
        
        
        %% Using the aligned spike times for all episodes, make a profile 
        % of the whisker response, in order to determine the maximum rate
        % of fire and the time of that peak ROF
        
        % create a vector of bins at the samplerate and count spikes
        whisk_profile_time_vect     = min(whiskwinedges):1/samplerate:max(whiskwinedges);
        whisk_profile_counts        = histc(whisk_spikes,whisk_profile_time_vect);
        
        % determine the number of samples for the smoothing kernel from the
        % smoothing window (s) and sample rate (Hz)
        n_smooth_samples            = profile_smoothing * samplerate;
        
        % for convolution, a kernel with an odd number of elements is
        % needed; add extra element if it is even
        if mod(n_smooth_samples,2) == 0
            n_smooth_samples        = n_smooth_samples + 1; 
        end
        
        kernel_vect                 = [-3:(6/n_smooth_samples):3]; % let vector run from -3(SDs) to +3(SDs)
        gauss_kernel                = normpdf(kernel_vect,0,1); % create gaussian using the normal probability density function
        
        % ensure that the kernel has a total area under the curve of 1
        gauss_kernel                = gauss_kernel / sum(gauss_kernel);
        
        % use kernel for convolution with the spike coutns binned at sample rate
        whisk_profile               = conv(whisk_profile_counts,gauss_kernel,'same');
        
        whisk_peak_test_vect        = whisk_profile(whisk_profile_time_vect > whiskwin(1) & whisk_profile_time_vect < whiskwin(2));
        whisk_peak_time_vect        = whisk_profile_time_vect(whisk_profile_time_vect > whiskwin(1) & whisk_profile_time_vect < whiskwin(2));
        
        whisk_peak_height           = max(whisk_peak_test_vect); % how to convert this back to spike rate equivalent?
        whisk_peak_ind              = find(whisk_peak_test_vect == whisk_peak_height,1,'first');
        whisk_peak_time             = whisk_peak_time_vect(whisk_peak_ind);
        
        whisk_peak_rate             = (whisk_peak_height * samplerate) / n_episodes;
        % now determine peak rate
        
        %%
        
        whiskwinsize                = (whiskwin(2:end) - whiskwin(1:end-1));
        
        whiskcount                  = histc(whisk_spikes, whiskwin);
        whiskcount                  = whiskcount(1:end-1);
        
        channels(a).conditions(b).whisk_spike_count   	= whiskcount;
        channels(a).conditions(b).whisk_spike_rate      = (whiskcount ./ whiskwinsize') / n_episodes;
        channels(a).conditions(b).whisk_spike_rel   	= channels(a).conditions(b).whisk_spike_rate / mean(channels(a).spontspikerate);
        channels(a).conditions(b).whisk_profile         = whisk_profile;
        channels(a).conditions(b).whisk_profile_time_vect   = whisk_profile_time_vect;
        channels(a).conditions(b).whisk_peak_time       = whisk_peak_time;
        channels(a).conditions(b).whisk_peak_height     = whisk_peak_height;
        channels(a).conditions(b).whisk_peak_rate       = whisk_peak_rate;
        
        %%
        
        whiskwinsizes       = (whiskwinedges(2:end) - whiskwinedges(1:end-1));
        
        whiskhist          	= histc(whisk_spikes, whiskwinedges);
        whiskhist        	= whiskhist(1:end-1);
        
        channels(a).conditions(b).whisk_win_counts      = whiskhist;
        channels(a).conditions(b).whisk_win_rates       = (whiskhist ./ whiskwinsizes') / n_episodes;
        channels(a).conditions(b).whisk_rel_rates      	= channels(a).conditions(b).whisk_win_rates / mean(channels(a).spontspikerate);
        
        
        %% LED ON
        
        LED_spikes        	= condspikes - conditions(b,1);
        
        LEDwinsize        	= (LEDwin(2:end) - LEDwin(1:end-1));
        
        LEDcount          	= histc(LED_spikes, LEDwin);
        LEDcount          	= LEDcount(1:end-1);
        
        channels(a).conditions(b).LED_spike_count      	= LEDcount;
        channels(a).conditions(b).LED_spike_rate      	= (LEDcount ./ LEDwinsize') / n_episodes;
        channels(a).conditions(b).LED_spike_rel         = channels(a).conditions(b).LED_spike_rate / mean(channels(a).spontspikerate);
        
        %% LED OFF
        
        LED_OFF_spikes      = LED_spikes - conditions(b,3);
        LEDOFFcount        	= histc(LED_OFF_spikes, LEDwin);
        LEDOFFcount        	= LEDOFFcount(1:end-1);
        
        channels(a).conditions(b).LED_OFF_spike_count 	= LEDOFFcount;
        channels(a).conditions(b).LED_OFF_spike_rate  	= (LEDOFFcount ./ LEDwinsize') / n_episodes;
        channels(a).conditions(b).LED_OFF_spike_rel   	= channels(a).conditions(b).LED_OFF_spike_rate / mean(channels(a).spontspikerate);
        
        
        %% Sustained LED
        
        LED_sust_count        	= histc(LED_spikes, LED_sust_win);
        LED_sust_count        	= LED_sust_count(1:end-1);
        LED_sust_winsize        = (LED_sust_win(2:end) - LED_sust_win(1:end-1));
        
        channels(a).conditions(b).LED_sust_spike_count 	= LED_sust_count;
        channels(a).conditions(b).LED_sust_spike_rate  	= (LED_sust_count ./ LED_sust_winsize') / n_episodes;
        channels(a).conditions(b).LED_sust_spike_rel   	= channels(a).conditions(b).LED_sust_spike_rate / mean(channels(a).spontspikerate);
        
        
        %% Activity in LED windows (when LED ON as well as OFF, for post timing experiment comparison)
        
        LEDwinsizes        	= (LEDwinedges(2:end) - LEDwinedges(1:end-1));
        LEDhist             = histc(LED_spikes, LEDwinedges);
        LEDhist             = LEDhist(1:end-1);
        
        
        channels(a).conditions(b).LED_win_counts        = LEDhist;
        channels(a).conditions(b).LED_win_rates         = (LEDhist ./ LEDwinsizes') / n_episodes;
        channels(a).conditions(b).LED_rel_rates         = channels(a).conditions(b).LED_win_rates / mean(channels(a).spontspikerate);
              
    end
end




