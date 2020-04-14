function [ephys_data] = extract_ephys_data(datafolder,datafilenr,parameters)

data_prefix             = parameters.data_prefix;           % = '100';
spike_thresh            = parameters.spike_thresh;          % spike threshold in SDs
LED_conditions_res      = parameters.LED_conditions_res;    % resolution in ms for automatically extracting conditions from LED delays
whisk_conditions_res    = parameters.whisk_conditions_res;  % resolution in ms for automatically extracting conditions from whisker delays

trial_input_nr          = parameters.trial_channel;         % Which input channel has the trial TTL
stim_input_nr           = parameters.whisk_channel;         % Which input channel has the stim / whisk TTL
opto_input_nr        	= parameters.LED_channel;           % Which input channel has the LED TTL
switch_input_nr         = parameters.stim_switch_channel;  	% Which input channel switches between stimulators?

expt_type               = parameters.experiment_type;       % what type of experiment is this?

get_channels            = parameters.get_channels;        	% which channels to get in which order?

get_LFP                 = parameters.get_LFP;               % get LFP data?

trials_from_whisk       = parameters.trials_from_whisk;     % discard trial information from ADC channels and determine trials based on whisker instead?
whisk_buffer            = parameters.whisk_buffer;          % if using whisk stim to divide recording into trials (above), trials start whisk_buffer (in seconds) before the whisker stim onset, and end 2*whisk buffer after whisker stim ONSET

do_CAR                  = parameters.do_CAR;                % If true, do common average referencing

save_sync_chans         = parameters.save_sync_chans;       % save the synchronisation channel outputs
sync_chans_res          = parameters.sync_chans_res;        % resolution for the synch channel outputs

baseline_moving_window  = parameters.baseline_moving_window;

%% Hard-coded parameters

samplefreq              = 30000;
spike_filt_band         = [500 5000];
LFP_filt_band           = [1 300];
max_opto_freq           = 100;

%% some file I/O pre-work

n_channels              = length(get_channels);

data_contents           = dir(datafolder); % find what is in the data folder
filefolders             = data_contents([data_contents.isdir]); % get only directories
filefolders             = {filefolders.name}'; % get the names of the folders
filefolders             = filefolders(3:end); % remove '.' and '..' folders

pattern                 = '._\d+$'; % regular expression to search for file numbers below

[startinds, endinds]    = regexp(filefolders, pattern, 'start','end'); % find number pattern in the remaining filefolders
filenumbers             = [];

% generate a vector of file numbers from the files in the data folder
for i = 1:length(startinds)
    filenumbers(i)  = str2num(filefolders{i}(startinds{i}+2:endinds{i}));
end

fileind                 = find(filenumbers == datafilenr); % find index of target data folder
filefolder              = filefolders{fileind}; % this is the folder we're after

%% Start collecting events data

if switch_input_nr ~= 0
    trial_threshold     = 0.25; % For trial, signal goes up to 2.5V or less
    stim_threshold      = 2.7; % For whisking (always during trial), signal goes up to 2.75 - 5V
    opto_threshold    	= 0.0425; % LED is no longer TTL, voltage varies with power (with 5V representing max); 0.05V thresh will detect events above ~1% max power
    switch_threshold    = 2.5; % normal TTL logic - 0 to 5V
else
    trial_threshold     = 0.25; % normal TTL logic - 0 to 5V
    stim_threshold      = 2.5; % normal TTL logic - 0 to 5V
    opto_threshold    	= 0.0425; % LED is no longer TTL, voltage varies with power (with 5V representing max); 0.05V thresh will detect events above ~1% max power
    switch_threshold    = 2.5; % normal TTL logic - 0 to 5V
end



adc_channel_nrs        	= [trial_input_nr stim_input_nr opto_input_nr switch_input_nr];
adc_channel_thresholds 	= [trial_threshold stim_threshold opto_threshold switch_threshold];

for a = 1:4 % loop through the analog input channels
    
    if a == 1 || adc_channel_nrs(a) ~= adc_channel_nrs(a-1) % Don't reload data if trace is already loaded
        disp(['Loading ADC input channel ' num2str(adc_channel_nrs(a))])
        disp(['File ' datafolder filesep filefolder filesep data_prefix '_ADC' num2str(adc_channel_nrs(a)) '.continuous'])
        [thisTTL timestamps info] = load_open_ephys_data([datafolder filesep filefolder filesep data_prefix '_ADC' num2str(adc_channel_nrs(a)) '.continuous']);
        
        thisTTL         = thisTTL - min(smooth(thisTTL(:),101));
        
        %% Special case for fixing baseline of LED input channel:
        if a == 3
            
            % Moving minimum based on LED conditions res.
            % Operates on resampled 'thisTTL' with 1 sample per 10ms; if
            % vector left in original length this operation takes too long.
            moving_TTL_min  = movmin(thisTTL(1:30:end),baseline_moving_window);
            
            % smooth transitions with window
            smooth_TTL_min  = smooth(moving_TTL_min,7);
            smooth_TTL_min  = [smooth_TTL_min(:); smooth_TTL_min(end)]; % make sure this vector when resampled will be longer than thisTTL
            
            % bring back to full 30kHz sample rate
            resamp_smooth_TTL_min   = resample(smooth_TTL_min,30,1);
            
            % Correct original TTL
            corr_TTL        = thisTTL - resamp_smooth_TTL_min(1:length(thisTTL));
            
            baseline_wobble = range(smooth_TTL_min);
            
            if baseline_wobble < 0.025 % = up to 05%
                disp(['Baseline stable (' num2str(baseline_wobble/0.05) '% baseline wobble using ' num2str(baseline_moving_window) 'ms moving minimum)'])
            elseif baseline_wobble < 0.05 % = up to 1%
                beep
                warning(['Minor baseline instability (' num2str(baseline_wobble/0.05) '% baseline wobble using ' num2str(baseline_moving_window) 'ms moving minimum)'])
                disp(['Adjusting minimum threshold for event detection to 1.5% of range'])
                adc_channel_thresholds(a) = 0.075;
            elseif baseline_wobble < 0.1 % = up to 2%
                beep
                warning(['Moderate baseline instability (' num2str(baseline_wobble/0.05) '% baseline wobble using ' num2str(baseline_moving_window) 'ms moving minimum)'])
                disp(['Adjusting minimum threshold for event detection to 3 % of range'])
                adc_channel_thresholds(a) = 0.15;
            elseif baseline_wobble < 0.25 % = up to 5%
                beep
                warning(['High baseline instability (' num2str(baseline_wobble/0.05) '% baseline wobble using ' num2str(baseline_moving_window) 'ms moving minimum)'])
                disp(['Adjusting minimum threshold for event detection to 7.5% of range'])
                adc_channel_thresholds(a) = 0.375;
            elseif baseline_wobble > 0.25
                beep
                warning(['Severe baseline instability (' num2str(baseline_wobble/0.05) '% baseline wobble using ' num2str(baseline_moving_window) 'ms moving minimum)'])
                disp(['Make sure to check baseline moving minimum window'])
                disp(['Adjusting minimum threshold for event detection to 10% of range'])
                adc_channel_thresholds(a) = 0.5;
            end
            
            %             % Uncomment For debugging
            %             figure
            %             plot(thisTTL(1:100:end),'k-','LineWidth',2)
            %             hold on
            %             plot(corr_TTL(1:100:end),'b-')
            %             plot(resamp_smooth_TTL_min(1:100:end),'c-')
            %             plot([0 length(thisTTL)/100], [adc_channel_thresholds(a) adc_channel_thresholds(a)],'r-')
            %             keyboard
            %             % end debugging code
            
            % Apply moving minimum correction
            thisTTL         = corr_TTL;
        end
        
        starttime       = min(timestamps);              % find start time
        endtime         = max(timestamps);              % find end time
        timestamps      = (1:length(thisTTL)) / 30000;  % manually create new timestamps at 30kHz, openephys sometimes suffers from timestamp wobble even though data acquisition is spot on
        timestamps      = timestamps + starttime;       % add start time to the newly created set of timestamps
    end
    
    thisTTL_bool   	= thisTTL > adc_channel_thresholds(a); % find where the TTL signal is 'high'
    
    start_inds      = find(diff(thisTTL_bool) > 0.5); % find instances where the TTL goes from low to high
    end_inds        = find(diff(thisTTL_bool) < -0.5); % find instances where the TTL goes from high to low
    
    if length(start_inds)>length(end_inds)
        end_inds = [end_inds; length(thisTTL)];
    end
    
    start_times 	= timestamps(start_inds); % find the timestamps of start events
    end_times    	= timestamps(end_inds); % find the timestamps of end events
    
    if ~isempty(start_times) & ~isempty(end_times)  % Some channels may not have events (e.g. stim switch channel if only 1 stimulator used)
        end_times(end_times < start_times(1))       = []; % discard potential initial end without start
        start_times(start_times > end_times(end))   = []; % discard potential final start without end
    end
    
    switch a % this determines what the start and end timestamps should be assigned to: trial/trial, LED/opto stim or stim/whisk stim.
        case 1
            trial_starts    = start_times(:);
            trial_ends      = end_times(:);
        case 2
            stim_starts 	= start_times(:);
            stim_ends     	= end_times(:);
            
            % Determine stimulus amplitude from signal
            stim_amps       = NaN(size(start_inds));
            
            for i = 1:length(start_inds)
                stim_segment    = thisTTL(start_inds(i):end_inds(i));
                if isempty(stim_segment)
                    stim_segment = NaN;
                end
                stim_amps(i)   	= ((max(stim_segment) - 2.5) / 2.5) * 100; % Stimulus amplitude in % of max
            end
        case 3
            opto_starts      = start_times(:);
            opto_ends        = end_times(:);
            
            opto_powers      = NaN(size(start_inds));
            
            for i = 1:length(start_inds)
                stim_segment    = thisTTL(start_inds(i):end_inds(i));
                if isempty(stim_segment)
                    stim_segment = NaN;
                end
                opto_powers(i) 	= max(stim_segment) / 5 * 100; % Stimulus amplitude in % of max
            end
        case 4
            switch_up       = start_times(:);
            switch_down     = end_times(:);
    end
    if save_sync_chans
        resample_freq               = round(samplefreq / sync_chans_res);
        sync_chan_traces(1:length(thisTTL(1:resample_freq:end)),a)  = thisTTL(1:resample_freq:end);
        if a == 1
            sync_chan_timestamps 	= timestamps(1:resample_freq:end);
        end
    end
end


%% A lot of cleanup and repair from here

if trials_from_whisk
    % We are setting trial starts and ends based on the whisker stimulus;
    % discard previous trial data
    trial_starts    = [];
    trial_counter   = 1;
    for a = 1:length(stim_starts)
        this_stim_start = stim_starts(a);
        if a == 1
            % first stimulus, so this will set the first trial start (= this_stim_start - whisk_buffer)
            trial_starts(trial_counter) = this_stim_start - whisk_buffer;
            trial_counter   = trial_counter + 1; % keep track of which trial we are on
            continue
        elseif (stim_starts(a) - whisk_buffer) <= stim_starts(a-1)
            % the interval between this stim start and the previous one is
            % too small for them to have happened in different trials;
            % don't increment trial number, and investigate next stim start
            % time
            continue
        elseif (stim_starts(a) - whisk_buffer) > stim_starts(a-1)
            % interval between this stim start and the previous one is
            % large and so this stimulus is happening in a new trial;
            % set new trial start, and increment trial counter
            trial_starts(trial_counter) = this_stim_start - whisk_buffer;
            trial_counter = trial_counter + 1;
        end
    end
    trial_ends = trial_starts + 2 * whisk_buffer;
end

% determine median trial length
trial_times     = trial_ends - trial_starts;
trial_length    = round(median(trial_times),1);

% all trials should have the same length; trials with anomalous
% length are likely arduino startup floating voltage artefacts; get rid
% of anomalies
qtrial          = round(trial_times,1) == trial_length;

trial_starts    = trial_starts(qtrial);
trial_ends      = trial_ends(qtrial);

total_length 	= round(median(diff(trial_starts)),3);
trial_gap       = median(trial_starts(2:end)-trial_ends(1:end-1));

%% At the end, simply chuck trials with no events in them?

if isempty(stim_starts)
    stim_starts = [0 0.02];
    stim_ends   = [0.01 00.03]; % set some fake whisk stimuli outside of the trials
    stim_amps   = [1 1];
end


%%

allwhisks                   = stim_starts; % ?
allwhisk_ends               = stim_ends;

%% dealing with bursts of whisker stimuli

first_stim_inds             = find(diff(allwhisks) > trial_gap)+1;

if ~isempty(first_stim_inds)
    % Find which whisking onsets are the first of a trial, and which onsets
    % are the last of a trial
    allwhisk_firstvect          = allwhisks(first_stim_inds);
    allwhisk_lastvect           = allwhisks(first_stim_inds-1);% Stimulus amplitude
    
    whisk_ends               	= allwhisk_ends(first_stim_inds);
    whisk_ends                	= [allwhisk_ends(1); whisk_ends(:)];
    
    stim_amps                   = stim_amps(first_stim_inds);
    stim_amps                   = [stim_amps(1); stim_amps(:)];
else
    allwhisk_firstvect          = [];
    allwhisk_lastvect           = [];
    
    whisk_ends                  = allwhisk_ends;
end

whisk_starts            	= [allwhisks(1); allwhisk_firstvect(:)];
whisk_lasts              	= [allwhisk_lastvect(:); allwhisks(end)];

% Stimulus length and repeat frequency
whisk_lengths              	= whisk_ends - whisk_starts;
whisk_freqs              	= NaN(size(whisk_starts));

for a = 1:length(whisk_starts)
    this_whisk_start    = whisk_starts(a);
    this_whisk_end      = whisk_lasts(a);
    q_whisks            = allwhisks > this_whisk_start & allwhisks < this_whisk_end;
    
    this_whisk_freq   	= mean(round(1./diff(allwhisks(q_whisks))));
    if isempty(this_whisk_freq)
        this_whisk_freq = 99;
    elseif isnan(this_whisk_freq)
        this_whisk_freq = 99;
    end
    whisk_freqs(a)      = this_whisk_freq;
end

allwhisks             	= whisk_starts;

%% Dealing with bursts of opto stimuli

if isempty(opto_starts)
    opto_starts = [0 0.02];
    opto_ends   = [0.01 00.03]; % set some fake opto stimuli outside of the trials
    opto_powers = [1 1];
end

% Find instances where the difference between a previous and next stimulus
% is more than the gap between trials; this should be where one trial ends and the
% next one begins
first_opto_inds             = find(diff(opto_starts) > trial_gap)+1;
last_opto_inds              = find(diff(opto_starts) > trial_gap);

% make sure to include the first onset and the last offset, which will not
% be captured by the diff criterion (no gap before the start, no gap after
% the end)
first_opto_inds             = [1; first_opto_inds];
last_opto_inds              = [last_opto_inds; length(opto_starts)];

if ~isempty(first_opto_inds)
    opto_firsts                 = opto_starts(first_opto_inds);
    opto_lasts                  = opto_ends(last_opto_inds);
    
    opto_first_amps             = opto_powers(first_opto_inds);
    opto_last_amps              = opto_powers(last_opto_inds);
    opto_amps                   = max(opto_first_amps,opto_last_amps);
    
else
    opto_firsts                 = [];
    opto_lasts                  = [];
    
    opto_burst_ends             = opto_ends;
end

opto_freqs              	= NaN(size(opto_starts));

for a = 1:length(opto_firsts)
    this_opto_first     = opto_firsts(a);
    this_opto_last   	= opto_lasts(a);
    
    q_opto_burst        = opto_starts > this_opto_first & opto_starts < this_opto_last;
    
    this_opto_freq   	= mean(round(1./diff(opto_starts(q_opto_burst))));
    if isempty(this_opto_freq)
        this_opto_freq  = 99;
    elseif isnan(this_opto_freq)
        this_opto_freq  = 99;
    elseif this_opto_freq > max_opto_freq
        this_opto_freq  = 99;
    end
    
    opto_freqs(a)      = this_opto_freq;
    
end

%% Match events to trials
ntrials                 = length(trial_starts);

whisk_stim_onsets       = NaN(size(trial_starts));
whisk_stim_lengths      = NaN(size(trial_starts));
whisk_stim_freqs        = NaN(size(trial_starts));
whisk_stim_relay        = NaN(size(trial_starts));
whisk_stim_amplitudes   = NaN(size(trial_starts));

opto_onsets             = NaN(size(trial_starts));
opto_offsets            = NaN(size(trial_starts));
opto_current_levels     = NaN(size(trial_starts));
opto_freq               = NaN(size(trial_starts));

for a = 1:ntrials
    this_trial_start    = trial_starts(a);
    this_trial_end      = trial_ends(a);
    
    % see whether there was a whisker stimulus
    select_whisk_start 	= whisk_starts >= this_trial_start & whisk_starts <= this_trial_end;
    
    if sum(select_whisk_start) > 0
        if sum(select_whisk_start) > 1
            beep
            warning('Multiple whisker stimulus values found for this trial')
            
            first_whisk_start_ind                       = find(select_whisk_start,1);
            select_whisk_start                          = false(size(select_whisk_start));
            select_whisk_start(first_whisk_start_ind)   = true;
        end
        
        whisk_stim_onsets(a)            = stim_starts(select_whisk_start);
        whisk_stim_lengths(a)       = whisk_lengths(select_whisk_start);
        whisk_stim_freqs(a)         = whisk_freqs(select_whisk_start);
        whisk_stim_amplitudes(a)    = stim_amps(select_whisk_start);
        
        % Determine which stimulator is being used (relay up = stim 2, relay down = stim 1)
        stim_start_mat_temp         = repmat(whisk_stim_onsets(a),size(switch_up));
        is_switch_up               	= stim_start_mat_temp > switch_up & stim_start_mat_temp < switch_down;
        
        if any(is_switch_up)
            whisk_stim_relay(a)     = 2;
        else
            whisk_stim_relay(a)     = 1;
        end
    end
    
    % see whether there was an LED on / offset here
    select_opto_start       = opto_firsts >= this_trial_start & opto_firsts <= this_trial_end;
    select_opto_end         = opto_lasts >= this_trial_start & opto_lasts <= this_trial_end;
    
    if sum(select_opto_start) == 1 && sum(select_opto_end) == 1
        opto_onsets(a)           = opto_firsts(select_opto_start);
        opto_offsets(a)          = opto_lasts(select_opto_end);
        opto_current_levels(a)   = max([opto_amps(select_opto_start) opto_amps(select_opto_end)]);
        opto_freq(a)             = opto_freqs(select_opto_start);
    elseif sum(select_opto_start) > 1 || sum(select_opto_end) > 1
        warning('Multiple LED stimulus values found for this trial')
        opto_onsets(a)           = min(opto_firsts(select_opto_start));
        opto_offsets(a)          = max(opto_lasts(select_opto_end));
        opto_current_levels(a)   = max(opto_powers(select_opto_start));
    elseif sum(select_opto_start) ~= sum(select_opto_end)
        warning('Mismatch in number of detected LED onsets and offsets for this trial')
    end
    
end

%%

switch expt_type % for each experiment, make sure not to split conditions by other conditions - NEEDS WORK
    case 'Velocity'
        binvec                  = [0:0.0001:2];
        [pks, locs]             = findpeaks(smooth(histc(whisk_stim_lengths,binvec),3),'MinPeakHeight',3);
        length_vals             = binvec(locs);
        whisk_stim_lengths     	= interp1(length_vals,length_vals,whisk_stim_lengths,'nearest','extrap');
    otherwise
        median_whisk_length     = nanmedian(whisk_stim_lengths);
        whisk_stim_lengths    	= repmat(median_whisk_length,size(whisk_stim_lengths));
end

whisk_stim_amplitudes       = round(whisk_stim_amplitudes / 5) * 5; % round to nearest 5%

%% opto current level fixing -- make sure current levels are grouped together into different levels sensibly

opto_current_levels(isnan(opto_current_levels)) = 0; % Set NaN values (= no LED detected) to 0

% Cluster the range of values together and set each value to the mean of its
% nearest cluster so we don't create spurious conditions
opto_current_levels     = contract_to_cluster_mean(opto_current_levels,20,'percentage');

opto_current_levels     = round(opto_current_levels); % round to nearest integer

%% Done with clean-up and event extraction; now determine the different conditions

% recover LED delays
LED_delays                  = round((opto_onsets(:) - trial_starts(:)) / LED_conditions_res,3) * LED_conditions_res;

% recover whisking delays
whisk_delays                = round((whisk_stim_onsets(:) - trial_starts(:)) / whisk_conditions_res,3) * whisk_conditions_res;

% recover LED durations
LED_ontimes                 = opto_offsets - opto_onsets;
LED_durations               = round(LED_ontimes(:) / LED_conditions_res,3) * LED_conditions_res;

% reconstruct trial matrix
trial_conditions         	= [LED_delays(:) whisk_delays(:) LED_durations(:) whisk_stim_freqs(:) round(1./whisk_stim_lengths(:)) whisk_stim_relay(:) whisk_stim_amplitudes(:) opto_current_levels(:) opto_freq(:)];

condition_headers           = {'LED start time' 'Whisk start time' 'LED duration' 'Whisk frequency' 'Whisk velocity' 'Whisk stim number' 'Whisk amplitude' 'LED Power' 'Opto_frequency'};

trial_conditions(isnan(trial_conditions)) = 999; % pass numerical flag for missing values / absent stimuli, 'unique' doesn't work well with NaNs (NaN ~= NaN)

% extract different conditions from trial matrix
[conditions, cond_inds, cond_vect]  = unique(trial_conditions,'rows');

conditions(conditions == 999)       = NaN; % replace flag with NaN again so it is clear which stimuli are absent for certain conditions

cond_nrs        = 1:size(conditions,1);

%% Get trace data (filter for spikes using 500 - 5000 Hz bandpass; can get LFP filtering e.g. with a 1-300Hz pass)

[filt_b,filt_a]       	= butter(2, spike_filt_band/(samplefreq/2));

[LFPfilt_b, LFPfilt_a] 	= butter(2, LFP_filt_band/(samplefreq/2));
LFPtraces            	= [];
LFPtimestamps           = [];

spikes(1:n_channels) 	= struct('times',[]);

%%

data_files  = [];
for a = 1:n_channels
    data_files{a}       = dir(fullfile(datafolder, filefolder, sprintf('*CH%d.continuous', get_channels(a)) ));
end

nSamples    = 1024;  % fixed to 1024 for now!

fid         = cell(n_channels, 1);

for j = 1:n_channels
    fid{j}             = fopen(fullfile(datafolder, filefolder, data_files{j}.name));
    % discard header information
    fseek(fid{j}, 1024, 0);
end

%% Common average referencing; OpenEphys reading routine partly adapted from CortexLab/Kilosort
nsamps = 0;
flag = 1;
while flag
    disp(['Processed ' num2str(round(nsamps/samplefreq)) ' seconds of ' num2str(n_channels) ' channel data...' ])
    channel_traces = zeros(nSamples * 1000, n_channels, 'int16');
    for j = 1:n_channels
        collectSamps    = zeros(nSamples * 1000, 1, 'int16');
        
        rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');
        
        nbatches        = ceil(numel(rawData)/(nSamples+6));
        for s = 1:nbatches
            
            collect_inds    = (s-1) * (nSamples + 6) + 6 + [1:nSamples];
            
            if max(collect_inds > length(rawData))
                collect_inds    = min(collect_inds:max(collect_inds));
            end
            
            rawSamps = rawData(collect_inds);
            
            
            collectSamps((s-1)*nSamples + [1:length(rawSamps)]) = rawSamps;
            
        end
        channel_traces(1:length(collectSamps),j)         = collectSamps;
    end
    
    if nbatches<1000 % we have reached the end of the file
        flag = 0;
        channel_traces = channel_traces(1:s*nSamples, :);
    end
    
    % Generate vector of timestamps based on sample number
    timestamps              = (nsamps + (1:size(channel_traces,1))) / samplefreq;
    
    % Add start time of the protocol
    timestamps              = timestamps + starttime;
    
    % Increment sample number
    nsamps                	= nsamps + size(channel_traces,1);
    
    % Subtract median of each channel from each channel
    channel_traces      	= bsxfun(@minus,channel_traces, median(channel_traces));
    
    % Convert to double
    channel_traces      	= double(channel_traces);
    
    if do_CAR
        
        CAR_channel_traces 	= bsxfun(@rdivide,channel_traces, std(channel_traces)); % Divide each trace by its standard deviation
        
        common_average      = mean(CAR_channel_traces,2);
        
        %%
        b_coeff     = NaN(1,n_channels);
        for i = 1:n_channels
            [b_coeff(i)]    = regress(channel_traces(:,i),common_average);
        end
        
        CAR_traces          = bsxfun(@times, common_average, b_coeff);
        
        channel_traces      = channel_traces - CAR_traces;
    end
    
    LFP_by_chan         = [];
    for a = 1:n_channels
        
        % Do simple homebrew spike detection by 1) filtering with bandpass, 2) smoothing and 3) thresholding
        spiketrace        	= filtfilt(filt_b,filt_a,channel_traces(:,a)); % filter data with butterworth bandpass filter to get spike traces
        
        sigma_n             = median(abs(spiketrace))/0.6745; % Robust standard deviation estimation adopted from Quian Quiroga et al. (2004)
        
        q_threshold         = (-spiketrace) > (spike_thresh * sigma_n); % determine standard deviation to determine threshold, detect threshold crossings (negative)
        spike_bool          = diff(q_threshold) == 1; % Determine instances of threshold being crossed
        these_spike_times 	= timestamps(spike_bool); % Get the timestamps of these instances
        
        spikes(a).times     = [spikes(a).times; these_spike_times(:)]; %
        
        if get_LFP
            LFPtrace        = filtfilt(LFPfilt_b,LFPfilt_a,channel_traces(:,a));
            LFPtrace        = LFPtrace(1:30:end); % resample at 1000Hz
            these_LFP_times = timestamps(1:30:end); % resample at 1000Hz
            
            LFP_by_chan    	= [LFP_by_chan LFPtrace(:)];
        end
        
    end
    LFPtraces               = [LFPtraces; LFP_by_chan];
    LFPtimestamps           = [LFPtimestamps; these_LFP_times(:)];
end

for j = 1:n_channels
    fclose(fid{j});
end


%% Sort spikes by trial and condition

cond_counters   = zeros(size(conditions,1),n_channels);
channels        = struct;   % initialise output variable

for a = 1:n_channels
    for b = 1:max(cond_nrs)
        channels(a).conditions(b).timings = conditions(b,:);
    end
end

for a = 1:n_channels
    for b = 1:length(trial_starts)
        qspiketimes                 = spikes(a).times >= trial_starts(b) & spikes(a).times < trial_ends(b);
        thesespiketimes             = spikes(a).times(qspiketimes);
        thesespiketimes             = thesespiketimes - trial_starts(b);
        
        thiscond                    = cond_vect(b);
        cond_counters(thiscond,a)   = cond_counters(thiscond,a) + 1;
        
        ephys_data.conditions(thiscond).spikes(a,cond_counters(thiscond,a),1:length(thesespiketimes(:)))  = thesespiketimes(:);
        
        % add trial time data
        if a == 1
            ephys_data.conditions(thiscond).trial_starts_timestamps(cond_counters(thiscond,a))	= trial_starts(b);
            ephys_data.conditions(thiscond).trial_ends_timestamps(cond_counters(thiscond,a))   	= trial_ends(b);
            ephys_data.conditions(thiscond).trial_starts(cond_counters(thiscond,a))             = trial_starts(b) - starttime;
            ephys_data.conditions(thiscond).trial_ends(cond_counters(thiscond,a))               = trial_ends(b) - starttime;
        end
        
        % LFP processing. Is it better to simply do this by condition, and
        % then have a matrix of data x episodes x data length?
        if get_LFP
            qLFPtimestamps              = LFPtimestamps >= trial_starts(b) & LFPtimestamps < trial_ends(b);
            ephys_data.conditions(thiscond).LFP_trace(a,cond_counters(thiscond,a),1:length(LFPtraces(qLFPtimestamps,a))) = LFPtraces(qLFPtimestamps,a);
            ephys_data.LFP_timestamps   = 1:sum(qLFPtimestamps) / 30000;
        end
        if a == 1 && save_sync_chans
            q_sync_timestamps           = sync_chan_timestamps >= trial_starts(b) & sync_chan_timestamps < trial_ends(b);
            for i = 1:size(sync_chan_traces,2)
                ephys_data.conditions(thiscond).sync_chan_traces(i,cond_counters(thiscond,a),1:sum(q_sync_timestamps)) = sync_chan_traces(q_sync_timestamps,i);
            end
        end
    end
end

for i = 1:length(ephys_data.conditions)
    ephys_data.conditions(i).values             = conditions(i,:);
    ephys_data.conditions(i).n_trials           = size(ephys_data.conditions(i).spikes,2);
    ephys_data.conditions(i).whisk_onset        = conditions(i,2);
    ephys_data.conditions(i).whisk_stimulator 	= conditions(i,6);
    ephys_data.conditions(i).whisk_amplitude    = conditions(i,7);
    ephys_data.conditions(i).whisk_velocity     = conditions(i,5);
    ephys_data.conditions(i).whisk_frequency    = conditions(i,4);
    ephys_data.conditions(i).LED_onset          = conditions(i,1);
    ephys_data.conditions(i).LED_power          = conditions(i,8);
    ephys_data.conditions(i).LED_duration       = conditions(i,3);
    ephys_data.conditions(i).opto_freq          = conditions(i,9);
    ephys_data.conditions(i).spikes(ephys_data.conditions(i).spikes == 0) = NaN; % replace empty values for spike matrix (default to 0) with NaN
end

ephys_data.condition_values     = conditions;
ephys_data.condition_names      = condition_headers;
ephys_data.trial_length         = trial_length;
ephys_data.trial_interval       = total_length;
ephys_data.block_length         = total_length * length(ephys_data.conditions);
ephys_data.protocol_duration    = total_length * length(trial_starts);
ephys_data.parameters           = parameters;
ephys_data.data_folder          = filefolder;
ephys_data.channelmap           = get_channels;
ephys_data.start_time           = starttime;
ephys_data.end_time             = endtime;
ephys_data.rec_length           = endtime - starttime;


