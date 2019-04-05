function [sync_data] = get_stim_sync_data(datafolder,events_chans,trials_from_whisk, whisk_buffer)
% function [sync_data] = get_stim_sync_data(DATAFOLDER,EVENTS_CHANS,TRIALS_FROM_WHISK, WHISK_BUFFER)
% 
% Get the PulsePal stimulus synchronisation data for an openephys experiment; 
% trial onset and offset times, whisk stimulus on and offset, opto stimulus 
% on and offset, amplitudes, frequencies and durations, etc.
% 
% SYNC_DATA:
% A struct with fairly self-explanatory field names;
% SYNC_DATA.conditions additionally provides the trial sync data organised 
% by the set of stimulus conditions.
% 
% 
% DATAFOLDER: Full path to a folder containing the raw openephys data for 
% this experiment (a set of '.continuous' files).
% 
% 
% EVENTS_CHANS:
% a 4-element vector of input channel numbers in the following order:
% [TRIAL_CHAN STIM_ON_CHAN OPTO_CHAN STIM_NR_CHAN]
% TRIAL_CHAN: the ADC channel that keeps track of trials (high = trial);
% STIM_ON_CHAN: the ADC channel that keeps track of when the whisker stimulator is on (high = on);
% OPTO_CHAN: the ADC channel that keeps track of whether the opto stimulation (LED or LASER) is on (high = on);
% STIM_NR_CHAN: the ADC channel that keeps track of which stimulator is being used; 1 or 2.
% 
% 
% TRIALS_FROM_WHISK: When trial sync data is incomplete or unsatisfactory,
% determine 'trials' based on a period of time around the whisker stimulus,
% the period of time is set by WHISK_BUFFER.
% 
% 
% WHISK_BUFFER: If TRIALS_FROM_WHISK == true, then the script will generate
% trials that start at whisk_onset - WHISK_BUFFER and last until whisk_onset
% + WHISK_BUFFER. Whisk_buffer is also used to determine whether consecutive
% whisks should be considered part of the same burst within a trial, or whether 
% they should constitute a new trial;
% if inter-whisk_interval > whisk_buffer --> new trial.
% 

if nargin < 3
    trials_from_whisk   = false;
end
if nargin < 4
    whisk_buffer        = 2; % default is to collect 2 secs from the onset of whisk stimulus
end

% Unpack 
trial_input_nr          = events_chans(1);         % Which input channel has the trial TTL
stim_input_nr           = events_chans(2);         % Which input channel has the stim / whisk TTL
opto_input_nr       	= events_chans(3);           % Which input channel has the LED TTL
switch_input_nr         = events_chans(4);  	% Which input channel switches between stimulators?

% Constants?
opto_conditions_res   	= 1;    % resolution in ms for automatically extracting conditions from LED delays
whisk_conditions_res    = 10;  % resolution in ms for automatically extracting conditions from whisker delays

%% Start collecting events data

if switch_input_nr ~= 0
    trial_threshold     = 0.25; % For trial, signal goes up to 2.5V or less
    stim_threshold      = 2.7; % For whisking (always during trial), signal goes up to 2.75 - 5V
    opto_threshold    	= 0.1; % normal TTL logic - 0 to 5V
    switch_threshold    = 2.5; % normal TTL logic - 0 to 5V
else
    trial_threshold     = 0.25; % normal TTL logic - 0 to 5V
    stim_threshold      = 2.5; % normal TTL logic - 0 to 5V
    opto_threshold    	= 0.1; % LED is no longer TTL, voltage varies with power (with 5V representing max); 0.05V thresh will detect events above ~1% max power
    switch_threshold    = 2.5; % normal TTL logic - 0 to 5V
end

adc_channel_nrs        	= [trial_input_nr stim_input_nr opto_input_nr switch_input_nr];
adc_channel_thresholds 	= [trial_threshold stim_threshold opto_threshold switch_threshold];

%% Some code to find file prefix for data and sync files

data_files 	= dir(fullfile(datafolder,'*.continuous'));
data_files	= {data_files.name}';

is_adc_file = ~cellfun(@isempty,regexp(data_files,'ADC')); % Remove ADC files

chan_files  = data_files(is_adc_file); %
chan_inds   = cell2mat(regexp(data_files,'CH')); % find index of occurrence of 'CH'

chan_ind    = chan_inds(1) - 1; % We only need to look at the prefix for one file
chan_file   = chan_files{1}; %

data_prefix = chan_file(1:chan_ind);

%% 
for a = 1:length(adc_channel_nrs) % loop through the analog input channels
    
    if a == 1 || adc_channel_nrs(a) ~= adc_channel_nrs(a-1) % Don't reload data if trace is already loaded
        disp(['Loading ADC input channel ' num2str(adc_channel_nrs(a))])
        disp(['File ' datafolder filesep data_prefix 'ADC' num2str(adc_channel_nrs(a)) '.continuous'])
        [thisTTL timestamps info] = load_open_ephys_data([datafolder filesep data_prefix 'ADC' num2str(adc_channel_nrs(a)) '.continuous']);
        
        % thisTTL         = thisTTL - min(thisTTL(:));
        
        starttime       = min(timestamps);              % find start time
        endtime         = max(timestamps);              % find end time
        timestamps      = (1:length(thisTTL)) / 30000;  % manually create new timestamps at 30kHz, openephys sometimes suffers from timestamp wobble even though data acquisition is spot on
        timestamps      = timestamps + starttime;       % add start time to the newly created set of timestamps
    end
    
    thisTTL_bool   	= thisTTL > adc_channel_thresholds(a); % find where the TTL signal is 'high'
    
    start_inds      = find(diff(thisTTL_bool) > 0.5);   % find instances where the TTL goes from low to high
    end_inds        = find(diff(thisTTL_bool) < -0.5);  % find instances where the TTL goes from high to low
    
    start_times 	= timestamps(start_inds);   % find the timestamps of start events
    end_times    	= timestamps(end_inds);     % find the timestamps of end events
    
    if ~isempty(start_times) % Some channels may not have events (e.g. stim switch channel if only 1 stimulator used)
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
                stim_amps(i)   	= ((median(stim_segment) - 2.5) / 2.5) * 100; % Stimulus amplitude in % of max
            end
        case 3
            opto_starts      = start_times(:);
            opto_ends        = end_times(:);
            
            opto_powers      = NaN(size(start_inds));
            for i = 1:length(start_inds)
                stim_segment    = thisTTL(start_inds(i):end_inds(i));
                opto_powers(i) 	= median(stim_segment) / 5 * 100; % Stimulus amplitude in % of max
            end
        case 4
            switch_up       = start_times(:);
            switch_down     = end_times(:);
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

total_length 	= round(median(diff(trial_starts)),1);

%% Work out the velocity (length) of a whisk, and the frequency:

% Find which whisking onsets are the first of a trial, and which onsets
% are the last of a trial
allwhisk_firstvect          = stim_starts(find(diff(stim_starts) > (trial_length/2))+1);
allwhisk_lastvect           = stim_starts(find(diff(stim_starts) > (trial_length/2)));

first_stim_amps             = stim_amps(find(diff(stim_starts) > (trial_length/2))+1);
first_stim_amps             = [stim_amps(1); first_stim_amps(:)];

whisk_starts            	= [stim_starts(1); allwhisk_firstvect(:)];
whisk_lasts              	= [allwhisk_lastvect(:); stim_starts(end)];

whisk_ends               	= stim_ends(find(diff(stim_starts) > (trial_length/2))+1);
whisk_ends                	= [stim_ends(1); whisk_ends(:)];

whisk_lengths              	= whisk_ends - whisk_starts;
whisk_freqs              	= NaN(size(whisk_starts));

for a = 1:length(whisk_starts)
    this_whisk_start    = whisk_starts(a);
    this_whisk_end      = whisk_lasts(a);
    q_whisks            = stim_starts > this_whisk_start & stim_starts < this_whisk_end;
    
    this_whisk_freq   	= mean(round(1./diff(stim_starts(q_whisks))));
    if isempty(this_whisk_freq)
        this_whisk_freq = 99;
    elseif isnan(this_whisk_freq)
        this_whisk_freq = 99;
    end
    whisk_freqs(a)      = this_whisk_freq;
end

stim_starts             = whisk_starts;
stim_amps               = first_stim_amps;

%% Match events to trials
ntrials                 = length(trial_starts);

whisk_stim_onsets       = NaN(size(trial_starts));
whisk_stim_lengths      = NaN(size(trial_starts));
whisk_stim_freqs        = NaN(size(trial_starts));
whisk_stim_relay        = NaN(size(trial_starts));
whisk_stim_amplitudes   = NaN(size(trial_starts));

opto_onsets              = NaN(size(trial_starts));
opto_offsets             = NaN(size(trial_starts));
opto_current_levels      = NaN(size(trial_starts));

for a = 1:ntrials
    this_trial_start    = trial_starts(a);
    this_trial_end      = trial_ends(a);
    
    % see whether there was a whisker stimulus
    select_whisk_start 	= whisk_starts >= this_trial_start & whisk_starts <= this_trial_end;
    
    if sum(select_whisk_start) == 1
        whisk_stim_onsets(a)        = stim_starts(select_whisk_start);
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
    elseif sum(select_whisk_start) > 1
        error('Multiple whisker stimulus values found for this trial')
    end
    
    % see whether there was an LED on / offset here
    select_opto_start   = opto_starts >= this_trial_start & opto_starts <= this_trial_end;
    select_opto_end 	= opto_ends >= this_trial_start & opto_ends <= this_trial_end;
    
    if sum(select_opto_start) == 1 && sum(select_opto_end) == 1
        opto_onsets(a)           = opto_starts(select_opto_start);
        opto_offsets(a)          = opto_ends(select_opto_end);
        opto_current_levels(a)   = opto_powers(select_opto_start);
    elseif sum(select_opto_start) > 1 || sum(select_opto_end) > 1
        error('Multiple opto stimulus values found for this trial')
    elseif sum(select_opto_start) ~= sum(select_opto_end)
        error('Mismatch in number of detected opto onsets and offsets for this trial')
    end
    
end

%% Whisker stim length

% Find whisker stim lengths; make histogram of all stim length values, find
% the peaks in the histogram, and then get rid of jitter in timing data by 
% rounding everything to those peak values
binvec                      = [0:0.0001:2];
[pks, locs]                 = findpeaks(smooth(histc(whisk_stim_lengths,binvec),3),'MinPeakHeight',3);
length_vals                 = binvec(locs);
if numel(length_vals) == 1
    % a single whisker stimulation length; just grab median to remove jitter
    median_whisk_length     = nanmedian(whisk_lengths);
    whisk_stim_lengths   	= repmat(median_whisk_length,size(whisk_stim_lengths));
else
    % multiple stimulus velocities for this experiment; set each value to 
    % its nearest peak / 'local median' value to get rid of jitter
    whisk_stim_lengths          = interp1(length_vals,length_vals,whisk_stim_lengths,'nearest','extrap');
end

whisk_stim_amplitudes       = round(whisk_stim_amplitudes);     % round to nearest 1% to remove jitter
opto_current_levels         = round(opto_current_levels);       % round to nearest 1% to remove jitter

%% Done with clean-up and event extraction; now determine the different conditions

% recover LED delays
opto_delays                 = round((opto_onsets(:) - trial_starts(:)) / opto_conditions_res,3) * opto_conditions_res;

% recover whisking delays
whisk_delays                = round((whisk_stim_onsets(:) - trial_starts(:)) / whisk_conditions_res,3) * whisk_conditions_res;

% recover LED durations
opto_durations           	= opto_offsets - opto_onsets;
opto_durations              = round(opto_durations(:) / 5,3) * 5;

% reconstruct trial matrix
trial_conditions         	= [whisk_delays(:) whisk_stim_relay(:)  whisk_stim_amplitudes(:)  whisk_stim_freqs(:) round(1./whisk_stim_lengths(:)) opto_delays(:) opto_durations(:) opto_current_levels(:)];

trial_conditions(isnan(trial_conditions))   = 999; % pass numerical flag for missing values / absent stimuli, 'unique' doesn't work well with NaNs (NaN ~= NaN)

% extract different conditions from trial matrix
[conditions, cond_inds, cond_vect]  = unique(trial_conditions,'rows');

conditions(conditions == 999)               = NaN; % replace flag with NaN again so it is clear which stimuli are absent for certain conditions
trial_conditions(trial_conditions == 999)   = NaN;

%% We've got all the info; start constructing the output variable:

% Things that vary by condition
for a = 1:length(conditions)
    q_cond_trials                           = cond_vect == a; 
    sync_data.conditions(a).trial_starts  	= trial_starts(q_cond_trials);
    sync_data.conditions(a).trial_ends    	= trial_ends(q_cond_trials);
    
    sync_data.conditions(a).whisk_starts 	= whisk_stim_onsets(q_cond_trials);
    sync_data.conditions(a).opto_onsets  	= opto_onsets(q_cond_trials);
    sync_data.conditions(a).opto_offsets  	= opto_offsets(q_cond_trials);
    
    sync_data.conditions(a).whisk_start   	= conditions(a,1);
    sync_data.conditions(a).whisk_stim_nr 	= conditions(a,2);
    sync_data.conditions(a).whisk_amp     	= conditions(a,3);
    sync_data.conditions(a).whisk_freq    	= conditions(a,4);
    sync_data.conditions(a).whisk_velocity	= conditions(a,5);
    
    sync_data.conditions(a).opto_onset    	= conditions(a,6);
    sync_data.conditions(a).opto_duration 	= conditions(a,7);
    sync_data.conditions(a).opto_power    	= conditions(a,8);
    
    sync_data.conditions(a).n_trials       	= sum(q_cond_trials);
end

% Things that are true across the entire experiment:
sync_data.data_folder     	= datafolder; % full file path to original raw data location
sync_data.rec_start_time  	= starttime;
sync_data.rec_end_time    	= endtime;
sync_data.rec_duration    	= endtime - starttime;
sync_data.trial_length   	= trial_length;
sync_data.trial_interval  	= total_length;
sync_data.block_length   	= total_length * length(conditions);
sync_data.protocol_duration = total_length * length(trial_starts);
sync_data.trial_starts      = trial_starts(:);
sync_data.trial_ends        = trial_ends(:);
sync_data.whisk_starts      = whisk_delays(:);
sync_data.whisk_stim_nr 	= whisk_stim_relay(:);
sync_data.whisk_amps        = whisk_stim_amplitudes(:);
sync_data.whisk_burst_freq  = whisk_freqs(:);
sync_data.whisk_veloc       = round(1./whisk_stim_lengths(:));
sync_data.opto_delays       = opto_delays(:);
sync_data.opto_durations    = opto_durations(:);
sync_data.opto_power        = opto_current_levels(:);

