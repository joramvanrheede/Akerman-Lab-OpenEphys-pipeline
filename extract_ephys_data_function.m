function [channels] = extract_ephys_data_function(datafolder,datafilenr,qspike_detection,parameters)

data_prefix             = parameters.data_prefix;           % = '100';
q_digital_events        = parameters.digital_events;        % 1/0% events determined using openephys event detection
q_channelmap            = parameters.channelmap;            % 1/0 have channels been mapped in openephys (1) or do we need to remap ourselves (0)?
spike_smoothwin         = parameters.spike_smoothwin;       % smoothing window for spike detection
spike_thresh            = parameters.spike_thresh;          % spike threshold in SDs
LED_conditions_res      = parameters.LED_conditions_res;    % resolution in ms for automatically extracting conditions from LED delays
whisk_conditions_res    = parameters.whisk_conditions_res;  % resolution in ms for automatically extracting conditions from whisker delays
spontwin                = parameters.spontwin;              % window for spont spikes e.g. [0 200] ms

trial_input_nr          = parameters.trial_channel;         % Which input channel has the trial TTL
stim_input_nr           = parameters.whisk_channel;         % Which input channel has the stim / whisk TTL
LED_input_nr            = parameters.LED_channel;           % Which input channel has the LED TTL
switch_input_nr         = parameters.stim_switch_channel;  	% Which input channel switches between stimulators?
q_override              = parameters.override_conds;        % Override TTL-based condition structure?
n_conds                 = parameters.n_conds;               % Number of conditions if using override

qmorse                  = parameters.morse;                 % use trial start/end morse code?

expt_type               = parameters.experiment_type;       % what type of experiment is this?

%% relatively fixed variables

% IMPORTANT - this is the channel order for a 16-channel NeuroNexus probe
% sent through the A16 --> Omnetics32 adapter and acquired with the openephys
% headstage: [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]
get_channels            = [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]; % IMPORTANT [20 4 29 13 18 2 30 14 17 1 32 16 31 15 19 3]

morse_start             = [3 1 3 1 3]; % morse code for start of trial; short = 1, long = 3;
morse_stop              = [1 1 1 3 1 3]; % morse code for end of trial; short = 1, long = 3;

%% some file I/O pre-work

n_channels              = length(get_channels);

data_contents           = dir(datafolder); % find what is in the data folder
filefolders             = data_contents([data_contents.isdir]); % get only directories
filefolders             = {filefolders.name}'; % get the names of the folders
filefolders             = filefolders(3:end); % remove '.' and '..' folders

pattern                 = '._\d+$'; % regular expression to search for file numbers below

[startinds, endinds]    = regexp(filefolders, pattern, 'start','end'); % find number pattern in the remaining filefolders
filenumbers             = {};

% generate a vector of file numbers from the files in the data folder
for i = 1:length(startinds)
    filenumbers{i} = filefolders{i}(startinds{i}+2:endinds{i});
end

fileind                 = find(strcmp(filenumbers,num2str(datafilenr))); % find index of target data folder
filefolder              = filefolders{fileind}; % this is the folder we're after

%% Start collecting events data

channels            = struct;   % initialise output variable
cond_counters       = [];       % initialise counter per condition per channel;

if q_digital_events % do we use digital events recorded by openephys? if so use 'all_channels.events' file
    [events timestamps info] = load_open_ephys_data([datafolder filesep filefolder filesep 'all_channels.events']);
    
    % get events for trial, stim, LED. input nr - 1 because of sensible
    % vs. programmer counting
    trial_starts    = timestamps(events == trial_input_nr - 1 & info.eventId == 1);
    trial_ends      = timestamps(events == trial_input_nr - 1 & info.eventId == 0);
    stim_starts  	= timestamps(events == stim_input_nr - 1 & info.eventId == 1);
    stim_ends       = timestamps(events == stim_input_nr - 1 & info.eventId == 0);
    LED_starts      = timestamps(events == LED_input_nr - 1 & info.eventId == 1);
    LED_ends        = timestamps(events == LED_input_nr - 1 & info.eventId == 0);
    switch_up       = timestamps(events == switch_input_nr - 1 & info.eventId == 1);
    switch_down     = timestamps(events == switch_input_nr - 1 & info.eventId == 0);
    
    % some cleaning up 
    trial_starts(trial_starts > trial_ends(end))   = [];
    trial_ends(trial_ends < trial_starts(1))       = [];
    
else % events need to be extracted manually from the analog input signal. 
    
    if switch_input_nr ~= 0
        trial_threshold     = 5/3; % For trial, signal goes up to 2.5V
        stim_threshold      = (5/3) * 2; % For whisking (always during trial), signal goes up to 5V
        LED_threshold       = 2.5; % normal TTL logic - 0 to 5V
        switch_threshold    = 2.5; % normal TTL logic - 0 to 5V
    else
        trial_threshold     = 2.5; % normal TTL logic - 0 to 5V
        stim_threshold      = 2.5; % normal TTL logic - 0 to 5V
        LED_threshold       = 2.5; % normal TTL logic - 0 to 5V
        switch_threshold    = 2.5; % normal TTL logic - 0 to 5V
    end
    
    adc_channel_nrs        	= [trial_input_nr stim_input_nr LED_input_nr switch_input_nr];
    adc_channel_thresholds 	= [trial_threshold stim_threshold LED_threshold switch_threshold];
    
    for a = 1:4 % loop through the analog input channels
        
        if a == 1 || adc_channel_nrs(a) ~= adc_channel_nrs(a-1) % Don't reload data if trace is already loaded
            disp(['Loading ADC input channel ' num2str(adc_channel_nrs(a))])
            [thisTTL timestamps info] = load_open_ephys_data([datafolder filesep filefolder filesep data_prefix '_ADC' num2str(adc_channel_nrs(a)) '.continuous']);
            
            starttime       = min(timestamps); % find start time
            timestamps      = (1:length(thisTTL)) / 30000; % manually create new timestamps at 30kHz, openephys sometimes suffers from timestamp wobble even though data acquisition is spot on
            timestamps      = timestamps + starttime; % add start time to the newly created set of timestamps
        end
        
        thisTTL_bool   	= thisTTL > adc_channel_thresholds(a); % find where the TTL signal is 'high'
        
        start_inds      = find(diff(thisTTL_bool) > 0.5); % find instances where the TTL goes from low to high
        end_inds        = find(diff(thisTTL_bool) < -0.5); % find instances where the TTL goes from high to low
        
        start_times 	= timestamps(start_inds); % find the timestamps of start events
        end_times    	= timestamps(end_inds); % find the timestamps of end events
        
        end_times(end_times < start_times(1))       = []; % discard potential initial end without start
        start_times(start_times > end_times(end))   = []; % discard potential final start without end
        
        switch a % this determines what the start and end timestamps should be assigned to: trial/trial, LED/opto stim or stim/whisk stim.
            case 1
                trial_starts    = start_times(:);
                trial_ends      = end_times(:);
            case 2
                stim_starts 	= start_times(:);
                stim_ends     	= end_times(:);
            case 3
                LED_starts      = start_times(:);
                LED_ends        = end_times(:);
            case 4
                switch_up       = start_times(:);
                switch_down     = end_times(:);
        end
    end
end

% if we're using a morse signal to detect start and end of the protocol,
% detect morse code for start and end of trials
if qmorse
    trial_lengths 	= trial_ends - trial_starts;
    morse_lengths   = round(trial_lengths * 10);
    
    start_ind       = strfind(morse_lengths(:)',morse_start) + length(morse_start);
    end_ind         = strfind(morse_lengths(:)',morse_stop) - 1;
    
    trial_starts    = trial_starts(start_ind:end_ind);
    trial_ends      = trial_ends(start_ind:end_ind);
end

%% A lot of cleanup and repair from here

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

% Correct for odd TTL pulses
trial_start     = min(trial_starts(:));
trial_end       = max(trial_ends(:));

% ensure we only look for events happening during trials;
% make sure all starts have matching ends and  vice versa
LED_starts(LED_starts < trial_start)     	= [];
LED_ends(LED_ends < trial_start)           	= [];

stim_starts(stim_starts < trial_start)    	= [];
stim_ends(stim_ends < trial_start)        	= [];

LED_starts(LED_starts > trial_end)        	= [];
LED_ends(LED_ends > trial_end)           	= [];

stim_starts(stim_starts > trial_end)       	= [];
stim_ends(stim_ends > trial_end)          	= [];

stim_starts(stim_starts > max(stim_ends))	= [];
stim_ends(stim_ends < min(stim_starts))   	= [];


% more cleaning up; assuming similar trial durations, see if there are any
% exceptionally long or short trials at the start or end (can happen
% because of floating voltages on the arduino pins)
qtrial2         = round(diff(trial_starts),1) == total_length; % compare 
anomalind       = find(~qtrial2); % indices of anomalous lengths

anomalearly    	= anomalind(anomalind < (.5 * length(trial_starts))); % find which anomalous indices are happenig early on in the trial sequence
anomalate      	= anomalind(anomalind > (.5 * length(trial_starts))); % find which are happening late

anomalate       = anomalate + 1;

anomalous       = [anomalearly(:); anomalate(:)];

trial_starts(anomalous)     = []; % remove anomalous trial starts
trial_ends(anomalous)       = []; % remove anomalous trial ends

% 
allwhisks                   = stim_starts;
allwhisk_ends               = stim_ends;

%% work out the velocity (length) of a whisk, and the frequency:

% Find which whisking onsets are the first of a trial, and which onsets
% are the last of a trial
allwhisk_firstvect          = allwhisks(find(diff(allwhisks) > (trial_length/2))+1);
allwhisk_lastvect           = allwhisks(find(diff(allwhisks) > (trial_length/2)));

whisk_starts            	= [allwhisks(1); allwhisk_firstvect(:)];
whisk_lasts              	= [allwhisk_lastvect(:); allwhisks(end)];

whisk_ends               	= allwhisk_ends(find(diff(allwhisks) > (trial_length/2))+1);
whisk_ends                	= [allwhisk_ends(1); whisk_ends(:)];

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

allwhisks                                       = whisk_starts;

allwhisks(allwhisks < min(trial_starts))      = [];
allwhisks(allwhisks > max(trial_ends))        = [];
whisk_freqs(allwhisks < min(trial_starts))    = [];
whisk_freqs(allwhisks > max(trial_ends))      = [];
whisk_lengths(allwhisks < min(trial_starts))  = [];
whisk_lengths(allwhisks > max(trial_ends))    = [];

LED_starts(LED_starts < min(trial_starts))    = [];
LED_starts(LED_starts > max(trial_ends))      = [];
LED_ends(LED_ends < min(trial_starts))        = [];
LED_ends(LED_ends > max(trial_ends))          = [];

%% fill in any missing trial start and end events

trial_interval      = median(diff(trial_starts));
qtrialgap           = [0; diff(trial_starts) > trial_interval + 1];

while any(qtrialgap)
    trial_ind       = find(qtrialgap,1);
    trial_starts  = [trial_starts(1:trial_ind-1); trial_starts(trial_ind-1)+trial_interval;  trial_starts(trial_ind:end)];
    qtrialgap       = [0; diff(trial_starts) > trial_interval + 1];
end

qtrialgap = [0; diff(trial_ends) > trial_interval + 1];

while any(qtrialgap)
    trial_ind       = find(qtrialgap,1);
    trial_ends    = [trial_ends(1:trial_ind-1); trial_ends(trial_ind-1)+trial_interval;  trial_ends(trial_ind:end)];
    qtrialgap       = [0; diff(trial_ends) > trial_interval + 1];
end

%% Check for trials without corresponding events and remove

% find the nr of trials
ntrials         = length(trial_starts);

% Check for trials without corresponding events
for j = 1:length(trial_starts)
    this_start  = trial_starts(j);
    this_end    = trial_ends(j);
    
    qLEDstart  	= any(LED_starts > this_start & LED_starts < this_end);
    qwhisk      = any(allwhisks > this_start & allwhisks < this_end);
    
    if ~qLEDstart & ~qwhisk
        trial_starts(j)   = NaN;
        trial_ends(j)     = NaN;
        warning(['No LED or whisker stim found on trial ' num2str(j)])
    end
    if ~qLEDstart & qwhisk
        allwhisks(j)        = NaN;
        whisk_freqs(j)      = NaN;
        whisk_lengths(j)    = NaN;
        trial_starts(j)   = NaN;
        trial_ends(j)     = NaN;
        warning(['No LED found on trial ' num2str(j)])
    end
    if qLEDstart & ~qwhisk
        trial_starts(j)   = NaN;
        trial_ends(j)     = NaN;
        LED_starts(j)       = NaN;
        LED_ends(j)         = NaN;
        warning(['No whisker stim found on trial ' num2str(j)])
    end
end

allwhisks(isnan(allwhisks))             = [];
whisk_freqs(isnan(whisk_freqs))         = [];
whisk_lengths(isnan(whisk_lengths))     = [];
trial_starts(isnan(trial_starts))   = [];
trial_ends(isnan(trial_ends))       = [];
LED_starts(isnan(LED_starts))           = [];
LED_ends(isnan(LED_ends))               = [];

switch expt_type % for each experiment, make sure not to split conditions by other conditions - NEEDS WORK
    case 'Timing'
        median_whisk_length     = nanmedian(whisk_lengths);
        whisk_lengths          	= repmat(median_whisk_length,size(whisk_lengths));
    case 'Velocity'
    case 'Drive'
        median_whisk_length     = nanmedian(whisk_lengths);
        whisk_lengths         	= repmat(median_whisk_length,size(whisk_lengths));
    case 'Frequency'
        median_whisk_length     = nanmedian(whisk_lengths);
        whisk_lengths         	= repmat(median_whisk_length,size(whisk_lengths));
end


%% Done with clean-up and event extraction; now determine the different conditions

% Which stimulator was used?

stim_nr                     = ones(size(trial_starts)); % initialise stimulator nr to '1' for all trials
[vals,up_stim_trials,ib]    = intersect(trial_starts, switch_up); % find instances of the stim switch going up at the same time as trial starts - those are the trials where stimulator 2 is being used
stim_nr(up_stim_trials)     = 2; % set those instances to stimulator 2

% recover LED delays
LED_delays      = round((LED_starts(:) - trial_starts(:)) / LED_conditions_res,3) * LED_conditions_res;

% recover whisking delays
whisk_delays    = round((allwhisks(:) - trial_starts(:)) / whisk_conditions_res,3) * whisk_conditions_res;

% recover LED durations
LED_ontimes     = LED_ends - LED_starts;
LED_durations   = round(LED_ontimes(:) / 10,3) * 10;

% reconstruct trial matrix
trial_timings   = [LED_delays(:) whisk_delays(:) LED_durations(:) whisk_freqs(:) round(1./whisk_lengths(:)) stim_nr(:)];

% extract different conditions from trial matrix
[conditions, cond_inds, cond_vect]   = unique(trial_timings,'rows');

if q_override
    conditions  = zeros(n_conds,5);
    cond_vect   = repmat([1:n_conds]',length(cond_vect)/n_conds,1);
end

cond_nrs        = 1:size(conditions,1);

if isempty(cond_counters)
    cond_counters   = zeros(size(conditions,1),n_channels);
end

for a = 1:n_channels
    for b = 1:max(cond_nrs)
        channels(a).conditions(b).timings = conditions(b,:);
    end
end

%% Get trace data (filter for spikes using 300 - 6000 Hz bandpass; can get LFP filtering e.g. with a 1-300Hz pass)
% initialise butterworth bandpass filter
[filt_b,filt_a]           = butter(2, [300 6000]/(30000/2));

if qspike_detection % 'manual' spike detection in matlab

    for a = 1:n_channels
        disp(['Loading channel ' num2str(a)]);
        [thistrace timestamps info] = load_open_ephys_data([datafolder filesep filefolder filesep data_prefix '_CH' num2str(get_channels(a)) '.continuous']);
        
        starttime          	= min(timestamps);
        
        timestamps          = ((1:length(thistrace)) / 30000); % generate timestamps
        
        % Do simple homebrew spike detection by 1) filtering with bandpass, 2) smoothing and 3) thresholding
        thistrace           = filter(filt_b,filt_a,thistrace); % filter data with butterworth bandpass filter to get spike traces
        thistrace           = smooth(thistrace, spike_smoothwin); % smooth data with
        
        q_threshold         = thistrace > (spike_thresh * std(thistrace)); % determine standard deviation to determine threshold
        spike_bool          = diff(q_threshold) == 1; % Determine instances of threshold being crossed 
        these_spike_times 	= timestamps(spike_bool); % Get the timestamps of these instances
        
        spikes(a).times     = these_spike_times + starttime; % 
        spikes(a).waveforms = [];
        spikes(a).info      = [];
    end
else  % use spikes detected from openephys .spikes file instead
    spikes  = []; % initialise spikes. Will become a 1xmaxchannelnr struct
    for a = 1:n_channels
        if q_channelmap
            [spike_waves these_spike_times info] = load_open_ephys_data([datafolder filesep filefolder filesep 'SE' num2str(a-1) '.spikes']);
        else
            [spike_waves these_spike_times info] = load_open_ephys_data([datafolder filesep filefolder filesep 'SE' num2str(get_channels(a)-1) '.spikes']);
        end
        
        spikes(a).times     = these_spike_times;
        spikes(a).waveforms = spike_waves;
        spikes(a).info      = info;
    end
end


%% Sort spikes by trial and condition
for a = 1:n_channels
    
    chan_spike_times    = [];
    for b = 1:length(trial_starts)
        qspiketimes                 = (spikes(a).times >= trial_starts(b)) & (spikes(a).times < trial_ends(b));
        thesespiketimes             = spikes(a).times(qspiketimes);
        thesespiketimes             = thesespiketimes - trial_starts(b);
        
        thiscond                    = cond_vect(b);
        cond_counters(thiscond,a)   = cond_counters(thiscond,a) + 1;
        
        channels(a).conditions(thiscond).episodes(cond_counters(thiscond,a)).spikes  = thesespiketimes(:);
        chan_spike_times            = [chan_spike_times; thesespiketimes(:)];
    end
    
    %% spont spike rate
    spontspikes     = chan_spike_times(chan_spike_times >= (spontwin(1) / 1000) & chan_spike_times < (spontwin(2) / 1000));
    spontwinsize    = (spontwin(2) - spontwin(1)) / 1000;           % size of window in seconds
    channels(a).spontspikerate   = length(spontspikes) / spontwinsize / ntrials; % spontaneous spike rate by channel
end
