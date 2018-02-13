function [channels] = extract_ephys_data_function(datafolder,datafilenr,qspike_detection,parameters)

data_prefix             = parameters.data_prefix;           % = '100';
q_digital_events        = parameters.digital_events;        % 1/0% events determined using openephys event detection
q_channelmap            = parameters.channelmap;            % 1/0 have channels been mapped in openephys (1) or do we need to remap ourselves (0)?
spike_smoothwin         = parameters.spike_smoothwin;       % smoothing window for spike detection
spike_thresh            = parameters.spike_thresh;          % spike threshold in SDs
LED_conditions_res      = parameters.LED_conditions_res;    % resolution in ms for automatically extracting conditions from LED delays
whisk_conditions_res    = parameters.whisk_conditions_res;  % resolution in ms for automatically extracting conditions from whisker delays
spontwin                = parameters.spontwin;              % window for spont spikes e.g. [0 200] ms

episode_input_nr        = parameters.trial_channel;         % Which input channel has the episode TTL
piezo_input_nr          = parameters.whisk_channel;         % Which input channel has the piezo / whisk TTL
LED_input_nr            = parameters.LED_channel;           % Which input channel has the LED TTL
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
    
    % get events for episode, piezo, LED. input nr - 1 because of sensible
    % vs. programmer counting
    episode_starts  = timestamps(events == episode_input_nr - 1 & info.eventId == 1);
    episode_ends    = timestamps(events == episode_input_nr - 1 & info.eventId == 0);
    piezo_starts  	= timestamps(events == piezo_input_nr - 1 & info.eventId == 1);
    piezo_ends      = timestamps(events == piezo_input_nr - 1 & info.eventId == 0);
    LED_starts      = timestamps(events == LED_input_nr - 1 & info.eventId == 1);
    LED_ends        = timestamps(events == LED_input_nr - 1 & info.eventId == 0);
    
    % some cleaning up 
    episode_starts(episode_starts > episode_ends(end))   = [];
    episode_ends(episode_ends < episode_starts(1))       = [];
    
else % events need to be extracted manually from the analog input signal. 
    for a = 1:3 % loop through the analog input channels
        disp(['Loading ADC input channel ' num2str(a)])
        [thisTTL timestamps info] = load_open_ephys_data([datafolder filesep filefolder filesep data_prefix '_ADC' num2str(a) '.continuous']);
        
        starttime       = min(timestamps); % find start time
        timestamps      = (1:length(thisTTL)) / 30000; % manually create new timestamps at 30kHz, openephys sometimes suffers from timestamp wobble even though data acquisition is spot on
        timestamps      = timestamps + starttime; % add start time to the newly created set of timestamps
        
        thisTTL         = thisTTL > 2.5; % find where the TTL signal is 'high'
        
        start_inds      = find(diff(thisTTL) > 0.5); % find instances where the TTL goes from low to high
        end_inds        = find(diff(thisTTL) < -0.5); % find instances where the TTL goes from high to low
        
        start_times 	= timestamps(start_inds); % find the timestamps of start events
        end_times    	= timestamps(end_inds); % find the timestamps of end events
        
        end_times(end_times < start_times(1))       = []; % discard potential initial end without start
        start_times(start_times > end_times(end))   = []; % discard potential final start without end
        
        switch a % this determines what the start and end timestamps should be assigned to: trial/episode, LED/opto stim or piezo/whisk stim.
            case episode_input_nr
                episode_starts  = start_times;
                episode_ends    = end_times;
            case LED_input_nr
                LED_starts      = start_times;
                LED_ends        = end_times;
            case piezo_input_nr
                piezo_starts    = start_times;
                piezo_ends      = end_times;
        end
    end
end

% if we're using a morse signal to detect start and end of the protocol,
% detect morse code for start and end of trials
if qmorse
    ep_lengths      = episode_ends - episode_starts;
    morse_lengths   = round(ep_lengths * 10);
    
    start_ind       = strfind(morse_lengths(:)',morse_start) + length(morse_start);
    end_ind         = strfind(morse_lengths(:)',morse_stop) - 1;
    
    episode_starts  = episode_starts(start_ind:end_ind);
    episode_ends    = episode_ends(start_ind:end_ind);
end

%% A lot of cleanup and repair from here

% determine median episode length
episode_times   = episode_ends - episode_starts;
episode_length  = round(median(episode_times),1);

% all episodes should have the same length; episodes with anomalous
% length are likely arduino startup floating voltage artefacts; get rid
% of anomalies
qepisode        = round(episode_times,1) == episode_length;

episode_starts  = episode_starts(qepisode);
episode_ends    = episode_ends(qepisode);

total_length    = round(median(diff(episode_starts)),1);

% Correct for odd TTL pulses
ep_start    = min(episode_starts(:));
ep_end      = max(episode_ends(:));

% ensure we only look for events happening during episodes;
% make sure all starts have matching ends and  vice versa
LED_starts(LED_starts < ep_start)               = [];
LED_ends(LED_ends < ep_start)                   = [];

piezo_starts(piezo_starts < ep_start)           = [];
piezo_ends(piezo_ends < ep_start)               = [];

LED_starts(LED_starts > ep_end)                 = [];
LED_ends(LED_ends > ep_end)                     = [];

piezo_starts(piezo_starts > ep_end)             = [];
piezo_ends(piezo_ends > ep_end)                 = [];

piezo_starts(piezo_starts > max(piezo_ends))    = [];
piezo_ends(piezo_ends < min(piezo_starts))      = [];

% more cleaning up; assuming similar trial durations, see if there are any
% exceptionally long or short trials at the start or end (can happen
% because of floating voltages on the arduino pins)
qepisode2       = round(diff(episode_starts),1) == total_length; % compare 
anomalind       = find(~qepisode2); % indices of anomalous lengths

anomalearly    	= anomalind(anomalind < (.5 * length(episode_starts))); % find which anomalous indices are happenig early on in the trial sequence
anomalate      	= anomalind(anomalind > (.5 * length(episode_starts))); % find which are happening late

anomalate       = anomalate + 1;

anomalous       = [anomalearly(:); anomalate(:)];

episode_starts(anomalous)   = []; % remove anomalous trial starts
episode_ends(anomalous)     = []; % remove anomalous trial ends

% 
allwhisks                   = piezo_starts;
allwhisk_ends               = piezo_ends;

%% work out the velocity (length) of a whisk, and the frequency:

% Find which whisking onsets are the first of a trial, and which onsets
% are the last of a trial
allwhisk_firstvect          = allwhisks(find(diff(allwhisks) > (episode_length/2))+1);
allwhisk_lastvect           = allwhisks(find(diff(allwhisks) > (episode_length/2)));

whisk_starts            	= [allwhisks(1); allwhisk_firstvect(:)];
whisk_lasts              	= [allwhisk_lastvect(:); allwhisks(end)];

whisk_ends               	= allwhisk_ends(find(diff(allwhisks) > (episode_length/2))+1);
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

allwhisks(allwhisks < min(episode_starts))      = [];
allwhisks(allwhisks > max(episode_ends))        = [];
whisk_freqs(allwhisks < min(episode_starts))    = [];
whisk_freqs(allwhisks > max(episode_ends))      = [];
whisk_lengths(allwhisks < min(episode_starts))  = [];
whisk_lengths(allwhisks > max(episode_ends))    = [];

LED_starts(LED_starts < min(episode_starts))    = [];
LED_starts(LED_starts > max(episode_ends))      = [];
LED_ends(LED_ends < min(episode_starts))        = [];
LED_ends(LED_ends > max(episode_ends))          = [];

%% fill in any missing trial start and end events

trial_interval      = median(diff(episode_starts));
qtrialgap           = [0; diff(episode_starts) > trial_interval + 1];

while any(qtrialgap)
    trial_ind       = find(qtrialgap,1);
    episode_starts  = [episode_starts(1:trial_ind-1); episode_starts(trial_ind-1)+trial_interval;  episode_starts(trial_ind:end)];
    qtrialgap       = [0; diff(episode_starts) > trial_interval + 1];
end

qtrialgap = [0; diff(episode_ends) > trial_interval + 1];

while any(qtrialgap)
    trial_ind       = find(qtrialgap,1);
    episode_ends    = [episode_ends(1:trial_ind-1); episode_ends(trial_ind-1)+trial_interval;  episode_ends(trial_ind:end)];
    qtrialgap       = [0; diff(episode_ends) > trial_interval + 1];
end

%% Check for trials without corresponding events and remove

% find the nr of episodes
ntrials         = length(episode_starts);

% Check for trials without corresponding events
for j = 1:length(episode_starts)
    this_start  = episode_starts(j);
    this_end    = episode_ends(j);
    
    qLEDstart  	= any(LED_starts > this_start & LED_starts < this_end);
    qwhisk      = any(allwhisks > this_start & allwhisks < this_end);
    
    if ~qLEDstart & ~qwhisk
        episode_starts(j)   = NaN;
        episode_ends(j)     = NaN;
        warning(['No LED or whisker stim found on trial ' num2str(j)])
    end
    if ~qLEDstart & qwhisk
        allwhisks(j)        = NaN;
        whisk_freqs(j)      = NaN;
        whisk_lengths(j)    = NaN;
        episode_starts(j)   = NaN;
        episode_ends(j)     = NaN;
        warning(['No LED found on trial ' num2str(j)])
    end
    if qLEDstart & ~qwhisk
        episode_starts(j)   = NaN;
        episode_ends(j)     = NaN;
        LED_starts(j)       = NaN;
        LED_ends(j)         = NaN;
        warning(['No whisker stim found on trial ' num2str(j)])
    end
end

allwhisks(isnan(allwhisks))             = [];
whisk_freqs(isnan(whisk_freqs))         = [];
whisk_lengths(isnan(whisk_lengths))     = [];
episode_starts(isnan(episode_starts))   = [];
episode_ends(isnan(episode_ends))       = [];
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

% recover LED delays
LED_delays      = round((LED_starts(:) - episode_starts(:)) / LED_conditions_res,3) * LED_conditions_res;

% recover whisking delays
whisk_delays    = round((allwhisks(:) - episode_starts(:)) / whisk_conditions_res,3) * whisk_conditions_res;

% recover LED durations
LED_ontimes     = LED_ends - LED_starts;
LED_durations   = round(LED_ontimes(:) / 10,3) * 10;

% reconstruct trial matrix
trial_timings   = [LED_delays(:) whisk_delays(:) LED_durations(:) whisk_freqs(:) round(1./whisk_lengths(:))];

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
    trace_data = [];
    for a = 1:n_channels
        disp(['Loading channel ' num2str(a)]);
        [thistrace timestamps info] = load_open_ephys_data([datafolder filesep filefolder filesep data_prefix '_CH' num2str(get_channels(a)) '.continuous']);
        
        starttime          	= min(timestamps);
        
        timestamps          = ((1:length(thistrace)) / 30000); % generate timestamps
        
        % Do simple homebrew spike detection by 1) filtering with bandpass, 2) smoothing and 3) thresholding
        thistrace           = filter(filt_b,filt_a,thistrace); % filter data with butterworth bandpass filter to get spike traces
        thistrace           = smooth(thistrace, spike_smoothwin); % smooth data with
        
        trace_data(a,1:length(thistrace))   = thistrace;
        
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


%% Sort spikes by episode and condition
for a = 1:n_channels
    
    chan_spike_times    = [];
    for b = 1:length(episode_starts)
        qspiketimes                 = (spikes(a).times >= episode_starts(b)) & (spikes(a).times < episode_ends(b));
        thesespiketimes             = spikes(a).times(qspiketimes);
        thesespiketimes             = thesespiketimes - episode_starts(b);
        
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
