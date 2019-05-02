function ephys_data = get_ephys_data(datafolder, events_chans, channel_map, spike_thresh)


% Load ADC channels
sync_data       = get_stim_sync_data(datafolder, events_chans);



%% Some code to find file prefix for data and sync files

data_files      = dir(fullfile(datafolder,'*.continuous'));
data_files      = {data_files.name}';

is_adc_file     = ~cellfun(@isempty,regexp(data_files,'ADC')); % Remove ADC files

chan_files      = data_files(~is_adc_file); %
chan_inds       = cell2mat(regexp(data_files,'CH')); % find index of occurrence of 'CH'

chan_ind        = chan_inds(1) - 1; % We only need to look at the prefix for one file
chan_file       = chan_files{1}; %
data_prefix     = chan_file(1:chan_ind);

n_channels      = length(channel_map);

%% Filter settings
samplefreq = 30000;

[filt_b,filt_a]       	= butter(2, [500 6000]/(samplefreq/2));
[LFPfilt_b, LFPfilt_a] 	= butter(2, [1 300]/(samplefreq/2));

LFPtraces            	= [];
LFPtimestamps         	= [];

for a = 1:n_channels
    disp(['Loading channel ' num2str(a)]);
    disp(['File ' datafolder filesep data_prefix 'CH' num2str(channel_map(a)) '.continuous'])
    
    [thistrace timestamps info] = load_open_ephys_data([datafolder filesep data_prefix 'CH' num2str(channel_map(a)) '.continuous']);
    
    starttime          	= min(timestamps); % get original offset
    
    timestamps          = ((1:length(thistrace)) / 30000); % generate timestamps (openephys output t)
    
    % Do simple homebrew spike detection by 1) filtering with bandpass, 2) smoothing and 3) thresholding
    spiketrace        	= filter(filt_b,filt_a,thistrace); % filter data with butterworth bandpass filter to get spike traces
    
    q_threshold         = (-spiketrace) > (spike_thresh * std(spiketrace)); % determine standard deviation to determine threshold, detect threshold crossings (negative)
    spike_bool          = diff(q_threshold) == 1; % Determine instances of threshold being crossed
    these_spike_times 	= timestamps(spike_bool); % Get the timestamps of these instances
    
    spikes(a).times     = these_spike_times + starttime; % put original offset back
    
    %% Consider adding spike amplitude measure?
    
    
    LFPtrace        = filter(LFPfilt_b,LFPfilt_a,thistrace);
    LFPtrace        = LFPtrace(1:30:end); % resample at 1000Hz
    LFPtimestamps   = timestamps(1:30:end) + starttime; % resample to 1000Hz
    
    LFPtraces       = [LFPtraces LFPtrace(:)];
    
end


%% Put it all together

ephys_data  = sync_data;

for a = 1:length(ephys_data)
    original_rec_start_time     = ephys_data(a).rec_start_time;
    for b = 1:length(ephys_data(a).conditions)
        trial_starts    = ephys_data(a).conditions(b).trial_starts;
        trial_ends      = ephys_data(a).conditions(b).trial_ends;
        for c = 1:length(trial_starts)
            for d = 1:n_channels
                
                qspiketimes     = spikes(d).times >= trial_starts(c) & spikes(d).times < trial_ends(c);
                thesespiketimes = spikes(d).times(qspiketimes);
                thesespiketimes = thesespiketimes - trial_starts(c);
                
                %
                %                 % Make relative to Kilosort concatenated data time, not openephys recording time stamp
                %                 trial_start     = trial_starts(c) - original_rec_start_time;
                %                 trial_end       = trial_ends(c) - original_rec_start_time;
                %
                %                 q_trial         = spike_times >= trial_start & spike_times < trial_end;
                %
                %                 % get relevant spike times and make them relative to trial onset
                %                 channel_trial_spikes   = spike_times(q_channel & q_trial) - trial_start;
                %
                
                ephys_data(a).conditions(b).spikes(d,c,1:length(thesespiketimes))  = thesespiketimes;
                
                
                
                % LFP processing. Is it better to simply do this by condition, and
                % then have a matrix of data x episodes x data length?
                
                qLFPtimestamps              = LFPtimestamps >= trial_starts(c) & LFPtimestamps < trial_ends(c);
                ephys_data.conditions(b).LFP_trace(d,c,1:length(LFPtraces(qLFPtimestamps,d))) = LFPtraces(qLFPtimestamps,d);
                ephys_data.LFP_timestamps   = [1:length(qLFPtimestamps)] / 30000;
            end
        end
        
        % set empty values to NaN instead of 0
        ephys_data(a).conditions(b).spikes(ephys_data(a).conditions(b).spikes == 0) = NaN;
    end
end

