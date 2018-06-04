
function oe2mda_multi(oe_data_loc,save_loc,channels,range)

% detect_threshold            = 5;
% channels                    = [1:4 13:20 29:32];
% save_loc                    = '/Users/Joram/Mountainsort_test/mda';
% file_loc                    = '/Users/Joram/Documents/In vivo data/2018_04_24/ChR2-YFP-ON_2018-04-24_13-32-29_1';

file_prefix               	= '100_CH'; % hardcoded in as this doesn't change in our data set
file_ext                    = '.continuous'; % as above

[filt_b,filt_a]             = butter(2, [300 6000]/(30000/2));

chan_data   = NaN(length(channels),length(range));
for a = 1:length(channels)
    this_channel                = channels(a);
    
    % convert channel nr to string
    chan_nr                     = num2str(this_channel);
    
    % make file name string of openephys .continuous file to be read
    read_file                   = [oe_data_loc filesep file_prefix chan_nr file_ext];
    
    % load openephys data file
    [data, timestamps, info]    = load_open_ephys_data(read_file);
    size(data(range))
    chan_data(a,:)              = data(range);% filter(filt_b,filt_a,data(range));
    size(chan_data)
end

% make file name string
write_file                  = [save_loc filesep 'All_16_chan.mda'];

% make directory for mda files if it does not yet exist
if ~isdir(save_loc)
    mkdir(save_loc)
end

%     % write file as mda (note that data has to go in as a row)
disp(['Writing mda file ' write_file])
disp(size(chan_data))
writemda16i(chan_data,write_file)

