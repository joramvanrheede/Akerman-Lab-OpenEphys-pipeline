
function oe2mda(oe_data_loc,save_loc,this_channel)

% detect_threshold            = 5;
% channels                    = [1:4 13:20 29:32];
% save_loc                    = '/Users/Joram/Mountainsort_test/mda';
% file_loc                    = '/Users/Joram/Documents/In vivo data/2018_04_24/ChR2-YFP-ON_2018-04-24_13-32-29_1';

file_prefix               	= '100_CH'; % hardcoded in as this doesn't change in our data set
file_ext                    = '.continuous'; % as above

% convert channel nr to string
chan_nr                     = num2str(this_channel);

% make file name string of openephys .continuous file to be read
read_file                   = [oe_data_loc filesep file_prefix chan_nr file_ext];

% load openephys data file
[data, timestamps, info]    = load_open_ephys_data(read_file);

% make file name string
write_file                  = [save_loc filesep file_prefix chan_nr '.mda'];

% make directory for mda files if it does not yet exist
if ~isdir(save_loc)
    mkdir(save_loc)
end

%     % write file as mda (note that data has to go in as a row)
disp(['Writing mda file ' write_file])

writemda(data',write_file)

