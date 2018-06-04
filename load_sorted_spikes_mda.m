% load_sorted_spikes 

mda_spikes_file = '/Users/Joram/Mountainsort_test/data/firings30.mda.prv';


%% use prv file to locate original file
fid             = fopen(mda_spikes_file);
prv_string      = fscanf(fid,'%s');
fclose(fid)

path_preamble	= 'original_path":"'; 
path_postamble  = '.mda';

preamb_length   = length(path_preamble);
postamb_length  = length(path_postamble);

preamb_start    = strfind(prv_string,path_preamble);
filename_start  = preamb_start + preamb_length;

postamb_start   = strfind(prv_string,path_postamble);
filename_end    = postamb_start + postamb_length - 1;

mda_spikes_file = prv_string(filename_start:filename_end);

%%

spike_data      = readmda(mda_spikes_file);
spike_data(2,:) = spike_data(2,:) / 30000;



