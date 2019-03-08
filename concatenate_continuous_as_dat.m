function concatenate_continuous_as_dat(filename,data_dirs,n_channels)
% function concatenate_continuous_as_dat(FILENAME,DATA_DIRS,N_CHANNELS)
% 
% Takes openephys .continuous data folders specified in DATA_DIRS and
% concatenates them into a single .dat file with filename specified in
% FILENAME.
% 
% FILENAME      = full file name (including path and '.dat' extension) of the .dat file to write to
% DATA_DIRS 	= cell array of openephys data directories to be concatenated
% N_CHANNELS  	= number of channels of openephys .continuous data
% 
% Adapted from Kilosort's convertOpenEphysToRawBInary by Joram van Rheede - March 2019

% How many openephys data sets to concatenate?
n_data_sets     = length(data_dirs);

% Open new .dat file for writing
fidout          = fopen(filename, 'w');

% loop over number of openephys data sets to be concatenated
for a = 1:n_data_sets
    
    disp(['Processing data folder ' data_dirs{a} '...'])
    
    % Get full filenames of .continuous files for each channel
    clear fs
    for j = 1:n_channels
        fs{j} = dir(fullfile(data_dirs{a}, sprintf('*CH%d.continuous', j) ));
    end
    
    nblocks = cellfun(@(x) numel(x), fs);
    if numel(unique(nblocks))>1
        error('different number of blocks for different channels!')
    end
    
    %
    nBlocks     = unique(nblocks);
    nSamples    = 1024;  % fixed to 1024 for now!
    
    fid = cell(n_channels, 1); % one file identifier per channel
    
    for k = 1:nBlocks
        
        for j = 1:n_channels
            % open .continuous file for this channel
            fid{j}             = fopen(fullfile(data_dirs{a}, fs{j}(k).name));
            % skip header information
            fseek(fid{j}, 1024, 0);
        end
        

        flag = 1;
        while flag
            samples = zeros(nSamples * 1000, n_channels, 'int16');
            for j = 1:n_channels % loop over the data files for individual channels
                collectSamps    = zeros(nSamples * 1000, 1, 'int16'); % preallocate space for the data
                
                rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b'); % reading in data from the jth channel
                
                nbatches        = ceil(numel(rawData)/(nSamples+6)); % Doing this in 'batches'? What's up with the recurring 6 anyway? Looks like we are skipping 6 samples after every batch, is this why this needs to happen in batches?
                for s = 1:nbatches
                    rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]); % s == 1 --> 0 + 6 + 1:nsamples = 7:nsamples+6; s == 2 --> 1 * (nsamples + 6) + 6
                    collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps; % s == 1 --> 0 + 1:nsamples
                end
                samples(:,j)         = collectSamps;
            end
            
            if nbatches<1000 % presumably we have reached the end of the file
                flag = 0;
                samples = samples(1:s*nSamples, :);
            end
            
            samples         = samples'; % 
            
            fwrite(fidout, samples, 'int16'); % Writing to .dat file happens here
            
        end
        
        % Close the open .continuous files
        for j = 1:n_channels
            fclose(fid{j});
        end
        
    end
end

fclose(fidout); % close .dat file here
