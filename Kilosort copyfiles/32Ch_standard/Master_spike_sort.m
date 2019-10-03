% Run config file
run('StandardConfig.m')

tic; % start timer

if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
rezToPhy(rez, cd);

% remove temporary file
delete(ops.fproc);
