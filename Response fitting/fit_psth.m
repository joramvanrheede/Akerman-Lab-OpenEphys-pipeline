% fit_psth
% A script that will fit a double polynomial to a PSTH from ephys_data.
% Controlled by a number of input parameters at the top of the script.
% The double exponential fit really only works when you start the fit from 
% the stimulus time or even shortly after (curr_offset) 

%% Input parameters
file_name   = '/Volumes/Akermanlab/Joram/Preprocessed data/AVK RBSN rAAV PoM/Drive/2019_04_01/2019_04_01-15-Drive.mat'; % full path to preprocessed ephys_data file
cond_nr     = [2]; % Which condition to look at within this data file 

fit_win     = [0 .3]; % window over which to fit the function (relative to stimulus onset). Start at 0  is recommended.
bin_size    = [0.001]; % Post-stimulus time histogram bin size
channels    = [1:32]; % Which channels to include

smooth_win  = 1; % smoothing window (1 = no smoothing); increasing may help get a sensible fit if the PSTH is too messy / sparse for fit to work
curr_offset = 0.000; % Offset for the fit; use 0 unless the current appears delayed and it causes problems with the fit

%% Code execution starts here

load(file_name)

% Get the onset of the 'current' / the shape of interest in the PSTH
curr_onset      = ephys_data.conditions(cond_nr).whisk_onset + curr_offset;
spikes          = ephys_data.conditions(cond_nr).spikes(channels,:,:);

% Make PSTH and get spike counts and time bin edges
figure
[handle, psth_counts, edges] = psth(spikes-curr_onset,bin_size);

% select only the part of the PSTH that is needed for the fit
q_fit           = edges > fit_win(1) & edges <= fit_win(2);
edges           = edges(q_fit);
psth_counts     = psth_counts(q_fit);

% apply smoothing; this can help for very sparse / messy PSTHs
smooth_psth     = smooth(psth_counts,smooth_win);

% Function that handles the fitting
fit_model       = fit_double_exp(smooth_psth,edges);

% Plotting for visualisation of fit, on top of PSTH:
hold on
plot(edges, fit_model(edges),'g:','LineWidth',4)
fixplot
xlim(fit_win)

% Uncomment the following to look at the difference between the fit and the psth:
% plot(edges, smooth_psth - fit_model(edges),'r:','LineWidth',4)

