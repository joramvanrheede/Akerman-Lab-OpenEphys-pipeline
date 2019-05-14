function fit_model = fit_double_exp(in_vals,t_vals, sign)
% function fit_model = fit_double_exp(in_vals,t_vals)
% 
% Fits a double exponential to in_vals; useful for postsynaptic currents
% or for fitting the profile of a post-stimulus time histogram.
% 
% It uses the following formula:
% y(t) = a * exp(-exp(-((t-b)/c))-((t-b)/c)+1.0)
% 
% 'a' determines the height of the peak and is initialised at the max value
% of IN_VALS
% 
% 'b' determines the horizontal offset and is initialised at the time of the max
% value of IN_VALS
% 
% 'c' determines the rise time of the function
% 
% 'd' determines the decay
% 
% IN_VALS is a vector with y-values you want to fit, e.g. current or spike 
% count.
% 
% T_VALS should be a vector the same length as IN_VALS with the timestamps 
% corresponding to the y-values.
% 
% SIGN: if sign is -1 (or any negative number), the function will fit to a 
% current with a negative peak. If not provided or positive, the function
% will default to a current with a positive peak.
%

if nargin < 3 || sign > 0 % default, or positive fit direction requested
    % Find max and peak time of the data to initialise parameters
    [psth_peak,peak_ind]    = max(in_vals);
    peak_time               = t_vals(peak_ind);
else % negative fit direction requested
    % Find min and trough time of the data to initialise parameters
    [psth_peak,peak_ind]    = min(in_vals);
    peak_time               = t_vals(peak_ind);
end

% Find max and peak time of the data to initialise parameters
[psth_peak,peak_ind]    = max(in_vals);
peak_time               = t_vals(peak_ind);

% The double exponential function that is used for the fit:
double_exp      = 'a * exp(-exp(-((x-b)/c))-((x-b)/d)+1.0)';

% Create a fit object according to this function
fit_obj     	= fittype(double_exp);

% Perform the fit
fit_model       = fit(t_vals(:), in_vals(:), fit_obj,'Start',[psth_peak peak_time peak_time peak_time],'Weight',length(in_vals)+1-(1:length(in_vals))); %'Start',[psth_peak 0.5 peak_time],'Weight',fliplr(1:length(testpsth)))









