function [out_values, cluster_nr] = contract_to_cluster_mean(in_values, tolerance, tolerance_mode);
% [out_values, cluster_nr]    = contract_to_cluster_mean(in_values, tolerance, tolerance_mode);
%
% /!\ WARNING /!\ - THIS IS A DUMB ALGORITHM CREATED FOR A VERY SPECIFIC PURPOSE.
% 
% TAKES A RANGE OF VALUES THAT SHOULD ALREADY BE TIGHTLY CLUSTERED AROUND A
% NUMBER OF SPECIFIC VALUES. OUTPUT WILL BECOME GARBAGE WHEN THERE IS NO
% EXTREMELY CLEAN SEPARATION BETWEEN CLUSTER RANGES.
% 
% WRITTEN MAINLY TO PULL ANALOG PULSEPAL D-TO-A OUTPUT VALUES TO A SENSIBLE 
% GROUPED MEAN VALUE TO RECOVER THE DIGITALLY REQUESTED VALUES
% 
% It finds an initial lower bound value, and then groups all values that fall
% within a range of TOLERANCE upwards from this lower bound together, and sets
% all those values to the mean of all values within this range.
% 
% OUTPUTS:
% OUT_VALUES: a vector the same length as in_values, in which the values have been
% replaced with the mean of the cluster of values that the value falls into.
% 
% CLUSTER_NR: a vector the same length as in_values, indicating in which cluster
% number each value fell (clusters numbered in ascending order of value)
%
% 
% IN_VALUES: a 1d vector of values that are tightly clustered around a number
% of mean values (e.g. the actual analog output values of a digital-to-analog
% converter, tightly clustered around the digitally requested values).
% 
% TOLERANCE: a tolerance (either as percentage or absolute, see below) around 
% each cluster. A full cluster should comfortably fall within this range.
% 
% TOLERANCE_MODE: 'percentage' or 'absolute'.
% 
% When 'percentage' is requested, the tolerance around a cluster is in a 
% percentage of the value of the cluster; e.g. a tolerance of 10% measn that
% for a cluster with a lower bound of 80 the algorighm will include values that
% fall within a range of 8 upwards from the lower bound; but for a cluster 
% around a value of 8 it will expect the values to fall within a range of only
% .8 upwards from the lower bound.
% 
% When 'absolute' is requested, the tolerance will not scale with values and
% will therefore be the same regardless of the value of the cluster.
% 
% 'absolute' is the default.
% 

%% Input checks

% Complain if no tolerance is provided
if nargin < 2
    error('Please specify tolerance')
end

% Default tolerance mode is 'absolute'
if nargin < 3
    tolerance_mode  = 'absolute';
end

% Set perc_mode flag depending on tolerance_mode, complain if input is invalid
if strcmpi(tolerance_mode, 'absolute')
    perc_mode = false;
elseif strcmpi(tolerance_mode,'percentage')
    perc_mode = true;
else
    error('unrecognised tolerance_mode type')
end

%% Main algorithm starts here

% Sort the input values
[sort_values, sort_inds]   	= sortrows(in_values);

% Recover the indices needed to invert the sorting / redistribute contracted values
% according to input value order later
[~, inverse_sort_indices] 	= sort(sort_inds);

% Preallocate array for cluster numbers
cluster_nr                  = zeros(size(sort_values));

% Some starting values for variables to go into the loop below:
this_cluster              	= 1; % Start with cluster 1
low_bound                   = sort_values(1); % set lower bound of first cluster to first (lowest) value
cluster_start_ind        	= 1; % Start at first value

% go through all values (sorted in ascending order)
for a = 1:length(sort_values)
    
    % Find current value
    this_level  = sort_values(a);
    
    % set current tolerance depending on tolerance mode
    if perc_mode
        this_tolerance  = (low_bound * tolerance) / 100; % set tolerance as a percentage of lower bound of current cluster
    else
        this_tolerance  = tolerance; % simply use absolute tolerance value
    end
    
    % If we have reached a value that exceeds (low_bound + this_tolerance)
    if this_level > (low_bound + this_tolerance)
        cluster_end_ind   	= a-1; % the
        
        % Set all cluster_nr values to the target cluster number
        cluster_nr(cluster_start_ind:cluster_end_ind) = this_cluster;
        
        % Now set new values for following cluster:
        low_bound           = sort_values(a); % current value becomes new lower bound
        this_cluster      	= this_cluster + 1; % cluster number increments by one
        cluster_start_ind 	= a; % the new cluster starts with the index of the current value
    end
    
end

% use inverse sort indices to rearrange cluster_nr so that its cluster_nrs
% correspond to the input values
cluster_nr  = cluster_nr(inverse_sort_indices);
uniq_bins   = unique(cluster_nr);

% pre-allocate output variable
out_values  = zeros(size(in_values));
for a = 1:length(uniq_bins)
    this_cluster        = uniq_bins(a);
    q_bin               = cluster_nr == this_cluster;
    
    out_values(q_bin)   = mean(in_values(q_bin));

end
