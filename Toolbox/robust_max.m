function max_val    = robust_max(indata, percentile, dim)
% function max_val    = robust_max(indata, percentile, dim)
% 
% 
% Get the maximum of INDATA bar the top PERCENTILE % of data points.
% Useful for e.g. setting colormap ranges that are not overly affected by 
% individual outliers.
% 
% Operates on N-D arrays along dimension specified by DIM
% Alternatively, if DIM = 'all' it will give a single max value for all
% data points.


% if DIM == 'all', function will treat INDATA as a single vector

if nargin < 2
    percentile = 5;
end

if nargin < 3
    dim = 1;
end

if dim == 'all'
    indata = indata(:);
    dim = 1;
end

% computes the max along 
sorted_indata   = sort(indata,dim,'descend');

sort_dim_size   = size(sorted_indata,dim);

cutoff          = sort_dim_size / 100 * percentile;
cutoff          = ceil(cutoff); % Round up so minimum is 1

n_dims          = ndims(sorted_indata);

indx    = [];
for a = 1:n_dims
    if a == dim
        indx{a}     = (cutoff+1):size(sorted_indata,dim);
    else
        indx{a}     = ':';
    end
    
end

trimmed_indata  = sorted_indata(indx{:}); % 

max_val         = max(trimmed_indata,[],dim);


