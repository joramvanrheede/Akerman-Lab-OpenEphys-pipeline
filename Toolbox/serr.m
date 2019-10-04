function stderror = serr(input, dim)

% function stderror = serr(input,dim)
% 
% Quick and dirty standard error function.
% Ignores NaNs.
% Works on vectors or on first non-singleton dimension of input unless 'dim'
% is specified.
% If 'dim' is 'all', will operate on all elements of input.

if nargin <2
    % find standard deviation of input
    instd       = nanstd(input);
    % determine n (non-NaN)
    n           = sum(~isnan(input));
else
    % find standard deviation of input
    instd       = nanstd(input,[],dim);
    % determine n (non-NaN)
    n           = sum(~isnan(input),dim);
end



% divide by square root of n
stderror    = instd ./ sqrt(n);

