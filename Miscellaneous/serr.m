function stderror = serr(input)

% function stderror = serr(input)
% 
% Quick and dirty standard error function.
% Ignores NaNs.
% Works on vectors or on first non-singleton dimension of input.

% find standard deviation of input
instd       = nanstd(input);

% determine n (non-NaN)
n           = sum(~isnan(input));

% divide by square root of n
stderror    = instd ./ sqrt(n);
