function[p,h] = Chi2_test(n1,N1,n2,N2,alpha);
%% Performs a simple Chi2 test on two proportions n1 events in N1 total observations, and n2 events in N2 observations. 
% resturns h  = 0 if null hypothesis true and both proportions from the
% same distribution. h= 1 if null hypothesis rejected and p value < set
% alpha.
% 
     % Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2);
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p = 1 - chi2cdf(chi2stat,1);
          if p < alpha;
              h = 1;
          else
              h =0;
          end;
end