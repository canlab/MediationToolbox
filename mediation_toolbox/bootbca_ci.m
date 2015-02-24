function ci = bootbca_ci(alph,bootfun, bstat, stat, varargin)
    % ci = bootbca_ci(alph,bootfun, bstat, stat, [x], [other inputs to bootfun])
    % 
% conf.interval from corrected and accelerated percentile bootstrap
% Specific implementation tested for mean, works on each column of x 
% (May not work well/not tested for other bootfun choices).
%
% Usage:
% % pass in finished bootstrap:
% [p, z] = bootbca_pval(testvalue,bootfun, bstat, stat, x, [other inputs]) 
% OR do bootstrap here
% [p, z] = bootbca_pval(testvalue,bootfun, [], [], x)  
%
% Based on BCa confidence intervals in:
% Efron & Tibshirani 1994
% T.J. DiCiccio and B. Efron (1996), "Bootstrap Confidence Intervals",
% statistical science, 11(3)
%
% implemented by Martin Lindquist & Tor Wager
%
% stat = statistic calculated from original sample
% bstat = B x 1 vector of statistics calculated from bootstrap sample
% i.e.,  
% % stat = bootfun(varargin{:});
% % bstat = bootstrp(nboot,bootfun,varargin{:});
%
% Examples:
% cis = zeros(1000,2); haszero = zeros(1000,1);
% n = 100;
% for i = 1:1000
% x = randn(n,1);
% ci = bootbca_ci(.025,@mean, [], [], x);
% cis(i,:) = ci;
% haszero(i) = ci(1) < 0 & ci(2) >= 0;
% if mod(i,50) == 0,fprintf(1,'%03d ',i);, end
% end
% fpr = 1 - (sum(haszero) ./ 1000)
%
% cis = zeros(1000,2); haszero = zeros(1000,1);
% n = 200;
% for i = 1:1000
% x = randn(n,1);
% ci = bootci(1000,{@mean, x},'alpha',.05)';
% cis(i,:) = ci;
% haszero(i) = ci(1) < 0 & ci(2) >= 0;
% if mod(i,50) == 0,fprintf(1,'%03d ',i);, end
% end
% fpr = 1 - (sum(haszero) ./ 1000)


if isempty(bstat) || isempty(stat)
    stat = bootfun(varargin{:});
    bstat = bootstrp(1000,bootfun,varargin{:});
end

% Bias correction
[B, ncols] = size(bstat);  

% if mean, proportion of bootstrap means that are less than the sample mean
prop_less = sum(bstat < repmat(stat, B, 1)) ./ B;
z_0 = norminv(prop_less);
 
% acceleration value, see DiCiccio and Efron (1996)

jstat = jackknife(bootfun, varargin{:});
n = size(jstat,1);

score = -(jstat - repmat(mean(jstat), n, 1) ); % score function at stat;
skew = sum(score.^3)./(sum(score.^2).^1.5);  %skewness of the score function
a =  skew ./ 6;  % acceleration  % Yoni:  i THINK that this is the correct acceleration for each bstat, individually. so below i only use the correct one, not all of the a's.


zalpha1 = norminv(alph);
zalpha2 = norminv(1 - alph);


a1 = z_0 + ( (z_0 + zalpha1) ./ (1 - a .* (z_0 + zalpha1)) );
a1 = normcdf(a1);

a2 = z_0 + ( (z_0 + zalpha2) ./ (1 - a .* (z_0 + zalpha2)) );
a2 = normcdf(a2);

ci_lower = prctile(bstat, 100 .* a1);
ci_upper = prctile(bstat, 100 .* a2);

ci = [ci_lower ci_upper];

% below is added by Yoni 4/22, i am "overwriting" the ci value above, which was the previous value returned.
% here, i only use use the acceleration (a) which I believe to be correct
% for given parameter/variable, tho I am not 100% sure of correctness.
clear ci
for i=1:size(ci_lower, 1)
    ci(i,1) = ci_lower(i,i);
    ci(i,2) = ci_upper(i,i);
end

end
