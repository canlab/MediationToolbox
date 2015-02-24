function [p, p2, z_lowertail, p_2tail] = bootbca_pval(testvalue,bootfun, bstat, stat, varargin)
    % [p, p2, z, p_2tail] = bootbca_pval(testvalue,bootfun, bstat, stat, [x], [other inputs to bootfun])
    % 
% p-values from corrected and accelerated percentile bootstrap
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
% T.J. DiCiccio and B. Efron (1996), "Bootstrap Confidence Intervals",
% statistical science, 11(3)
%
% Adapted to return p-values rather than CIs for massive univariate analysis by Martin Lindquist
% implemented by Tor Wager
% P-value reflects prob. that stat is ~= testvalue (e.g. zero, for 1-sample
% t-test against Ho: stat = 0)
%
%
% stat = statistic calculated from original sample
% bstat = B x 1 vector of statistics calculated from bootstrap sample
% i.e.,  
% % stat = bootfun(varargin{:});
% % bstat = bootstrp(nboot,bootfun,varargin{:});
%
% Examples:
% x = randn(100,1) + .1;
% stat = mean(x);
% bootstat = bootstrp(1000, @mean, x);
% [p, z] = bootbca_pval(0, @mean, bootstat,stat, x);
%
% wmean = @(w,paths) diag(w'*paths)';
% [p, z] = bootbca_pval(0,wmean, [], [], x, wts);
%
% Simulated test data:
% % for j = 1:10
%     rand('twister',sum(100*clock))
% 
%     niter = 1000;
%     n = 20;
%     
%     stat = zeros(niter, 1); p = stat; p2 = stat; z = stat; p_2tail = stat;
%     for i = 1:niter, x = randn(n,1); stat(i) = mean(x); bstat = bootstrp(1000,@mean,x);
%         [p(i), p2(i), z(i), p_2tail(i)] = bootbca_pval_onetail(0,@mean,bstat,stat(i),x);
%         if mod(i,100) == 0,fprintf(1,'%03d ',i); end
%     end
%     fpr_pos(j) = sum(p < .05) ./ niter;
%     fpr_neg(j) = sum(p2 < .05) ./ niter;
%     fpr_2tail(j) = sum(p_2tail < .05) ./ niter;
% 
% end
% 
% fpr_pos
% fpr_neg
% fpr_2tail
% 
% mean(fpr_pos') ./ ste(fpr_pos')
% mean(fpr_neg') ./ ste(fpr_neg')
% mean(fpr_2tail') ./ ste(fpr_2tail')
%
% create_figure('mean_vs_Z'); plot(stat, z, 'k.');
% create_figure('pvals'); hist(p, 100);
% create_figure('zscores'); hist(z, 100);


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
a =  skew ./ 6;  % acceleration

% prctile of the distribution below nullvalue (2-tailed)
% pct = min( [sum(bstat < testvalue); sum(bstat > testvalue)] ) ./ B;
pct_lowertail = sum(bstat < testvalue) ./ B;
pct_uppertail = sum(bstat > testvalue) ./ B;

% lower tail is smaller for positive effect. is_lowertail == 1 effects
% should be positive.
is_lowertail = 1; %pct_lowertail < pct_uppertail;

% make sure pct is not 1 or 0 due to limited bootstrap samples
pct_lowertail = max(pct_lowertail, 1 ./ B);
pct_uppertail = min(pct_uppertail, 1 - 1./B);

pct_uppertail = max(pct_uppertail, 1 ./ B);
pct_lowertail = min(pct_lowertail, 1 - 1./B);
 
% % pct = min(sum([bstat < testvalue bstat > testvalue])) ./ B; %
% non-matrix form

% bias-corrected z-score of 
zpct_lowertail = norminv(pct_lowertail) - z_0;  % neg z for pos effect, and vice versa
zpct_uppertail = norminv(pct_uppertail) - z_0;

% acceleration adjustment
z_lowertail = ( zpct_lowertail .* (1 - a .* (z_0)) - z_0 ) ./ (1 + a .* zpct_lowertail);
z_uppertail = ( zpct_uppertail .* (1 - a .* (z_0)) - z_0 ) ./ (1 + a .* zpct_uppertail);

%z = z_uppertail;
%z(is_lowertail) = -z_lowertail(is_lowertail);

%p = normcdf(z); 
%p = normcdf(-z_lowertail) + (1 - normcdf(z_uppertail)); 

p = normcdf(z_lowertail);  % p-value for whether mean > test value  (pos effect)

p2 = normcdf(z_uppertail); % p-value for whether mean < test value (neg effect)

% reverse sign to be consistent with typical direction of z-score
z_lowertail = -z_lowertail;     
z_uppertail = -z_uppertail;

% same
% p = min([p; 1-p]);
% p2 = min([p2; 1-p2]);
% 
% p_2tail = p + p2;

p_2tail = min(p, 1 - p) + min(p2, 1 - p2);

%p = 2 .* min([p 1-p]);

% % % 
% % % p = normcdf(z_uppertail);
% % % p(is_lowertail) = normcdf(z_lowertail(is_lowertail));


% % z = [z_lowertail; z_uppertail];
% % 
% % p = [normcdf(z_lowertail); (normcdf(z_uppertail))];
% % 
% % % take min of lower and upper p-values
% % % if lower pval is smaller, + z-score, else - z-score
% % 
% % whichvals = p == repmat(min(p), 2, 1);
% % 
% % p = p(whichvals)';
% % z = z(whichvals)';


%fprintf(1,'Pcrtile p = %3.4f, Adj. p = %3.4f, diff = %3.4f\n',2*pct, p, 2*pct-p)
doplot = 0;
if doplot
    
create_figure('bstat'); hist(bstat, 100);
handle = plot_vertical_line(testvalue); set(handle, 'LineWidth', 2);
text(testvalue + .1, 20, sprintf('Raw: %3.0f%%', 100*pct_uppertail), 'FontSize', 24)
text(testvalue - .2, 20, sprintf('Raw: %3.0f%%', 100*pct_lowertail), 'FontSize', 24)

text(testvalue - .2, 15, sprintf('BCA: %3.0f%%', 100*normcdf(z_lowertail)), 'FontSize', 24, 'Color', 'r')
text(testvalue + .1, 15, sprintf('BCA: %3.0f%%', 100*normcdf(z_uppertail)), 'FontSize', 24, 'Color', 'r')

fprintf('Positive effect: z = %3.2f, p = %3.6f\n', z_lowertail, p);
fprintf('Negative effect: z = %3.2f, p = %3.6f\n', z_uppertail, p2);

end

end
