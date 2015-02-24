function [B95,Be,alphaaccept,p_upperbound] = Bneeded(alpha,targetu,varargin)
% Bootstrap samples needed (B)
% so that contribution to alpha estimation error of bootstrap is within utarget % of alpha 95% of the
% time
%
% B95 and Be depend on target alpha and target acceptable p-value variability (targetu) expressed as a percentage
% of alpha.  Because this computation reflects the contribution of the
% boostrapping procedure specifically, it does not reflect error in
% p-values due to the limited number of observations in the sample.
% Therefore, n is effectively infinite for the calculation of B95 and Be,
% and the bootstrap samples minimize the contribution of the bootstrapping
% itself to the error.  Other sources of error due to finite sampling may
% also increase uncertainty, as with any test.
%
% [B95,Be,p_upperbound] = Bneeded(alpha,targetu);
% [B95,Be,p_upperbound] = Bneeded(.01,.15);
%
% p_upperbound depends on variance of prctile at alpha, n in sample, and B
% number of bootstrap samples.  It reflects both sampling error and error
% due to the collection of a limited number of bootstrap samples.
% p_upperbound displays some strange nonlinear behavior with low alpha
% values and small to moderate samples. Thus, it may be 
%
% [B95,Be,p_upperbound] = Bneeded(alpha,targetu,n,B)
% [B95,Be,p_upperbound] = Bneeded(.01,.15,30,1000)
%
% Example
% [B95,Be,alphaaccept] = Bneeded(.001,.10);
% We need B95 boot samples to give us a p-value of alphaaccept or less 95% of
% the time with infinite sample size n. This is the error in the p-value due to bootstrapping.
% 
% [B95,Be,alphaaccept,p_upperbound] = Bneeded(.001,.10,Inf,B95)
% If we pass B95 in as the number of boot samples B and then get
% p_upperbound with n = Inf, we can reproduce alphaaccept
%
% Bottom line:
% with B95 samples and a true p-value of alpha, if we repeated the bootstrapping many times,
% we should get p-values less than or equal to alphaaccept 95% of the time.
% With Be samples, we should get p-values <= alphaaccept 84% of the time.
%
% tor wager, Feb 4, 2007

% q is prctile of distribution
q = norminv(alpha);

% g is pdf at q
g = normpdf(q);
g = max(g, .00001); % avoid zeros

a = ( alpha .* (1 - alpha) ) ./ g.^2;

alphaaccept = targetu .* alpha + alpha;
alphaaccept = min(alphaaccept, 1);  % can only accept fpr of 1 max.
% assume n is infinite, so var will come from B
% dividend changed from 1.96 to 1.64, cause this is two tailed...
b = ( (norminv(alphaaccept) - q) ./ 1.6450 ) .^ 2;

B95 = a ./ b;

b = ( (norminv(alphaaccept) - q) ./ 1 ) .^ 2;
Be = a ./ b;

% 19.7 in bootstrapping book
% as n and B approach Inf, approaches 0
%varq = (alpha .* (1 - alpha) ) ./ g .^ 2 .* (1 ./ n + 1 ./ B);

if nargout > 3
    if length(varargin) < 2, error('Incorrect usage: not enough inputs to compute all requested outputs.  See help notes.'); end
    n = varargin{1};
    B = varargin{2};
    
    %a = ( alpha .* (1 - alpha) ) ./ g.^2;
    varq = a .* (1 ./ n + 1 ./ B);
    
    p_upperbound = min(1,normcdf(q + 1.96 * sqrt(varq)));
    
end

end