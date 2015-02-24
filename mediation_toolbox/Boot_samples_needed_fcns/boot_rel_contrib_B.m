function [relB,propB,alphainflate,var_caused_by_B] = boot_rel_contrib_B(alpha,n,B)

% 19.7 in bootstrapping book
% as n and B approach Inf, approaches 0
%varq = (alpha .* (1 - alpha) ) ./ g .^ 2 .* (1 ./ n + 1 ./ B);

% q is prctile of distribution
q = norminv(alpha);

% g is pdf at q
g = normpdf(q);

c = (alpha .* (1 - alpha) ) ./ g .^ 2;

varq = c .* (1 ./ n + 1 ./ B);

relB = (c./n) ./ (c./B);  % = B / n

propB = (c ./ B) ./ varq; % = n / (n + B)

    CIalpha = [normcdf(q - 1.96 .* varq .^.5); normcdf(q + 1.96 .* varq .^.5)];
    wid = CIalpha(2,:) - CIalpha(1,:);

varqideal = c .* (1 ./ n);
CIideal = [normcdf(q - 1.96 .* varqideal .^.5); normcdf(q + 1.96 .* varqideal .^.5)];
widideal = CIideal(2) - CIideal(1);

alphainflate = (wid - widideal) ./ widideal;

var_caused_by_B = (c ./ B);

end