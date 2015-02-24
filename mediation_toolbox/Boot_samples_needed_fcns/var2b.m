function [B,min_n,Blowerlimit] = var2b(v,n,alpha)
% Bootstrap samples needed (B)
% for a target prctile of variance (v)
% and n and alpha

% q is prctile of distribution
q = norminv(alpha);

% g is pdf at q
g = normpdf(q);

bpart = v .* g.^2 ./ ( alpha .* (1 - alpha) );
B = 1 ./ ( bpart - (1 ./ n) );

% 19.7 in bootstrapping book
% as n and B approach Inf, approaches 0
%varq = (alpha .* (1 - alpha) ) ./ g .^ 2 .* (1 ./ n + 1 ./ B);

varlowerlimit = (alpha .* (1 - alpha) ) ./ g .^ 2 .* (1 ./ n + 0);

% variance can never be this low for this n
B(v <= varlowerlimit) = Inf;

min_n = ( (alpha .* (1 - alpha) ) ./ g .^ 2 ) ./ v;

% as n approaches infinity
Blowerlimit = 1 ./ ( bpart - (1 ./ Inf) );

end



