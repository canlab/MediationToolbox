function v = v_from_targetu(targetu,alpha)
% v = v_from_targetu(targetu,alpha)
%
% targetu is alpha inflation factor (% above alpha tolerable)
%
% q is prctile of distribution
q = norminv(alpha);

alphaaccept = targetu .* alpha + alpha;

alphaaccept = min(alphaaccept, 1);  % can only accept fpr of 1 max.

v = ( (norminv(alphaaccept) - q) ./ 1.96 ) .^ 2;

end



