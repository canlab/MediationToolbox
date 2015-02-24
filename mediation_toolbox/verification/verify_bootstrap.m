% Compare our bootstrap CI to the Matlab built-in version
% ----------------------------------------------------------------------
%
% 2/15/10 : Martin and Tor : mediation bootstrap CI looks good

nboot = 1000;

for j = 1:200
    
ci = bootci(nboot, bootfun, y);

stat = bootfun(y);
bstat = bootstrp(nboot,bootfun,y);
ci2 = bootbca_ci(.025,bootfun, bstat, stat, y);

diffs(:, j) = ci - ci2';

end

[h, p, ci, stat] = ttest(diffs')

create_figure('test', 1, 2); hist(diffs(1, :)', 40); title('Lower bound'); subplot(1, 2, 2); hist(diffs(2, :)', 40), title('Upper bound');

%%
