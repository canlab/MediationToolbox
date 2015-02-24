% x and m share a feature in the data, and y and m share one
% x and m and y and m each have positive cov
% but there are really two different effects going on
% in this hypothetical example.
% the observed a and b are pos, and c is zero
% the mediation algorithm finds a negative direct link
% and a positive mediated link because of the pos. implied ab effect
% suggesting 'suppression', but this is just an artifact
%
% the 'implied correlation' between x and y is spurious, because it comes
% from two different features in the data (unmodeled causes)
% is allowing 'suppression' causing this instability, where there are
% two solutions-one with zero xy cov, and another with positive xy cov that is 'suppressed'?
% is there another model that could capture whether the a and b covariances
% are driven by the same features in the data?

f = ones(10, 1); %spm_hrf(1);
z = zeros(size(f));
%z = z .* mean(f);
x = [f;z;z;z];
m = [f;z;f;z];
y = [z;z;f;z];

x = x + rand(size(m)) .* .5;
y = y + rand(size(m)) .* .5;
m = m + rand(size(m)) .* .5;


corrcoef([x y m])
create_figure('Data plot', 3, 1); 
plot(x,'o-', 'MarkerFaceColor', [.5 .5 .5]); 
title('X variable')
subplot(3,1,2); plot(m, 'o-', 'MarkerFaceColor', [.5 .5 .5]); 
title('M variable')
subplot(3,1,3);plot(y,'o-', 'MarkerFaceColor', [.5 .5 .5]), 
title('Y variable')
xlabel('Observations')

[paths, stats] = mediation(x,y,m, 'stats', 'boot', 'verbose', 'plots');

% % disp(['Correlation in residuals'])
% % corrcoef([stats.residm{1} stats.residy{1} stats.residy2{1}])
% % 
% % tor_fig(3,1); plot(stats.residm{1}); subplot(3,1,2); plot(stats.residy{1},'r'); subplot(3,1,3);plot(stats.residy2{1},'g')

xmx = scale(x) .* scale(m);
myx = scale(m) .* scale(y);

[r, p] = corrcoef(xmx, myx);
fprintf(1,'Correlation of x-m crossproducts and m-y crossproducts\n')
fprintf(1,'This tests whether the units that drive x-m covariance are the same units as those that drive m-y covariance \n')
r, p

dat = [x m y];
plot_parallel_coords(dat)
set(gca, 'XTick', 1:3, 'XTickLabel', {'X' 'M' 'Y'})
ylabel('Value')

