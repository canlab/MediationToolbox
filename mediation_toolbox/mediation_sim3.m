function [pvals, results_table, stats, randseed] = mediation_sim3(n, t, s_g, s_a, s_b, apop, bpop, abcov, l2m_acov, l2m_bcov, color, varargin)
% Simulated results and significance rates (power or false positives) for multilevel mediation.
% generates a dataset (in an encapusulated subfunction) and runs multiple iterations, 
% collecting stats on significance (power or false positive rate, depending on whether there 
% is a true effect or not for each parameter) across them. Overall summary
% stats for each effect are in results_table.
%
% mediation_sim3(n, t, s_g, s_a, s_b, apop, bpop, abcov, l2m_acov, l2m_bcov, color, [optional inputs])
%
% Uses embedded subfunction:
% [x, m, y, l2m] = mediation_gen_data1(n, t, s_g, s_a, s_b, apop, bpop, abcov, l2m_acov, l2m_bcov)
%
% Simulated data for multi-level mediation:
% Generate data for a set of n participants with t observations each
% - Fixed population Path a and b effects (apop, bpop)
% - 2nd-level moderator with population Var = 1, 
%   which may be correlated with individual differences in individual Path a and b
% - Gaussian random noise, no correlated errors across participants
% - Complete mediation (no direct effect)
%
% n         Num of participants (2nd level units)
% t         Num of 1st-level observations per participant
% s_g       Population variance (between-person, 2nd level) of paths a and b
% s_a       variance of within-person measurement error for Path a, sqrt(s_a) is std. deviation
% s_b       variance of within-person measurement error for Path b, sqrt(s_b) is std. deviation
% apop      Path a population fixed effect
% bpop      Path b population fixed effect
% abcov     Population covariance between a and b paths
% l2m_acov  Population covariance for Path a and level-2 moderator
% l2m_bcov  Population covariance for Path b and level-2 moderator
%
% Optional inputs:
% -------------------------------------------------------------------------
% 'niter', followed by number of iterations
% 'randseed', followed by random-number seed.  Use this if you want to
% replicate an exact set of random samples across different simulations.
% The default behavior is to re-set the seed each iteration.
%
% Simple case, no a*b cov, use all default params
% [pvals, stats] = mediation_sim3();
%
% Simple case, no a*b cov
% mediation_sim3(20, 40, .5, 2, 2, .3, .3, 0, 0, 0, 'k')
%
% Add a*b cov
% mediation_sim3(20, 40, .5, 2, 2, .3, .3, .5, 0, 0, 'k')
%
% Add level 2 moderator that is correlated
% mediation_sim3(20, 40, .5, 2, 2, .3, .3, .5, .5, .5, 'k')
%
% blue:  L2 cov is unrelated to mediation
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3, 0, 0, 'b');
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3,  0, 0, 'b');
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3,  0, 0, 'b');
%
% % red: with L2 cov that is related to a and b
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3, .6, .6, 'r');
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3, .6, .6, 'r');
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3, .6, .6, 'r');
%
% % green:  no a*b cov
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, 0, 0, 0, 'g');
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, 0,  0, 0, 'g');
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, 0,  0, 0, 'g');
%
% lt green:  no a*b cov, but stronger a and b
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .6, .6, 0, 0, 0, [.5 1 .5]);
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .6, .6, 0,  0, 0, [.5 1 .5]);
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .6, .6, 0,  0, 0, [.5 1 .5]);
%
% % orange:  no a*b cov, stronger a and b, L2 cov that is related
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .6, .6, 0, .5, .5, [1 .5 0]);
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .6, .6, 0,  .5, .5, [1 .5 0]);
% [pvals, results_table] = mediation_sim3(20, 40, .5, 2, 2, .6, .6, 0,  .5, .5, [1 .5 0]);
%
% Example: Replicate an analysis exactly
% [pvals, results_table, randseed] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3, .6, .6, 'r', 'niter', 100);
% [pvals2, results_table2] = mediation_sim3(20, 40, .5, 2, 2, .3, .3, .3, .6, .6, 'r', 'randseed', randseed, 'niter', 100);
% figure; plot_correlation_samefig(pvals(1:100, 5), pvals2(1:100, 5));

% DEFAULT SETTINGS
% -------------------------------------------------------------------------
omitl2m = 0;
niter = 200;

if nargin < 10
    % default case
    
    n = 20;
    t = 40;
    
    s_g = .5;  % sigma - between (2nd level)
    s_a = 2; % sigma - within on a
    s_b = 2;
    
    apop = .3;
    bpop = .3;
    abcov = .5;
    
    l2m_acov = 0;
    l2m_bcov = 0;
    
    color = [1 0 0];
    doplot = true;
    
end


% OPTIONAL INPUTS
% -------------------------------------------------------------------------
randseed = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'randseed'
                randseed = varargin{i+1};
                rng(randseed);              % use previous value
                
            case {'niter', 'iterations'}, niter = varargin{i+1};
                    
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(randseed)
    
    rng('default')  % restore rng in case we're used legacy before
    randseed = rng('shuffle'); % re-set random number generator
    rng(randseed)       % so that subsequent runs will be the same
end

% note: for speed increase, should eventually replace this with pre-defined
% pval and abcovs vars with correct size
clear pvals abcovs

fprintf('Iteration %4.0f', 0);

for i = 1:niter
    fprintf('\b\b\b\b%4.0f', i);
    [x, m, y, l2m] = mediation_gen_data1(n, t, s_g, s_a, s_b, apop, bpop, abcov, l2m_acov, l2m_bcov);
    
    if omitl2m
        [paths, stats] = mediation(x, y, m);
    else
        [paths, stats] = mediation(x, y, m, 'l2m', l2m);
    end
    
    pvals(i, :) = stats.p(1, :);
    tmp = corrcoef(stats.paths(:, 1), stats.paths(:, 2));
    abcovs(i, 1) = tmp(1,2);
end

fprintf('\n');

pvals = [pvals abcovs];
stats.names{end+1} = 'abcov';

if doplot
    
    ncols = size(pvals, 2);
    create_figure('Msim3 Results2', 1, ncols);
    
    for i = 1:ncols
        subplot(1, ncols, i);
        boxplot(pvals(:, i));
        title(stats.names{i});
    end
    
    create_figure('ab_hist_pvals', 1, 1, 1); hold on;
    hh = hist(pvals(:, 5), 0:.05:1);
    plot(0:.05:1, hh, '-', 'Color', color, 'LineWidth', 2); %rand(1,3));
    title('Histogram of P-values for a*b');
    
end

% Table
sig_rate_05 = sum(pvals < 0.05) ./ niter;

results_table = array2table(sig_rate_05, 'VariableNames', strrep(stats.names, '''', 'prime'));

end % main function


function [x, m, y, l2m] = mediation_gen_data1(n, t, s_g, s_a, s_b, apop, bpop, abcov, l2m_acov, l2m_bcov)
%
% [x, m, y, l2m] = mediation_gen_data1(n, t, s_g, s_a, s_b, apop, bpop, abcov, l2m_acov, l2m_bcov)
%
% Simulated data for multi-level mediation:
% Generate data for a set of n participants with t observations each
% - Fixed population Path a and b effects (apop, bpop)
% - 2nd-level moderator with population Var = 1, 
%   which may be correlated with individual differences in individual Path a and b
% - Gaussian random noise, no correlated errors across participants
% - Complete mediation (no direct effect)
%
% n         Num of participants (2nd level units)
% t         Num of 1st-level observations per participant
% s_g       Population variance (between-person, 2nd level) of paths a and b
% s_a       variance of within-person measurement error for Path a, sqrt(s_a) is std. deviation
% s_b       variance of within-person measurement error for Path b, sqrt(s_b) is std. deviation
% apop      Path a population fixed effect
% bpop      Path b population fixed effect
% abcov     Population covariance between a and b paths
% l2m_acov  Population covariance for Path a and level-2 moderator
% l2m_bcov  Population covariance for Path b and level-2 moderator

% tor: i *think* this may no longer be necessary, as matlab has updated its
% rng methods
% if isempty(randseed) % default
%     randn('state',sum(100*clock))
% else
%     % do not re-set
% end

% Population covariance matrix 
% Diagonals are between-person (rfx) variances of path a and b.
covmtx = [s_g abcov l2m_acov; abcov s_g l2m_bcov; l2m_acov l2m_bcov 1]; % path a, path b, level 2 mod

% Sample n "true" individual-participant path coefficients from this covariance matrix
% Means are apop and bpop for path a and b, and 0 for a level-2 moderator
tmp = mvnrnd([apop bpop 0], covmtx, n);
a = tmp(:, 1);
b = tmp(:, 2);
l2m = tmp(:, 3);

% s_a is variance of within-person measurement error, sqrt(s_a) is std. deviation.

for i = 1:n
    
    x{i} = randn(t, 1);
    m{i} = a(i) .* x{i} + sqrt(s_a) .* randn(t, 1);
    y{i} = b(i) .* m{i} + sqrt(s_b) .* randn(t, 1);
    
end

end

%%




