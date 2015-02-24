function [paths, varargout] = mediation_threepaths(X, Y, M1, M2, varargin)

% function [paths, varargout] = mediation_threepaths(X, Y, M1, M2, varargin)
%
% Usage: This function tests the three-path mediation effect (X -> M1 -> M2 -> Y).
%
% This is based on Tor Wager's original mediation function (mediation.m) and also 
% Taylor et al.(2007)'s three path mediation analysis. Now, this is working only 
% with the multilevel and bootstrap options. X, Y, M1, and M2 have to have
% each subject's data in a cell. default bootsamples = 1000.
% 
% The 'logistic' (or 'logit') option was added by Wani (01/13/2013). 
%       You can use this option when Y variable is categorical. 
%       Y can be binary (E.g., 1:true, 0:false) or multinomial 
%       (E.g., 1:Red,2:Blue,3:Green,4:Yellow,5:Black).
%       ref) Iacobucci (2012) J of Consumer Psych. 22, 592-594
%       
% The example is as follow: 
%
% Example)
% for i = 1:30
%   X{i} = rand(50,1); M1{i} = rand(50,1); M2{i} = rand(50, 1); 
%   Y{i} = rand(50,1); cov{i} = rand(50,2);
% end
% 
% [paths, stats] = mediation_threepaths(X, Y, M1, M2, 'covs', cov, ...
%       'verbose', 'plots', 'bootsamples', 10000);
% 
% Wani Woo, 11/04/2012
%
% Reference: 
% Taylor, MacKinnon, & Tein (2007). Tests of the Three-Path Mediated Effect. 
%     Organizational Research Methods, 11(2), 241?269. 


% -------------------------------------------------------------------------
% Setup inputs, print info to screen if verbose
% -------------------------------------------------------------------------
global mediation_covariates  % used in mediationfun and mediation_path_coefficients; empty = no covs

varargout = cell(1, nargout - 1);

% default
targetu = .20;              % proportion contribution of boot procedure to p-value
whpvals_for_boot = 6;   % indices of p-values, the min of which is used to determine boot samples needed
domultilev = 1;             % multi-level analysis (if N replications > 1)

if ~iscell(Y), error('The three-paths mediation is working only for multi-level analysis.');
else N = length(Y);
end

vnames = {'X' 'Y' 'M1' 'M2'};

% varargin setting
boottop = 1;                  % bootstrap instead of OLS
verbose = 0;                % verbose output
bootsamples = 1000;         % initial bootstrap samples
mediation_covariates = [];
doplots = 0;
X_2ndlevel = ones(N, 1);
logistic_Y = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'verbose', verbose = 1;
            case 'bootsamples', bootsamples = varargin{i+1};
            case {'covs', 'covariates', 'mediation_covariates'}, mediation_covariates = varargin{i+1};
            case 'names', vnames = varargin{i+1};
            case {'plot', 'plots'}, doplots = 1;
            case {'logit', 'logistic', 'logistic_Y'}, logistic_Y = 1;
        end
    end
end

if verbose
    nobs_tmp = size(Y, 1);
    fprintf('Mediation analysis\n\nObservations: %3.0f, Replications: %3.0f\n', nobs_tmp, N);
    fprintf('Predictor (X): %s, Outcome (Y): %s: Mediator (M): %s\n', vnames{1}, vnames{2}, vnames{3});
    
    nms = {'No' 'Yes'};
    
    fprintf('\nCovariates: %s', nms{~isempty(mediation_covariates) + 1})
    if ~isempty(mediation_covariates)
        fprintf(', %3.0f columns', size(mediation_covariates, 2));
    end
    fprintf('\n')
    fprintf('\nMulti-level analysis.\n')
    fprintf('Options:\n\tMultilevel weights: %s\n\tBootstrap: %s\n\tLogistic(Y): %s\n', nms{domultilev+1}, nms{boottop+1}, nms{logistic_Y+1});
end

if verbose, totalt = clock; end

npaths = 6; % b1, b2, b3, c, c', and b1b2b3

% Initialize outputs

paths = NaN .* zeros(N, npaths);
sterrs = Inf .* ones(N, npaths);
[b1betas, b2betas, b3betas, cpbetas, cbetas] = deal(cell(1, N));
[Vxm1, Vm1m2, Vm2y, Vxy, residm1, residm2, residy, residy2] = deal(cell(N, 1));

n = zeros(N, 1);
nrm = zeros(N, 1);

mediationfun = @fast_ols_ab;

% Weighted mean function:
% Used in getstats to get mean bootstrap coefficients
% Second level function for getting weights across subjects
% Faster than looping, and faster than using mean() for equal weights
% Also makes code simpler below.
wmean = @(paths, w) diag(w'*paths)';
w = ones(N) ./ N;  % start with weights equal

% =========================================================================
%
% * Single-level model, run for each replication N
%
% =========================================================================

MCOVS = mediation_covariates;

for i = 1:N

    x = X{i}; 
    y = Y{i}; 
    m1 = M1{i}; 
    m2 = M2{i}; 
    
    if ~isempty(MCOVS)
        mediation_covariates = MCOVS{i}; 
    end
    
    % --------------------------------------------
    % Remove NaNs and Check for bad data
    % --------------------------------------------
    % Remove NaNs from all vars, casewise
    [nanvec, x, y, m1, m2, mediation_covariates] = nanremove(x, y, m1, m2, mediation_covariates);
    
    % ***NOTE: Check to make sure no exact duplicates in x , y , m
    % ***check to make sure length valid obs > 1
    
    tol = .000001;
    isbad = isempty(x) || any(max(abs(m1)) < tol) || all(abs(m1(:,1) - m1(1)) < tol);
    isbad = isbad || any(max(abs(x)) < tol) || all(abs(x(:,1) - x(1)) < tol);
    isbad = isbad || any(max(abs(y)) < tol) || all(abs(y(:,1) - y(1)) < tol);
    isbad = isbad || any(max(abs(m2)) < tol) || all(abs(m2(:,1) - m2(1)) < tol);
    
    % Note: Added 2011 for additional error checking by Tor
    if ~isempty(mediation_covariates)
        vv = var(mediation_covariates);
        isbad = isbad || any(max(abs(mediation_covariates(:))) < tol) || any(vv < tol);
    end
    
    if isbad
        % No valid data; skip this subject
        if verbose
            fprintf('\n')
            warning('mediation:BadData',['No valid data for subject ' num2str(i) ': Skipping.']);
        end
        
        stats = get_ols_stats(paths(i,:) , sterrs(i,:), n(i));  % setup dummy stats structure
        stats = add_info_to_stats(stats);
        stats.analysisname = 'Empty mediation';
        
        if N > 1
            [stats.beta, stats.p, stats.Z] = deal( NaN .* zeros(size(X_2ndlevel, 2), size(stats.mean, 2)) );
            
            if i == 1, wistats = stats; wistats.mediation_covariates = mediation_covariates;  end
            wistats = collect_within_stats(wistats, stats, x, y, m1, m2, i);
        else
            [stats.beta, stats.p, stats.Z] = deal( NaN .* zeros(1, size(stats.mean, 2)) );
            
        end
        
        continue
    end
    
    
    % --------------------------------------------
    % Compute path coefficients (and standard errors in some cases)
    % If bootstrapping is 'on', skips standard errors
    % If bootstrapping is 'off', returns std. errors, intercept, and n
    % needed for OLS stats
    % --------------------------------------------
    
    
    [paths(i,:), b1betas{i}, b2betas{i}, b3betas{i}, cpbetas{i}, cbetas{i}, sterrs(i,:), intcpt, n(i), residm1{i}, residm2{i}, residy{i}, residy2{i},...
        Vxm1(i),Vm1m2(i),Vm2y(i),Vxy(i)] = mediation_path_coefficients_threepaths(x, y, m1, m2, domultilev, logistic_Y);
    
    % Multilevel only: save std. deviations for getting standardized paths
    if N > 1, nrm(i, 1) = std(x) ./ std(y); end
    
    stats = get_ols_stats(paths(i, :), sterrs(i, :), n(i));
    
    if exist('stats', 'var')
        stats.residm1 = residm1;
        stats.residm2 = residm2;
        stats.residy = residy;
        stats.residy2 = residy2;
    end
    
    
    % --------------------------------------------
    % multi-subject (multilevel)
    % collect within-subject stats; assumes independence of trials
    % --------------------------------------------
    if N > 1
        if i == 1, wistats = stats; end
        wistats = collect_within_stats(wistats, stats, x, y, m1, m2, i);
        wistats.mediation_covariates = mediation_covariates;
    end
    
end


% =========================================================================
%
% Done with single-level estimation.  Now finish multi-level stuff.
%
% =========================================================================


% --------------------------------------------
% summarize within-subjects bootstrapped stats
% if multi-level
% --------------------------------------------
if N > 1
    collect_within_stats_summary;
    varargout{2} = wistats;
end

if verbose, fprintf(' Done.\n'); end


% --------------------------------------------------
%
% * Check if bad data and exit gracefully if so
%
% --------------------------------------------------

if N > 1 && all(isnan(paths(:)))  %~any(wistats.df < 1)
    % NO valid data; exit gracefully without continuing computation
    
    if verbose, disp('No valid data for any subjects!'); end
    
    stats2 = getstats(paths);
    stats2 = add_info_to_stats(stats2);
    
    [stats2.beta, stats2.p, stats2.Z] = deal( NaN .* zeros(size(X_2ndlevel, 2), size(stats2.mean, 2)) );
    
    varargout{1} = stats2;
    varargout{2} = wistats;
    
    return
    
elseif all(isnan(paths(:)))
    % Single-level without valid data; exit gracefully
    
    if verbose, disp('No valid data for any subjects!'); end
    varargout{1} = stats;
end


% --------------------------------------------------
%
% * Second-level model, if N > 1
%
% --------------------------------------------------
% If N > 1, then get group stats (2nd level)
% --------------------------------------------------
if N > 1
    if boottop
        means = [];
        second_level_bootstrap(); % only works for bootstrap
        
%     elseif dosignperm
%         means = [];
%         second_level_permtest();
%         
%     else
%         % Check for and remove bad/missing 2nd level units
%         % -------------------------------------------------
%         whomit = any(w == 0, 2) | any(paths == 0, 2) | any(isnan(paths), 2) | any(isnan(X_2ndlevel), 2);
%         if verbose && any(whomit), warning('mediation:BadData', 'Some 2nd-level units have missing or bad values.'); end
%         isOK = ~whomit;
%         
%         if verbose && all(whomit), warning('mediation:BadData', 'No valid data in 2nd-level analysis.'); end
%         
%         % If all bad, run it with all data, and we'll get NaNs and the
%         % right kind of structure back so we don't crash
%         if all(whomit), isOK = whomit; end
%         % -------------------------------------------------
%         
%         % OLS case (2nd-level regression)
%         [b, s2between_ols, stats2] = scn_stats_helper_functions('gls', paths(isOK, :), w(isOK, :), X_2ndlevel(isOK, :));
%         
%         stats2.std = stats2.ste .* sqrt(N);
%         
%         % weight, and do it again with weights
%         get_weights_based_on_varcomponents()
%         
%         [b, s2between_ols, stats2] = scn_stats_helper_functions('gls', paths(isOK, :), w(isOK, :), X_2ndlevel(isOK, :));
%         
%         stats2.name = '2nd level statistics';
%         
%         % THis is still wrong
%         stats2.Z = abs(norminv(stats2.p ./ 2)) .* sign(stats2.beta);  %repmat(sign(stats2.beta), 1 + num_additionalM, 1);
    end
    
    stats2 = add_info_to_stats(stats2);
    stats2.analysisname = 'Multi-level model';
    
    varargout{1} = stats2;
    
    if doplots && (boottop || dosignperm)
        plot_hists(means, vnames);
    end
end


% --------------------------------------------------
%
% * Output
%
% --------------------------------------------------

% path and slope plots
% 
% if doplots && N > 1
%     % Multi-level ----------------------------------
%     
%     if dorobust
%         disp('No robust slope plots yet.')
%     else
%         % plots of slopes for each subject
%         mediation_plots(stats2, 'slopes', varargin);
%         mediation_scatterplots(stats2);
%         mediation_plots(stats2, 'abcov');
%         
%     end
%     
%     try mediation_path_diagram(stats2);
%     catch, disp('Error in mediation_path_diagram');
%     end
%     
%     if dolatent
%         try plot_hrf_in_latent_model(wistats);
%         catch, disp('Error in plot_hrf_in_latent_model');
%         end
%     end
%     
%     % Plot individual effects
%     create_figure('Individual Effects'); barplot_columns(paths, 'Individual Effects', [], 'nofig', 'plotout','number');
%     set(gca, 'XTickLabel', stats.names);
%     ylabel('Beta (slope) values');
%     
% elseif doplots && N == 1
%     % Single-level ----------------------------------
%     
%     try mediation_path_diagram(stats);
%     catch, disp('Error in mediation_path_diagram');
%     end
%     
%     mediation_scatterplots(stats);
% end

if verbose && exist('stats2', 'var')
    if exist('wistats', 'var') && ~isempty(wistats)
        print_outcome(stats2, wistats)
    else
        print_outcome(stats2)
    end
    
% elseif verbose && N == 1 && exist('stats', 'var')
%     print_outcome(stats)
end

if verbose, fprintf('Total time: %3.0f s\n', etime(clock, totalt));  end


% error checking:
% estab = stats2.mean(1) * stats2.mean(2) + 2 * cov(means(:,1:2)) % approx; cov should be weighted?

% if dosave
%     if doplots
%         save_figures(verbose);
%     end
%     
%     if verbose
%         diary off
%         fprintf('All finished: Closed file Mediation_Output.txt\n')
%     end
% end



% _________________________________________________________________________
%
%
%
% * Inline functions
%
%
%
%__________________________________________________________________________


% -------------------------------------------------------------------------
% Permutation test for 2nd level
% -------------------------------------------------------------------------
%     function second_level_permtest()
%         persistent permsign
%         
%         if ~persistent_perms
%             permsign = [];
%         end
%         
%         if verbose, fprintf('Nonparametric test with %3.0f permutations...', bootsamples); t12 = clock; end
%         
%         % start with weights all equal whether multilevel or not
%         
%         [p, Z, xbar, permsign, means] = permute_signtest(paths, bootsamples, w, permsign);
%         
%         stats2 = [];
%         stats2.name = '2nd level statistics';
%         %%%%stats2.names = {'a' 'b' 'c''' 'c' 'ab'};
%         stats2.mean = xbar;
%         stats2.ste = std(means);
%         stats2.std = std(paths);
%         stats2.Z = Z;
%         stats2.p = p;
%         
%         
%         % get weights, if multilevel
%         % returns w and stats2.w, etc.
%         if domultilev
%             if verbose, fprintf('Re-permuting with multi-level weights.\n'), end
%             get_weights_based_on_varcomponents();
%             % do full test based on new weights
%             [p, Z, xbar, permsign, means] = permute_signtest(paths, bootsamples, w, permsign);
%             stats2.mean = xbar;
%             
%             % *** could replace with weighted std; also, check this
%             stats2.ste = std(means);
%             stats2.std = std(paths);
%             stats2.Z = Z;
%             stats2.p = p;
%         end
%         
%         
%         if verbose, fprintf(' Done in %3.0f s \n', etime(clock, t12)); end
%     end


% -------------------------------------------------------------------------
% Bootstrap test for 2nd level
% -------------------------------------------------------------------------
    function second_level_bootstrap()
        whnan = nanremove(paths, w);
        whgood = ~whnan;
        
        if all(whnan), error('No valid subjects for 2nd level bootstrap...check data.'); end
        
        if verbose, fprintf('Bootstrapping %3.0f samples...', bootsamples); end
        
        % set up boot samples; make sure all are valid (otherwise
        % warnings/problems with categorical predictors + small samples)
        % and initialize random number generator
        %  subjindx = (1:size(X_2ndlevel, 1))';  % indices of subjects (2nd level units), to get bootstrap samples, etc.
        subjindx = (1:size(X_2ndlevel(whgood), 1))';  % indices of subjects (2nd level units), to get bootstrap samples, etc.
        bootsam = setup_boot_samples(subjindx, bootsamples);
        
        
        if verbose && any(whnan), fprintf('Warning! %3.0f 2nd-level observations are invalid/have missing data.\n', sum(whnan)); end
        
        % start with weights all equal whether multilevel or not
        if sum(whgood) < 2
            means = zeros(bootsamples, npaths);
        else
            means = bootstrp_havesamples(bootsam, wmean, paths(whgood,:), w(whgood,:));
        end
        stats2 = getstats(means);
        
        stats2.name = '2nd level statistics';
        
        % std is needed for weighted mean option below.
        stats2.std = stats2.ste .* sqrt(N);
        
        % Get boot samples needed
        % Add samples as needed to ensure p-value is reasonably accurate
        % Be: Boot samples needed so that Expected (average) p-value
        % error due to limited-sample bootstrap is targetu*100 %
        
        [Be, alphaaccept] = get_boot_samples_needed(stats2.p, whpvals_for_boot, targetu, bootsamples, verbose); % uses whpvals_for_boot, returns Be
        
        if Be > 0 % if/else for naming convention, to keep code clearer
            bootsam_plusBe = [bootsam setup_boot_samples(subjindx, Be)];
        else
            bootsam_plusBe = bootsam;
        end
        
        if verbose && boottop && (Be > 0 || domultilev), fprintf('Bootstrapping...'); end
        if verbose, t12 = clock; end
        
        % get weights, if multilevel
        % returns w and stats2.w, etc.
        if domultilev
            get_weights_based_on_varcomponents();
            
            % get full bootstrap sample based on new weights, whether or
            % not we've added samples
            if sum(whgood) > 2
                means = bootstrp_havesamples(bootsam_plusBe, wmean, paths(whgood,:), w(whgood,:));
            end
        else
            % no weighting; summary stats.
            % add additional boot samples only if necessary, and only
            % additional samples -- to save time
            if Be > 0 && sum(whgood) > 2
                means = [means; bootstrp_havesamples(bootsam_plusBe(:, bootsamples+1:end), wmean, paths(whgood,:), w(whgood,:))];
            end
        end
        
        stats2 = getstats(means, stats2); % last input preserves existing info
        
        % use original weighted mean, not mean of bootstrap samples
        stats2.mean = wmean(paths(whgood,:), w(whgood,:));
        
        % std is needed for weighted mean option below.
        stats2.std = stats2.ste .* sqrt( sum(whgood) );
        
        % bias correction for final bootstrap samples
        % [p, z] = bootbca_pval(testvalue, bootfun, bstat, stat, [x], [other inputs to bootfun])
        stats2.prctilep = stats2.p; % percentile method, biased
        [stats2.p, stats2.z] = bootbca_pval(0, wmean, means, wmean(paths(whgood,:), w(whgood,:)), paths(whgood,:), w(whgood,:));
        
        %[dummy, dummy, stats2.z, stats2.p] = bootbca_pval_onetail(0, wmean, means, wmean(paths(whgood,:), w(whgood,:)), paths(whgood,:), w(whgood,:));
        
        stats2.biascorrect = 'BCa bias corrected';
        stats2.alphaaccept = alphaaccept;
        
        if verbose, fprintf(' Done in %3.0f s \n', etime(clock, t12)); end
    end



% -------------------------------------------------------------------------
% Add bootstrap samples to run -- FIRST LEVEL
% -------------------------------------------------------------------------
    function add_boot_samples_needed()
        
        % Be: Boot samples needed so that Expected (average) p-value error due to limited-sample bootstrap is targetu*100 %
        
        [Be, alphaaccept] = get_boot_samples_needed(stats.p, whpvals_for_boot, targetu, bootsamples, verbose);
        
        bootsam_Be = setup_boot_samples(x, Be);
        
        if Be > 0
            if verbose, fprintf(' Adding %3.0f samples ', Be), t12 = clock; end
            if dorobust
                if isempty(mediation_covariates)
                    bootpaths = [bootpaths; bootstrp_havesamples(bootsam_Be, @fast_robust_ab, x, y, m)];   % twice as fast as long version
                else
                    bootpaths = [bootpaths; bootstrp_havesamples(bootsam_Be, @fast_robust_ab, x, y, m, mediation_covariates)];
                end
                
            else
                if isempty(mediation_covariates)
                    bootpaths = [bootpaths; bootstrp_havesamples(bootsam_Be, @fast_ols_ab, x, y, m, intcpt)];
                else
                    bootpaths = [bootpaths; bootstrp_havesamples(bootsam_Be, @fast_ols_ab, x, y, m, intcpt, mediation_covariates)];
                end
            end
            
            stats = getstats(bootpaths);
            stats.vnames = vnames;
            stats.bootsamples = size(bootpaths, 1);
            stats.alphaaccept = alphaaccept;
            
            if verbose, fprintf(' Done adding in %3.0f s \n', etime(clock, t12)); end
        end
    end

% -------------------------------------------------------------------------
% For multilevel model: Get weights
% 1st-level standard errors & variance components
% -------------------------------------------------------------------------
    function get_weights_based_on_varcomponents()
        if verbose && boottop, fprintf('Estimating variance components based on %3.0f bootstrap samples\n', bootsamples); end
        
        % estimates of variance components -- naive est.
        stats2.wi_var_est = (wistats.avg_wi_std .^2); % Use STD, not ste of param est. for individual
        stats2.total_var_est = (stats2.std.^2);
        
        % get weights : work in progress. Sep. for each equation.
        
        wh_omit = false(1, length(Vxm1));
        
        for ii = 1:size(Vxm1, 1)
            if isempty(Vxm1{ii}) || isempty(Vm1m2{ii}) || isempty(Vm2y{ii}) || isempty(Vxy{ii})
                wh_omit(ii) = 1;
            end
        end

        Vxm1(wh_omit) = []; % can have multiple cols, one per mediator (one per x->m eq.)
        Vm1m2(wh_omit) = []; % can have multiple cols, one per mediator (one per x->m eq.)
        Vm2y(wh_omit) = [];
        Vxy(wh_omit) = [];
        
        %if multiple mediators, we have to reorganize a bit,
        %bec. of how things are stored in mediation_path_coefficients.
        tmp_b1betas = cat(2, b1betas{:});
        tmp_b2betas = cat(2, b2betas{:});
        tmp_b3betas = cat(2, b3betas{:});
        
        my_b1betas{1} = tmp_b1betas;
        my_b2betas{1} = tmp_b2betas;
        my_b3betas{1} = tmp_b3betas;
            
        % Note: do not use EBayes params for inference; they are individual
        % slopes, but they have been shrunk towards the group mean...not
        % confident in validity of inference.  use betas.
        % A effect 1
        % ------------------------
        mycol = 1;
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(cat(2, b1betas{:}), Vxm1(:, 1));
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_std(mycol) = sqrt(sigma2_b_est(1));
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        stats2.ebayes_bstar(:, mycol) = b_star(1, :)';
        
        mycol = 2;  % output columns storing b effects for all mediators
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(cat(2, b2betas{:}), Vm1m2(:, 1));
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_std(mycol) = sqrt(sigma2_b_est(1));
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        stats2.ebayes_bstar(:, mycol) = b_star(1, :)';

        mycol = [3 4];  % output column for cp effect
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(cat(2, b3betas{:}), Vm2y(:, 1));
        stats2.wls_mean(mycol) = gam_hat(1:2);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1:2);
        stats2.btwn_std(mycol) = sqrt(sigma2_b_est(1:2));
        stats2.btwn_var_est(mycol) = sigma2_b_est(1:2);
        stats2.ebayes_bstar(:, mycol) = b_star(1:2, :)';
        
        mycol = 5;
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(cat(2, cbetas{:}), Vxy);
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_std(mycol) = sqrt(sigma2_b_est(1));
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        stats2.ebayes_bstar(:, mycol) = b_star(1, :)';
        
        % A*B effects
        % ------------------------
        mycol = 6;
        % need to figure out ab var...
        % ****NOTE: this is a fix for the case in which some subjects have missing data.
        wh_include = find(~isnan(paths(:, mycol)));
        for i = 1:length(wh_include), Vab{i} = sterrs(wh_include(i), mycol) .^ 2; end   % get V within subj for each ab path
            
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(paths(wh_include, 6)', Vab);
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_std(mycol) = sqrt(sigma2_b_est(1));
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        stats2.ebayes_bstar(:, mycol) = b_star(1, :)';
        
        
        % make sure we deal with missing data subjects
        if length(wh_include) < N
            w = zeros(size(paths));
            w(wh_include, :) = stats2.w;
            stats2.w = w;
        end
        
        w = stats2.w;
        
        whnan = nanremove(w);
        if all(whnan), error('No valid subjects after reweighting...invalid estimates for at least one subject...check data.'); end
        
    end % get_weights subfunction




% -------------------------------------------------------------------------
% Collect within-subject stats to save as output
% -------------------------------------------------------------------------
    function collect_within_stats_summary()
        wistats.name = '1st level statistics: summary';
        wistats.vnames = vnames;
        wistats.avg_wi_std = nanmean(wistats.std);
        wistats.avg_wi_ste = nanmean(wistats.ste);
        wistats.avg_wi_d = mean(wistats.mean ./ max(wistats.std, eps)); % effect size (save)
        
        goodobs = ~isinf(wistats.sx);
        if ~any(goodobs)
            wistats.avg_sx = NaN;
            wistats.avg_sy = NaN;
            wistats.avg_r = NaN;
        else
            wistats.avg_sx = mean(wistats.sx( goodobs ));
            wistats.avg_sy = mean(wistats.sy( goodobs ));
            wistats.avg_r = mean(wistats.r( goodobs ));
        end
        
        if exist('residy', 'var')
            wistats.residy = residy;
            wistats.residy2 = residy2;
        end
        
        % make sure
    end


% -------------------------------------------------------------------------
% Update top-level stats structure and save final path estimates and
% input options
% -------------------------------------------------------------------------
    function stats_struct = add_info_to_stats(stats_struct)
        ynstr = {'No' 'Yes'};
        
        inputOptions = struct('N', N, 'npaths', npaths, 'initial_bootsamples', bootsamples, 'targetu', targetu);
        
        inputOptions.vnames = vnames; % do we need this?
        %%%%inputOptions.names = stats_struct.names;
        
        % save data
        inputOptions.X = X;
        inputOptions.Y = Y;
        inputOptions.M1 = M1;
        inputOptions.M2 = M2;
        inputOptions.mediation_covariates = MCOVS; %mediation_covariates;
        
        % save options
        % inputOptions.bootfirstlevel = ynstr{dobootfirstlevel + 1};
        inputOptions.multilevel = ynstr{domultilev + 1};
        inputOptions.bootstrap = ynstr{boottop + 1};
        
        if boottop || dosignperm
            if exist('means', 'var') % for multi-level
                inputOptions.final_bootsamples = size(means, 1);
                
            elseif exist('bootpaths', 'var') % kludge for single-level model
                inputOptions.final_bootsamples = size(bootpaths, 1);
            end
        end
        
        stats_struct.paths = paths;
        
        stats_struct.w = w;
        
        % standardized paths
        stats_struct.stdpaths = paths .* repmat(nrm, 1, npaths);
        stats_struct.b1betas = b1betas;
        stats_struct.b2betas = b2betas;
        stats_struct.b3betas = b3betas;
        stats_struct.cpbetas = cpbetas;
        stats_struct.cbetas = cbetas;
        
        stats_struct.inputOptions = inputOptions;
        
        stats_struct = get_stats_names(stats_struct);
        
    end
end





% _________________________________________________________________________
%
%
%
% * Sub-functions
%
%
%
%__________________________________________________________________________






% -------------------------------------------------------------------------
% Figure out how many more bootstrap samples to run (general)
% -------------------------------------------------------------------------
function [Be, alphaaccept] = get_boot_samples_needed(p, whpvals_for_boot, targetu, bootsamples, verbose)
minp = min(p(whpvals_for_boot));
[B95, Be, alphaaccept] = Bneeded(max(.005, minp), targetu);
Be = max(0, ceil(Be - bootsamples));   % additional number to run
if verbose, fprintf(' Min p-value is %3.6f. Adding %3.0f samples\n ', minp, Be), end
end


% -------------------------------------------------------------------------
% Get statistic structure from bootstrapped samples, including p-values and conf. intervals
% -------------------------------------------------------------------------
function stats = getstats(bp, varargin)

if ~isempty(varargin), stats = varargin{1}; end

stats.names = {'b1' 'b2' 'b3' 'c''' 'c' 'b1b2b3'};

stats.mean = mean(bp);
stats.ste = std(bp);

stats.p = 2.*(min(sum(bp<=0), sum(bp>=0)) ./ size(bp, 1));

% avoid exactly-zero p-values
% replace with 1/# bootstrap samples
stats.p = max(stats.p, 1./size(bp, 1));

end

% -------------------------------------------------------------------------
% Get statistic structure from OLS regression, including p-values and conf. intervals
% -------------------------------------------------------------------------
function stats = get_ols_stats(bp, sterrs, n)
stats.names = {'b1' 'b2' 'b3' 'c''' 'c' 'b1b2b3'};

stats.mean = bp;
stats.ste = sterrs;
stats.t = bp ./ sterrs;
stats.df = [n-1 n-2 n-3 n-3 n-1 n-3];   % not sure about last one ***check***
stats.p = min(1, (2 .* (1 - tcdf(abs(stats.t), stats.df))));

% don't do to save time
% stats.ci95 = [prctile(bp, 97.5); prctile(bp, 2.5)];
% stats.ci90 = [prctile(bp, 95); prctile(bp, 5)];
% stats.ci99 = [prctile(bp, 99.5); prctile(bp, .5)];
end


% -------------------------------------------------------------------------
% Collect structure of within-subject statistics
% -------------------------------------------------------------------------
function wistats = collect_within_stats(wistats, stats, x, y, m1, m2, i)
if i == 1
    % only for bootstrapped
    % %         if isfield(wistats, 'ci95'), wistats = rmfield(wistats, 'ci95'); end
    % %         if isfield(wistats, 'ci90'), wistats = rmfield(wistats, 'ci90'); end
    % %         if isfield(wistats, 'ci99'), wistats = rmfield(wistats, 'ci99'); end
    wistats = rmfield(wistats, 'p');
end

wistats.mean(i,:) = stats.mean;
wistats.ste(i,:) = stats.ste;
wistats.std(i,:) = stats.ste .* sqrt(size(x, 1));

if isempty(y) || all(y == y(1)), wistats.sy(i, 1) = Inf; else wistats.sy(i, 1) = std(y); end
if isempty(x) || all(x == x(1)), wistats.sx(i, 1) = Inf; else wistats.sx(i, 1) = std(x); end
if isempty(m1) || all(m1 == m1(1)), wistats.sm1(i, 1) = Inf; else wistats.sm1(i, 1) = std(m1); end
if isempty(m2) || all(m2 == m2(1)), wistats.sm2(i, 1) = Inf; else wistats.sm2(i, 1) = std(m2); end

tmp = corrcoef(x, y);
wistats.r(i, 1) = tmp(1, 2);

end


%__________________________________________________________________________
%
%
% *** Path Computation Support ***
% see also mediation_path_coefficients.m
%__________________________________________________________________________


% -------------------------------------------------------------------------
% Compute critical OLS paths in a faster way, for bootstrap
% -------------------------------------------------------------------------
function paths = fast_ols_ab(x, y, m, intcpt, varargin)

%%% **** action item: make subfunction; persistent px pmx
% if search y, save px, pmx, a, c
% if search x, save nothing
% if search m, save px, c
% Check whether grouping px * [m y] is faster
% check how regress, glmfit calculate betas, sterrs --

nmediators = size(m, 2);
mediation_covariates = [];  % need this in this format for bootstrapping
if ~isempty(varargin), mediation_covariates = varargin{1}; end

xx = [x intcpt mediation_covariates];

% notes on speed:
% X \ y is 10 x faster than pinv(X)*y, and inv(X'X)X'y is 3 x faster than
% pinv

px = inv(xx'*xx)*xx'; % X beta-forming matrix
mx = [m xx];
pmx = inv(mx'*mx)*mx'; % M+X beta-forming matrix

% Is X related to the mediator?
% Eq. 2, a = M ~ X
apaths = px * m;    % predictors x mediators
a = apaths(1,:);   % slopes for x -> all M

% Are M (b path) and X (c' path) independently related to Y?
% Eq. 3, [b, c'] = Y ~ X + M
bpaths = pmx * y;

b(1,:) = bpaths(1:nmediators)';       % slopes for each b effect
cp = bpaths(nmediators + 1);            % slope for cp effect

% is X related to Y without mediator?
% Eq. 1, Y ~ X
cpaths = px * y;
c = cpaths(1);


ab = a .* b;            % vector of ab for each mediator
paths = [a(1) b(1) cp c ab(1)];

for i = 2:nmediators
    paths = [paths a(i) b(i) ab(i)];
end

end


% -------------------------------------------------------------------------
% Compute critical robust IRLS paths in a faster way, for bootstrap
% -------------------------------------------------------------------------
function [paths] = fast_robust_ab(x, y, m, intcpt, varargin)
% varargin is for intercept, but not used; so we can use the same code
% above for OLS and robust with anon. f. handle mediationfun

%%% *** action item: re-do faster version of robustfit
%%% *** need to update to handle multiple mediators!!

mediation_covariates = [];
if ~isempty(varargin), mediation_covariates = varargin{1}; end

% Is X related to the mediator?
% Eq. 2, a = M ~ X
a = robustfit([x mediation_covariates], m);
a = a(2);

% Are M (b path) and X (c' path) independently related to Y?
% Eq. 3, [b, c'] = Y ~ X + M
b = robustfit([m x mediation_covariates], y);
cp = b(3);
b = b(2);

ab = a .* b;

c = robustfit([x mediation_covariates], y);
c = c(2);

paths = [a b cp c ab];
end



%__________________________________________________________________________
%
%
% *** Plotting and Display Support ***
%
%__________________________________________________________________________


% -------------------------------------------------------------------------
% Plot histograms of bootstrapped path coefficients
% -------------------------------------------------------------------------
function plot_hists(bootpaths, vnames)

tmp = bootpaths(:,4:5);  % c ones
tmp = tmp(:);
nbins = max(10, round(length(tmp) ./ 100));
nbins = min(nbins, 100);

[h, xx] = hist(tmp, nbins);

myperc = mean(xx) .* .2;
cxlimit = [min(min(xx)-myperc, 0 - myperc) max(xx)+myperc];
b1 = bootpaths(:,1);
b2 = bootpaths(:,2);
b3 = bootpaths(:,3);
cp = bootpaths(:,4);
c = bootpaths(:,5);
b1b2b3 = bootpaths(:,6);

fh = create_figure('Histogram_Plot', 1, 6);
widen_figure(fh);

% %     fh = findobj('Tag', 'Histogram_Plot');
% %     if isempty(fh)
% %         fh = tor_fig(1, 5);
% %         set(fh, 'Tag', 'Histogram_Plot');
% %         widen_figure(fh);
% %     else
% %         set(0, 'CurrentFigure', fh);
% %         hold on;
% %         arrayfun(@(axh)cla(axh), findobj(get(fh, 'Children'), 'Type', 'axes'))
% %     end

subplot(1, 6, 1);
shaded_hist(b1);
title(['b1: ' vnames{1} '->' vnames{3}]);
xlimit = [min(min(b1), 0 - .2*nanmean(b1)) max(b1)];
set(gca, 'XLim', xlimit);

subplot(1, 6, 2);
shaded_hist(b2);
title(['b2: ' vnames{3} '->' vnames{4}]);
xlimit = [min(min(b2), 0 - .2*nanmean(b2)) max(b2)];
set(gca, 'XLim', xlimit);

subplot(1, 6, 3);
shaded_hist(b3);
title(['b3: ' vnames{4} '->' vnames{2}]);
xlimit = [min(min(b3), 0 - .2*nanmean(b3)) max(b3)];
set(gca, 'XLim', xlimit);

subplot(1, 6, 4);
shaded_hist(cp, xx);
title(['c'':' vnames{1} '->' vnames{2}]);
set(gca, 'XLim', cxlimit);

subplot(1, 6, 5);
shaded_hist(c, xx);
title(['c: ' vnames{1} '->' vnames{2}]);
set(gca, 'XLim', cxlimit);

subplot(1, 6, 6);
shaded_hist(b1b2b3);
title(['b1b2b3: ' vnames{1} '->' vnames{2}]);
xlimit = [min(min(b1b2b3), 0 - .2*nanmean(b1b2b3)) max(b1b2b3)];
set(gca, 'XLim', xlimit)

drawnow
end



% -------------------------------------------------------------------------
% Create a shaded gray-scale histogram with .05 2-tailed in darker gray
% -------------------------------------------------------------------------
function shaded_hist(a, xx)
if nargin < 2
    nbins = max(10, round(length(a) ./ 20));
    nbins = min(nbins, 100);
    [h, xx] = hist(a, nbins);
else
    h = hist(a, xx);
end
han = bar(xx, h);
set(han, 'FaceColor', [.7 .7 .7]);
set(han, 'EdgeColor', [.7 .7 .7]);
wh = (xx > prctile(a, 97.5) | xx < prctile(a, 2.5));
h(~wh) = 0;
if any(wh')
    hold on;
    han = bar(xx, h);
    set(han, 'FaceColor', [.3 .3 .3], 'EdgeColor', [.3 .3 .3]);
end
plot_vertical_line(0);
end

% % % % %
% % % % % % -------------------------------------------------------------------------
% % % % % % Plot individual slopes of regressions, using conf. interval for x range
% % % % % % -------------------------------------------------------------------------
% % % % % function plot_slopes(abetas, bbetas, cbetas, cpbetas, X, M, additionalM, vnames, stats2)
% % % % %     % abetas and others are matrices: rows = [slope, intercept], cols =
% % % % %     % subjects
% % % % %
% % % % %     N = stats2.inputOptions.N;
% % % % %     npaths = stats2.inputOptions.npaths;
% % % % %
% % % % %     if iscell(abetas)
% % % % %         whgood = true(N, 1);
% % % % %         for i = 1:length(abetas)
% % % % %             if isempty(abetas{i}) || isempty(bbetas{i}) || isempty(cbetas{i}) || isempty(cpbetas{i})
% % % % %                 whgood(i) = 0;
% % % % %             end
% % % % %         end
% % % % %
% % % % %         abetas = cat(2, abetas{:});
% % % % %         bbetas = cat(2, bbetas{:});
% % % % %         cbetas = cat(2, cbetas{:});
% % % % %         cpbetas = cat(2, cpbetas{:});
% % % % %     else
% % % % %         whgood = ~any(isnan([abetas' bbetas' cbetas' cpbetas']), 2);
% % % % %     end
% % % % %
% % % % %
% % % % %     % stats2 is varargin{1}
% % % % %     % w is [a b c' c ab] weights, N rows
% % % % %
% % % % %
% % % % %     w = ones(N, npaths);
% % % % %
% % % % %     % we have stats structure with (maybe) weights
% % % % %     if isfield(stats2, 'w')
% % % % %         w = stats2.w;
% % % % %         w = w .* N;
% % % % %     end
% % % % %
% % % % %     w = w(whgood,:);
% % % % %
% % % % %     ncols = 4;
% % % % %     if npaths == 5, nrows = 1; else nrows = ceil(npaths ./ ncols); end
% % % % %
% % % % %     fh = create_figure('Slope_Plot', nrows, ncols);
% % % % %
% % % % %     %     fh = findobj('Tag', 'Slope_Plot');
% % % % %     %     if isempty(fh)
% % % % %     %         fh = tor_fig(nrows, ncols);
% % % % %     %         set(fh, 'Tag', 'Slope_Plot');
% % % % %     %         widen_figure(fh);
% % % % %     %     else
% % % % %     %         set(0, 'CurrentFigure', fh);
% % % % %     %         arrayfun(@(axh)cla(axh), findobj(get(fh, 'Children'), 'Type', 'axes'))
% % % % %     %         hold on;
% % % % %     %     end
% % % % %
% % % % %     subplot(nrows, ncols, 1); hold on;
% % % % %     plot_individual_slopes(abetas, X, w(:,1));
% % % % %     xlabel(vnames{1}), ylabel(vnames{3});
% % % % %     title('a: X->M');
% % % % %
% % % % %     subplot(nrows, ncols, 2); hold on;
% % % % %     plot_individual_slopes(bbetas, M, w(:,2));
% % % % %     xlabel(vnames{3}), ylabel(vnames{2});
% % % % %     title('b: M->Y controlling X');
% % % % %
% % % % %     subplot(nrows, ncols, 3); hold on;
% % % % %     plot_individual_slopes(cpbetas, X, w(:,3));
% % % % %     xlabel(vnames{1}), ylabel(vnames{2});
% % % % %     title('c'': X->Y controlling M');
% % % % %
% % % % %     subplot(nrows, ncols, 4); hold on;
% % % % %     plot_individual_slopes(cbetas, X, w(:,4));
% % % % %     xlabel(vnames{1}), ylabel(vnames{2});
% % % % %     title('c: X->Y');
% % % % %
% % % % %     n_addM = size(additionalM, 2);
% % % % %     %%% ***add plots for multiple mediators!
% % % % % end
% % % % %
% % % % %
% % % % % function plot_individual_slopes(betas, X, w)
% % % % %
% % % % %     % sort so that we plot from lowest to highest weights
% % % % %     [w, sorti] = sort(w);
% % % % %     betas = betas(:,sorti);
% % % % %     if iscell(X), X = X(sorti); else X = X(:,sorti); end
% % % % %
% % % % %     % get colors based on weights
% % % % %     N = size(betas, 2);
% % % % %     minwt = .4;     % make sure all lines are visible.
% % % % %     w = rescale_range(w, [minwt 1]);
% % % % %
% % % % %     % line widths: median split, top half gets 2.
% % % % %     linew = (w(:,1) >= median(w(:,1))) + 1;
% % % % %
% % % % %     w = 1 - w;  % for colors, 0 is black
% % % % %     colors = repmat([1 1 1], N, 1) .* repmat(w, 1, 3);
% % % % %     colors(colors < 0) = 0; % include to remove rounding error
% % % % %
% % % % %     %all_mx = zeros(N, 1);
% % % % %     for i = 1:N
% % % % %         if iscell(X), x = X{i}; else x = X(:,i); end
% % % % %
% % % % %         % 95% conf. interval for x
% % % % %         [nanvec, x] = nanremove(x);
% % % % %         mx = mean(x);
% % % % %
% % % % %         %all_mx(i) = mx;
% % % % %         s = std(x) * tinv(.975, length(x)-1);
% % % % %         x = [mx - s mx + s];
% % % % %         y = betas(2, i) + betas(1, i) * x;
% % % % %         plot(x, y, '-', 'Color', colors(i,:), 'LineWidth', linew(i));
% % % % %     end
% % % % %
% % % % %     % mean line
% % % % %     wmean = @(paths, w) ((w ./ sum(w))'*paths)';
% % % % %     means = wmean(betas', w); % should be slope then intercept
% % % % %     x = get(gca, 'XLim');
% % % % %     y = means(2) + means(1) * x;
% % % % %     plot(x, y, '--', 'Color', 'r', 'LineWidth', 3);
% % % % %
% % % % %     % SE lines (bootstrapped; multilevel)
% % % % %     xx = linspace(x(1), x(2), 50);
% % % % %
% % % % %     boots = 1000;
% % % % %     means =  bootstrp(boots, wmean, betas', w);
% % % % %     yy = zeros(boots, 50);
% % % % %     for i = 1:boots
% % % % %         yy(i,:) = means(i, 1) * xx + means(i, 2);
% % % % %     end
% % % % %
% % % % %     ub = prctile(yy, 95);
% % % % %     plot(xx, ub, 'r', 'LineWidth', 1);
% % % % %     lb = prctile(yy, 5);
% % % % %     plot(xx, lb, 'r', 'LineWidth', 1);
% % % % %
% % % % %     % xb = mean(x); s = resid. std. stat.se(2) is se of slope
% % % % %     %se_mean = (stat.se(2).^2 * (xx-xb).^2 + s.^2/length(x)) .^ .5;
% % % % %
% % % % %     drawnow
% % % % % end
% % % % %
% % % % %
% % % % % function rx = rescale_range(x, y)
% % % % %     % re-scale x to range of y
% % % % %     m = range(y)./range(x);
% % % % %
% % % % %     if isinf(m)
% % % % %         % no range/do not rescale
% % % % %         rx = x;
% % % % %     else
% % % % %         b = y(1) - m * x(1);
% % % % %         rx = m*x + b;
% % % % %     end
% % % % % end


% -------------------------------------------------------------------------
% Plot estimated impulse responses (HRFs) in latent model
% -------------------------------------------------------------------------
function plot_hrf_in_latent_model(stats1)

N = size(stats1.mean, 1);

nplots = ceil(sqrt(N));

fh = findobj('Tag', 'HRF_Plot');
if isempty(fh)
    fh = tor_fig;
    set(fh, 'Tag', 'HRF_Plot');
else
    set(0, 'CurrentFigure', fh);
    clf(fh);
end


for i = 1:N
    subplot(nplots, nplots, i);
    plot(stats1.hrf_xmy{i});
    if i == 1, legend(stats1.vnames), end
    axis auto
end
drawnow
end


% -------------------------------------------------------------------------
% Special figure position for wide plots
% -------------------------------------------------------------------------
function widen_figure(fh)
s = get(0, 'ScreenSize');
pix = max(s) .* .8;
pix = min(1000, pix);
set(fh, 'Position', [100 s(4)-100 pix pix ./ (2*1.67)]);
end

% -------------------------------------------------------------------------
% Save all figures
% -------------------------------------------------------------------------
function save_figures(verbose)

if nargin == 0, verbose = 1; end

fignames = {'Path_Diagram' 'Slope_Plot' 'Histogram_Plot' 'HRF_Plot' 'Mediation_Scatterplots', 'Individual Effects', 'cov(a,b)'};

for i = 1:length(fignames)
    name = fignames{i};
    h = findobj('Tag', name);
    if ~isempty(h)
        
        if length(h) > 1
            if verbose, disp(['Warning: More than one plot with tag ' name '. Saving highest-numbered figure.']); end
            h = max(h);
        end
        
        figure(h)
        scn_export_papersetup(500);
        
        if verbose
            fprintf('Saving: %s%s\n', name, '.fig');
        end
        saveas(h, name, 'fig');
        
        if verbose
            fprintf('Saving: %s%s\n', name, '.png');
        end
        saveas(h, name, 'png');
    else
        if ~strcmp(name, 'cov(a,b)')  % this won't exist for any single-level model, so avoid annoying warning.
            disp(['Cannot find figure with tag ' name ' to save.  Maybe figure was not created, or it was closed?']);
        end
    end
end

end



% -------------------------------------------------------------------------
% Print summary of mediation analysis results to screen
% -------------------------------------------------------------------------
function print_outcome(stats, stats1)

if nargin < 2, stats1 = []; end

fprintf('\n________________________________________\n')
if isfield(stats, 'bootsamples')
    fprintf('Final bootstrap samples: %3.0f\n', stats.bootsamples)
end
if isfield(stats, 'alphaaccept')
    fprintf('Average p-value tolerance (average max alpha): %3.4f\n', stats.alphaaccept)
end

% convergence values
if ~isempty(stats1)
    if isfield(stats1, 'isconverged') && isfield(stats1, 'mean')
        N = size(stats1.mean, 1);
        mysum = sum(stats1.isconverged);
        myperc = 100*mysum ./ N;
        fprintf('Number converged: %3.0f, %3.0f%%\n', mysum, myperc)
    end
end

% this would be true for one-tailed, but not 2-tailed
%Z = abs(norminv(stats.p)) .* sign(stats.mean);
% Changed 4/9/08
if isfield(stats, 'beta')
    % compatible with scn_stats_helper
    Z = abs(norminv(stats.p ./ 2)) .* sign(stats.beta);
else
    % compatible with bootstrap/signperm
    Z = abs(norminv(stats.p ./ 2)) .* sign(stats.mean);
end

fprintf('\n%s\n', stats.analysisname)
fprintf('\t')
fprintf('%s\t', stats.inputOptions.names{:})
fprintf('\n')
print_line('Coeff', stats.mean(1, :));
print_line('STE', stats.ste(1, :))
if isfield(stats, 't')
    print_line('t (~N)', stats.t(1, :))
else
    print_line('t (~N)', stats.mean(1, :) ./stats.ste(1, :))
end
print_line('Z', Z(1, :))
print_line('p', stats.p(1, :), 4)
fprintf('\n')


fprintf('________________________________________\n')

if size(stats.p, 1) > 1 && isfield(stats, 'beta')
    % we have OLS 2nd level model (***beta is carryover)
    
    for i = 2:size(stats.p, 1)
        % Each effect: Intercept, then 2nd-level covs
        
        fprintf('\n________________________________________\n')
        fprintf('Second Level Moderator')
        
        fprintf('\n%s\n', stats.analysisname)
        fprintf('\t')
        fprintf('%s\t', stats.inputOptions.names{:})
        fprintf('\n')
        print_line('Coeff', stats.beta(i, :));
        print_line('STE', stats.ste(i, :))
        if isfield(stats, 't')
            print_line('t (~N)', stats.t(i, :))
        else
            print_line('t (~N)', stats.mean(i, :) ./stats.ste(i, :))
        end
        print_line('Z', Z(i, :))
        print_line('p', stats.p(i, :), 4)
        fprintf('\n')
        
        
        fprintf('________________________________________\n')
        
    end
end
end

function print_line(hdr, data, dec)
if nargin == 2
    dec = num2str(2);
else
    dec = num2str(dec);
end

fprintf('%s\t', hdr)
if iscell(data)
    eval(['fprintf(''%3.' dec 'f\t'', data{:})']);
else
    eval(['fprintf(''%3.' dec 'f\t'', data)']);
end
fprintf('\n');
end


function mystruct = get_stats_names(mystruct)
% add standard path names onto stats strucutre
mystruct.names = {'b1' 'b2' 'b3' 'c''' 'c' 'b1b2b3'};
mystruct.inputOptions.names = mystruct.names;

end

