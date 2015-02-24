% function [paths, toplevelstats, 1stlevelstats] = mediation(X, Y, M, [stats], [plots], [other optional args])
%
% Tor Wager, March 2006
%
% X, Y, M can be
% 1) vectors of observations
% 2) matrices of N columns (observations x N subjects)
% 3) cell arrays of length N (each cell is vector of obs. for one
% subject)
%
%
% columns of paths:
% -------------------------------------------------------------------------
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)
%
% 6 a2  X -> M2, if additional M's are entered
% 7 b2  M2 -> Y, if additional M's are entered
% 8 ab2 a2 * b2
%
% Optional arguments:
% -------------------------------------------------------------------------
% Data input
% 'M'               Followed by additional mediator(s), in same input
%                   format as X, Y, and M (vector, matrix for multilevel,
%                   or cells with vectors for multilevel)
% 'L2M'             Followed by second-level moderators (n subjects x k
%                   moderators)
%                   * NOTE: works only with OLS at 2nd level now.
%
% Output control
% 'names'           Followed by 3-cell vector of names for X, Y, and M
% 'verbose'         Verbose output
% 'noverbose'       Suppress output (default)
% 'plots'           Plot histograms of coefficients and regression fits
% 'dosave'          Save figures after plotting
%
% Bootstrapping/permutation options
% 'boottop'         In a single-level analysis, do bootstrap statistics (aka 'boot1')
% 'bootstrapfirst'  In a multi-level analysis, bootstrap first level stats
% 'bootsamples'     Number of bootstrap samples (default 1000)
% 'signperm'        Do sign permutation test at 2nd level (alternative to
%                   bootstrap)
%
% Multilevel options
% 'hierarchical'    multi-level model with weighing based on 1st level variance
% 'summarystats'    use summary-statistics approach with multi-level data;
% no weighting
%
% Timeseries modeling options
% 'arorder'         followed by order of AR(p) model for timeseries data
%                   used at 1st level only in multi-level analysis
% 'shiftrange'      followed by range in elements to shift first-level vectors, e.g., [-6 6]
%                   appropriate if first-level data are timeseries with
%                   unknown relative lags.  Y is the reference region, and
%                   others are shifted relative to Y.
% 'latent'          mediation is done on latent varibles after deconvolving
%                   hrf estimates from each region.  hrf estimates are
%                   constrained but vary somewhat from region to region
%                   they are returned in stats1 output
%                   this option is a mutually exclusive alternative to
%                   'shiftrange'
% Other options
% 'covs'            Followed by set of covariates
%                   controlled for in all regressions 
%                   (*now: 1st lev covariates in multi-level model do not work with SHIFT and LATENT models)
%
% 'persistent_perms' Re-use same permutations. Much faster in sign test,
%                   but not a good idea for independence across replications.
%                   May be OK for repeated estimates across brain image,
%                   but creates problems in simulation because replications
%                   are not all independent due to stable choice of perms.
%                   Default is 0 (no).
%
% Programmers' notes:
% 7/10/07, Tor added multiple mediator option
% need to update path slope plots, robust, other parts of
% mediation_path_coefficients, path plot,
%
% Examples:
% -------------------------------------------------------------------------
%
% X = rand(50, 1); Y = X + rand(50, 1); M = X + rand(50, 1);
% [paths, stats1, stats2] = mediation(X, Y, M, 'boottop', 'stats', 'plots');
%
% [paths, stats2] = mediation(X, Y, M, 'plots');
%
% [paths, stats2] = mediation(data.X, data.Y, data.M, 'plots', 'robust', 'verbose');
%
% [paths, stats2] = mediation(data.X, data.Y, data.M, 'plots', 'robust', 'verbose', 'names', {'Left V1' 'Left M1' 'Right V1'});
%
% [paths, stats2] = mediation(X, Y, M, 'plots', 'boot', 'shift', [-6 6]);
%
% [paths, stats] = mediation(x, y, m, 'covs', y2(:,4:10))
%
% Second-level moderator example
% [paths, stats2] = mediation(speechvsbase, hr, vmpfc.timeseries, 'plots',
% 'verbose', 'names', {'Speech-Base' 'HR increase' 'vmPFC'}, 'L2M', hrinc);
%
% Multiple mediator example:
%  x = randn(100, 1); m = x + randn(100, 1); m2 = x + randn(100, 1); y = m2 + randn(100, 1);
% [paths, stats] = mediation(x, y, m, 'M', m2, 'boot', 'plots', 'verbose', 'bootsamples', 10000);
%
% See also mediation_brain.m

% Programmer's notes:
% 6/8/08 : Worked on wls and empirical bayes from Raudenbush and Bryk book.
% This is needed to give us optimal weighting, and needs to be developed.
% 
% 
%

function [paths, varargout] = mediation(X, Y, M, varargin)
    % -------------------------------------------------------------------------
    % Setup inputs, print info to screen if verbose
    % -------------------------------------------------------------------------
    global mediation_covariates  % used in mediationfun and mediation_path_coefficients; empty = no covs

    varargout = cell(1, nargout - 1);

    [doplots, boottop, dorobust, verbose, vnames, N, wistats, bootsamples, dobootfirstlevel, domultilev, ...
        arorder, targetu, shiftrange, dolatent, whpvals_for_boot, dosave, dosignperm, persistent_perms, additionalM, ...
        X_2ndlevel] = ...
        setup_inputs(Y, varargin);

    if verbose, totalt = clock; end

    npaths = 5 + 3 * size(additionalM, 2);

    paths = NaN .* zeros(N, npaths);
    sterrs = Inf .* ones(N, npaths);
    abetas = cell(1, N);
    bbetas = cell(1, N);
    cpbetas = cell(1, N);
    cbetas = cell(1, N);

    residm = cell(1, N);
    residy = cell(1, N);
    residy2 = cell(1, N);

    n = zeros(N, 1);

    nrm = zeros(N, 1);

    % First level/Single level paths: set handle for mediation function to run
    if dorobust
        mediationfun = @fast_robust_ab;
    else
        mediationfun = @fast_ols_ab; 
    end

    % Weighted mean function:
    % Used in getstats to get mean bootstrap coefficients
    % Second level function for getting weights across subjects
    % Faster than looping, and faster than using mean() for equal weights
    % Also makes code simpler below.
    wmean = @(paths, w) diag(w'*paths)';

    if N > 1
        w = ones(N) ./ N;  % start with weights equal
    else
        w = ones(size(X, 1)) ./ size(X, 1);
    end

    if dolatent
        if ~isempty(mediation_covariates), warning('Covariates not implemented yet for latent model!'), end
        hrfparams = cell(1, N);
        hrf_xmy = cell(1, N);
    elseif any(shiftrange)
        if ~isempty(mediation_covariates), warning('Covariates not implemented yet for shift model!'), end
        totalsse = zeros(N, 1);
        delays = zeros(N, 2);
        isconverged = zeros(N, 1);
    end

    % =========================================================================
    %
    % * Single-level model, run for each replication N
    %
    % =========================================================================

    % make nicer code later for lines below***
    MCOVS = mediation_covariates;
    
    for i = 1:N

        if verbose && N > 1, fprintf('\b\b\b%03d', i); end

        if iscell(X), x = X{i}; else x = X(:,i); end
        if iscell(Y), y = Y{i}; else y = Y(:,i); end
        if iscell(M), m = M{i}; else m = M(:,i); end

        if ~isempty(additionalM)
            if iscell(additionalM), m = [m additionalM{i}]; else m = [m additionalM(:,i)]; end
        end

         if ~isempty(MCOVS)
             % Cell for multi-level, matrix for single-level
            if iscell(MCOVS), mediation_covariates = MCOVS{i}; else mediation_covariates = MCOVS; end
         end
        
        % --------------------------------------------
        % Remove NaNs and Check for bad data
        % --------------------------------------------
        % Remove NaNs from all vars, casewise
        [nanvec, x, y, m, mediation_covariates] = nanremove(x, y, m, mediation_covariates);

        %****NOTE: Check to make sure no exact duplicates in x , y , m
        % ***check to make sure length valid obs > 1
        
        if isempty(x) || all(x == x(1)) || all(y == y(1)) || all(m(:,1) == m(1))
            % No valid data; skip this subject
            if verbose, warning(['No valid data for subject ' num2str(i) ': Skipping.']); end

            stats = get_ols_stats(paths(i,:) , sterrs(i,:), n(i), additionalM);  % setup dummy stats structure
            if N > 1
                if i == 1, wistats = stats; wistats.mediation_covariates = mediation_covariates;  end
                wistats = collect_within_stats(wistats, stats, x, y, m, i);
            end

            continue
        end


        % --------------------------------------------
        % Compute path coefficients (and standard errors in some cases)
        % If bootstrapping is 'on', skips standard errors
        % If bootstrapping is 'off', returns std. errors, intercept, and n
        % needed for OLS stats
        % --------------------------------------------

        if dolatent
            % ***need to update to cells for multiple mediators
            [paths(i,:), abetas{i}, bbetas{i}, cpbetas{i}, cbetas{i}, sterrs(i,:), intcpt, n(i), residm{i}, residy{i}, residy2{i}, ...
                totalsse(i), hrfparams{i}, hrf_xmy{i}, isconverged(i)] = ...
                mediation_latent(x, y, m, 'ga', domultilev, dorobust, boottop);

        elseif any(shiftrange)
            [paths(i,:), abetas{i}, bbetas{i}, cpbetas{i}, cbetas{i}, sterrs(i,:), intcpt, n(i), residm{i}, residy{i}, residy2{i}, ...
                totalsse(i), delays(i,:), isconverged(i)] = ...
                mediation_shift(x, y, m, shiftrange, 'ga', domultilev, dorobust, boottop);
            
        else
            [paths(i,:), abetas{i}, bbetas{i}, cpbetas{i}, cbetas{i}, sterrs(i,:), intcpt, n(i), residm{i}, residy{i}, residy2{i}, ...
                Vxm{i}, Vmy{i}, Vxy{i}] = ...
                mediation_path_coefficients(x, y, m, domultilev, dorobust, boottop);
        end

        % Multilevel only: save std. deviations for getting standardized paths
        if N > 1, nrm(i, 1) = std(x) ./ std(y); end


        % --------------------------------------------
        % Bootstrap a*b [optional]
        % This will be done for single-level if boottop is specified
        % or multi-level if dobootfirst is specified
        % --------------------------------------------

        boot_this_pass = (boottop && N == 1) || (dobootfirstlevel && N > 1);
        if boot_this_pass
            if verbose, fprintf('\nBootstrapping: '), t1 = clock; end

            if isempty(mediation_covariates)
                bootpaths = bootstrp(bootsamples, mediationfun, x, y, m, intcpt);
            else
                bootpaths = bootstrp(bootsamples, mediationfun, x, y, m, intcpt, mediation_covariates);
            end

            stats = getstats(bootpaths, additionalM);
            stats.vnames = vnames;

            % check how many Boot samples we need, and get more if necessary
            add_boot_samples_needed();
            stats; 
            bootpaths;

            % bias correction for final bootstrap samples
            %[p, z] = bootbca_pval(testvalue, bootfun, bstat, stat, [x], [other inputs to bootfun])
            if isempty(mediation_covariates)
                %[stats.p, stats.z] = bootbca_pval(0, mediationfun, bootpaths, mediationfun(x, y, m, intcpt), x, y, m, intcpt);

                [dummy, dummy, stats.z, stats.p] = bootbca_pval_onetail(0, mediationfun, bootpaths, mediationfun(x, y, m, intcpt), x, y, m, intcpt);


            else
                %[stats.p, stats.z] = bootbca_pval(0, mediationfun, bootpaths, mediationfun(x, y, m, intcpt, mediation_covariates), x, y, m, intcpt, mediation_covariates);
                [dummy, dummy, stats.z, stats.p] = bootbca_pval_onetail(0, mediationfun, bootpaths, mediationfun(x, y, m, intcpt, mediation_covariates), x, y, m, intcpt, mediation_covariates);

            end
            stats.biascorrect = 'BCa bias corrected';

            if verbose, fprintf('Done in %3.0f (s) \n', etime(clock, t1)); end
        else
            % we still must collect stats for weighting in multi-level model
            stats = get_ols_stats(paths(i,:), sterrs(i,:), n(i), additionalM);
        end

        if exist('stats', 'var')
            stats.residm = residm;
            stats.residy = residy;
            stats.residy2 = residy2;
        end


        % --------------------------------------------
        % multi-subject (multilevel)
        % collect within-subject stats; assumes independence of trials
        % --------------------------------------------
        if N > 1
            if i == 1, wistats = stats; end
            wistats = collect_within_stats(wistats, stats, x, y, m, i);
            wistats.mediation_covariates = mediation_covariates;
        elseif exist('stats', 'var')
            % --------------------------------------------
            % single-level model
            % save stats stuff for output
            % --------------------------------------------
            % single-level model, return stats structure if we have one
            % either from bootstrapping or OLS

            % save input options and stuff; borrow 2nd-level code for
            % convenience
            %stats2 = stats;
            stats = add_info_to_stats(stats);
            stats.inputOptions.names = stats.names;  % make consistent with glmfit_general
            stats.analysisname = 'Single-level model';
            %stats = stats2;
            %clear stats2
            varargout{1} = stats;
        end

        % Plots and printing
        if doplots && boot_this_pass %boottop && i == 1
            plot_hists(bootpaths, vnames);
        elseif doplots && i == 1
            % Do nothing -- do not do all individual plots in multi-level
            % model
            % disp('Plotting for non-bootstrapped results not implemented yet.')
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
    % * Second-level model, if N > 1
    %
    % --------------------------------------------------

    if N > 1 && all(isnan(paths(:)))  %~any(wistats.df < 1)
        % NO valid data; exit gracefully without continuing computation

        if verbose, disp('No valid data for any subjects!'); end

        stats2 = getstats(paths, additionalM);
        stats2 = add_info_to_stats(stats2);

        varargout{1} = stats2;
        varargout{2} = wistats;

        return
    end


    % --------------------------------------------
    % If N > 1, then get group stats (2nd level)
    % --------------------------------------------
    if N > 1
        if boottop
            means = [];
            second_level_bootstrap();

        elseif dosignperm
            means = [];
            second_level_permtest();

        else
            % OLS case (2nd-level regression)
            [b, s2between_ols, stats2] = scn_stats_helper_functions('gls', paths, w, X_2ndlevel);

            stats2.std = stats2.ste .* sqrt(N);
            
            % weight, and do it again with weights
            get_weights_based_on_varcomponents()
            [b, s2between_ols, stats2] = scn_stats_helper_functions('gls', paths, w, X_2ndlevel);
            
            stats2.name = '2nd level statistics';
            stats2.names = {'a' 'b' 'c''' 'c' 'ab'};
            % dfe should be all the same, except in autocorr. case
            stats2.Z = abs(norminv(stats.p ./ 2)) .* sign(stats.mean);
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

    if doplots && N > 1
        % Multi-level ----------------------------------
        
        if dorobust
            disp('No robust slope plots yet.')
        else
            % plots of slopes for each subject
            mediation_plots(stats2, 'slopes', varargin);
            %plot_slopes(abetas, bbetas, cbetas, cpbetas, X, M, additionalM, vnames, stats2);
        end

        try mediation_path_diagram(stats2);
        catch disp('Error in mediation_path_diagram');
        end

        if dolatent
            try plot_hrf_in_latent_model(wistats);
            catch disp('Error in plot_hrf_in_latent_model');
            end
        end
        
        % Plot individual effects
        create_figure('Individual Effects'); barplot_columns(paths, 'Individual Effects', [], 'nofig', 'plotout','number');
        set(gca, 'XTickLabel', stats.names);
        ylabel('Beta (slope) values');
        
    elseif doplots && N == 1
        % Single-level ----------------------------------
        
        try mediation_path_diagram(stats);
        catch disp('Error in mediation_path_diagram');
        end

        mediation_scatterplots(stats);
    end

    if verbose && exist('stats2', 'var')
        if exist('wistats', 'var') && ~isempty(wistats)
            print_outcome(stats2, wistats)
        else
            print_outcome(stats2)
        end
    elseif verbose && N == 1 && exist('stats', 'var')
        print_outcome(stats)
    end

    if verbose, fprintf('Total time: %3.0f s\n', etime(clock, totalt));  end


    % error checking:
    % estab = stats2.mean(1) * stats2.mean(2) + 2 * cov(means(:,1:2)) % approx; cov should be weighted?

    if dosave
        if doplots
            save_figures(verbose);
        end

        if verbose
            diary off
            fprintf('All finished: Closed file Mediation_Output.txt\n')
        end
    end



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
    function second_level_permtest()
        persistent permsign

        if ~persistent_perms
            permsign = [];
        end

        if verbose, fprintf('Nonparametric test with %3.0f permutations...', bootsamples); t12 = clock; end

        % start with weights all equal whether multilevel or not

        [p, Z, xbar, permsign, means] = permute_signtest(paths, bootsamples, w, permsign);

        stats2 = [];
        stats2.name = '2nd level statistics';
        stats2.names = {'a' 'b' 'c''' 'c' 'ab'};
        stats2.mean = xbar;
        stats2.ste = std(means);
        stats2.std = std(paths);
        stats2.Z = Z;
        stats2.p = p;


        % get weights, if multilevel
        % returns w and stats2.w, etc.
        if domultilev
            if verbose, fprintf('Re-permuting with multi-level weights.\n'), end
            get_weights_based_on_varcomponents();
            % do full test based on new weights
            [p, Z, xbar, permsign, means] = permute_signtest(paths, bootsamples, w, permsign);
            stats2.mean = xbar;

            % *** could replace with weighted std; also, check this
            stats2.ste = std(means);
            stats2.std = std(paths);
            stats2.Z = Z;
            stats2.p = p;
        end


        if verbose, fprintf(' Done in %3.0f s \n', etime(clock, t12)); end
    end


    % -------------------------------------------------------------------------
    % Bootstrap test for 2nd level
    % -------------------------------------------------------------------------
    function second_level_bootstrap()
        if verbose, fprintf('Bootstrapping %3.0f samples...', bootsamples); end

        % initalize random number generator to new values; bootstrp uses this
        rand('twister', sum(100*clock))


        whnan = nanremove(paths, w);
        whgood = ~whnan;
        if verbose && any(whnan), fprintf('Warning! %3.0f 2nd-level observations are invalid/have missing data.\n', sum(whnan)); end

        % start with weights all equal whether multilevel or not
        if sum(whgood) < 2
            means = zeros(bootsamples, npaths);
        else
            means = bootstrp(bootsamples, wmean, paths(whgood,:), w(whgood,:));
        end
        stats2 = getstats(means, additionalM); 

        stats2.name = '2nd level statistics';

        % std is needed for weighted mean option below.
        stats2.std = stats2.ste .* sqrt(N);

        % Get boot samples needed
        % Add samples as needed to ensure p-value is reasonably accurate
        % Be: Boot samples needed so that Expected (average) p-value
        % error due to limited-sample bootstrap is targetu*100 %

        [Be, alphaaccept] = get_boot_samples_needed(stats2.p, whpvals_for_boot, targetu, bootsamples, verbose); % uses whpvals_for_boot, returns Be

        if verbose && boottop && (Be > 0 || domultilev), fprintf('Bootstrapping...'); end
        if verbose, t12 = clock; end

        % get weights, if multilevel
        % returns w and stats2.w, etc.
        if domultilev
            get_weights_based_on_varcomponents();

            % get full bootstrap sample based on new weights
            if sum(whgood) > 2
                means = bootstrp(bootsamples + Be, wmean, paths(whgood,:), w(whgood,:));
            end
        else
            % no weighting; summary stats.
            % add boot samples only if necessary
            if Be > 0 && sum(whgood) > 2
                means = [means; bootstrp(Be, wmean, paths(whgood,:), w(whgood,:))];
            end
        end

        stats2 = getstats(means, additionalM, stats2); % last input preserves existing info

        % use original weighted mean, not mean of bootstrap samples
        stats2.mean = wmean(paths(whgood,:), w(whgood,:));

        % std is needed for weighted mean option below.
        stats2.std = stats2.ste .* sqrt( sum(whgood) );

        % bias correction for final bootstrap samples
        %[p, z] = bootbca_pval(testvalue, bootfun, bstat, stat, [x], [other inputs to bootfun])
        stats2.prctilep = stats2.p; % percentile method, biased
        %[stats2.p, stats2.z] = bootbca_pval(0, wmean, means, wmean(paths, w), paths, w);
        [dummy, dummy, stats2.z, stats2.p] = bootbca_pval_onetail(0, wmean, means, wmean(paths(whgood,:), w(whgood,:)), paths(whgood,:), w(whgood,:));

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

        if Be > 0
            if verbose, fprintf(' Adding %3.0f samples ', Be), t12 = clock; end
            if dorobust
                if isempty(mediation_covariates)
                    bootpaths = [bootpaths; bootstrp(Be, @fast_robust_ab, x, y, m)];   % twice as fast as long version
                else
                    bootpaths = [bootpaths; bootstrp(Be, @fast_robust_ab, x, y, m, mediation_covariates)];
                end

            else
                if isempty(mediation_covariates)
                    bootpaths = [bootpaths; bootstrp(Be, @fast_ols_ab, x, y, m, intcpt)];
                else
                    bootpaths = [bootpaths; bootstrp(Be, @fast_ols_ab, x, y, m, intcpt, mediation_covariates)];
                end
            end

            stats = getstats(bootpaths, additionalM);
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
        % SHOULD THIS BE STD??

        stats2.total_var_est = (stats2.std.^2);

        %******************** Edits by ML ********************
        % Use SST = SSB + SSW.
        % Once you have calculated SSB, you can find MSB
        % N*(nobs-1) was changed for diff. n for each subject (TW)
        totaln = sum(n-1);
% % %         SST = totaln * (stats2.std.^2);
% % %         SSW = totaln * stats2.wi_var_est;
% % %         stats2.btwn_var_est = (SST-SSW)/(N-1);
% % % 
% % %         stats2.btwn_var_est = max(0, stats2.btwn_var_est);
% % %         %*****************************************************
% % % 
% % %         %stats2.std_beta = stats2.ste .* wistats.avg_sx ./ wistats.avg_sy;
% % %         %stats2.btwn_prop = stats2.std_beta ./ (stats2.std_beta + mean(wistats.ste));
% % %         stats2.btwn_std = stats2.btwn_var_est .^ .5;
% % % 
% % %         b_plus_w = ( wistats.std.^2 + repmat(stats2.btwn_var_est, N, 1) );
% % %         b_plus_w( b_plus_w < 10 * eps | isnan(b_plus_w) ) = Inf;
% % % 
% % %         stats2.w = 1  ./ b_plus_w;
% % %         % normalize by sum so weights sum to 1
% % %         stats2.w = stats2.w ./ repmat(sum(stats2.w), N, 1);
% % % 
% % %         w = stats2.w;
        
        % get weights : work in progress. Sep. for each equation.
        
        wh_omit = false(1, length(Vxm));
        
        for i = 1:length(Vxm)
            if isempty(Vxm{i}) | isempty(Vmy{i}) | isempty(Vxy{i}), wh_omit(i) = 1;
            end
        end
        Vxm(wh_omit) = [];
        Vmy(wh_omit) = [];
        Vxy(wh_omit) = [];
           
        mycol = 1;
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(cat(2, abetas{:}), Vxm);
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_std(mycol) = sqrt(sigma2_b_est(1));
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        
        mycol = [2 3];
        bb = cat(2, bbetas{:});
        bb = [bb(1, :); cat(2, cpbetas{:})];
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(bb, Vmy);
        stats2.wls_mean(mycol) = gam_hat(1:2);  stats2.wls_t(mycol) = gam_t(1:2); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1:2);
        stats2.btwn_var_est(mycol) = sigma2_b_est(1:2);
        
        mycol = 4;
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(cat(2, cbetas{:}), Vxy);
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        
        mycol = 5;
        % need to figure out ab var...
        wh_include = find(~isnan(paths(:, 5)));
        for i = 1:length(wh_include), Vab{i} = sterrs(wh_include(i), 5) .^ 2; end   % get V within subj for each ab path
        
        [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(paths(wh_include, 5)', Vab);
        
        % ****NOTE: fix for some subjects have missing data.
        
        stats2.wls_mean(mycol) = gam_hat(1);  stats2.wls_t(mycol) = gam_t(1); % 1st col is what we want
        stats2.w(:, mycol) = subject_w(:, 1);
        stats2.btwn_var_est(mycol) = sigma2_b_est(1);
        
        % make sure we deal with missing data subjects
        if length(wh_include) < N
            w = zeros(size(paths));
            w(wh_include, :) = stats2.w;
            stats2.w = w;
        end
        
        w = stats2.w;
    end




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

        if dolatent
            wistats.totalsse = totalsse;
            wistats.hrfparams = hrfparams;
            wistats.hrf_xmy = hrf_xmy;
            wistats.isconverged = isconverged;
        elseif any(shiftrange)
            wistats.totalsse = totalsse;
            wistats.delays = delays;
            wistats.isconverged = isconverged;
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

        inputOptions = struct('N', N, 'npaths', npaths, 'initial_bootsamples', bootsamples, 'arorder', arorder, 'targetu', targetu);

        inputOptions.vnames = vnames; % do we need this?
        inputOptions.names = stats_struct.names;

        % save data
        inputOptions.X = X;
        inputOptions.Y = Y;
        inputOptions.M = M;
        inputOptions.additionalM = additionalM;
        inputOptions.mediation_covariates = mediation_covariates;
        inputOptions.X_2ndlevel = X_2ndlevel;

        % save options
        inputOptions.bootfirstlevel = ynstr{dobootfirstlevel + 1};
        inputOptions.robust = ynstr{dorobust + 1};
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
        stats_struct.abetas = abetas;
        stats_struct.bbetas = bbetas;
        stats_struct.cpbetas = cpbetas;
        stats_struct.cbetas = cbetas;

        stats_struct.inputOptions = inputOptions;
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
% Setup inputs, print info to screen if verbose
% -------------------------------------------------------------------------
function [doplots, boottop, dorobust, verbose, vnames, N, wistats, bootsamples, dobootfirstlevel, domultilev, ...
        arorder, targetu, shiftrange, dolatent, whpvals_for_boot, dosave, dosignperm, persistent_perms, ...
        additionalM, X_2ndlevel] = ...
        setup_inputs(Y, varargin)

    varargin = varargin{1}; % do this because we passed in all varargin args from parent into varargin cell

    % Defaults

    targetu = .20;              % proportion contribution of boot procedure to p-value
    doplots = 0;                % make plots
    boottop = 0;                  % bootstrap instead of OLS
    dobootfirstlevel = 0;       % bootstrap first level in multi-level analysis
    dorobust = 0;               % robust IRLS
    verbose = 0;                % verbose output
    vnames = {'X' 'Y' 'M'};     % variable names
    bootsamples = 1000;         % initial bootstrap samples
    whpvals_for_boot = [3 5];   % indices of p-values, the min of which is used to determine boot samples needed
    % lower p-values require more boot samples
    % for the p-vals to be meaningful.
    domultilev = 1;             % multi-level analysis (if N replications > 1)
    arorder = 0;                % timeseries AR(p) model order at first level
    shiftrange = [];            % min and max to shift timeseries by, in elements (or [] for no shift)
    dolatent = 0;               % latent HRF model

    dosave = 0;                 % save figures at end

    dosignperm = 0;             % sign permutation (alt. to bootstrap)
    persistent_perms = 0;       % keep same permutation matrix across repeated calls to mediation

    additionalM = [];           % additional mediators, entered in same format as other inputs

    X_2ndlevel = [];
    
    global mediation_covariates
    mediation_covariates = [];

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'plots', doplots = 1;
                case {'dosave', 'save', 'saveplots'}, dosave = 1;

                case {'boottop', 'boot1' 'boot'}, boottop = 1;
                case 'robust', dorobust = 1;
                case 'verbose', verbose = 1;
                case 'noverbose', verbose = 0;

                case {'ar', 'arorder' 'armodel'}, arorder = varargin{i+1};
                case {'shift', 'shiftrange'}, shiftrange = varargin{i+1};
                case {'latent', 'dolatent'}, dolatent = 1;

                case 'names', vnames = varargin{i+1};
                case 'bootsamples', bootsamples = varargin{i+1};
                case {'pvals', 'whpvals_for_boot'}, whpvals_for_boot = varargin{i+1};

                case {'m', 'M'}, additionalM = varargin{i+1};

                    % specific for multi-level model
                case 'bootstrapfirst', dobootfirstlevel = 1;
                case {'multilevel', 'hierarchical'}, domultilev = 1;
                case 'summarystats', domultilev = 0;
                case {'dosignperm', 'signperm'}, dosignperm = 1;

                case {'persistent_perms', 'persistent'}, persistent_perms = 1;

                case {'covs', 'covariates', 'mediation_covariates'}, mediation_covariates = varargin{i+1};

                case {'X2','L2moderators', 'L2M'}, X_2ndlevel = varargin{i+1};
            end
        end
    end

    % start diary, if appropriate
    if dosave
        if verbose
            fprintf('Saving output in file Mediation_Output.txt\n')
        end
        diary Mediation_Output.txt

    end

    % now done separately for each subject, because there could be nans,
    % diff. n for each subject, etc.
    %[nobs, N] = size(Y);
    if iscell(Y), N = length(Y);
    else N = size(Y, 2);
    end
    
    wistats = [];


    % if N is 1, then only one level; if N > 1, first level of multi-level
    % analysis.
    if N == 1, domultilev = 0; end

    % 2nd-level design matrix
    % ***NOTE: used only in OLS option now (not w/boot or signperm)
    % ***so 2nd level moderators only used in OLS
    if iscell(X_2ndlevel)
        error('2nd-level covariates must be a matrix of N subjects by k predictors');
    end
    X_2ndlevel = [ones(N, 1) scale(X_2ndlevel, 1)];

    
    % if additional mediators, make sure we have names
    if ~isempty(additionalM)
        for i = 1 : 3 + size(additionalM, 2) - length(vnames)
            vnames{3 + i} = ['AddM' num2str(i)];
        end
    end

    if boottop && dosignperm
        error('Bootstrapping and sign permutation cannot both be requested.');
    end

    if dorobust
        warning('off', 'stats:statrobustfit:IterationLimit');
        if verbose, disp('Note: Turning off iteration limit warning for robustfit.'); end
    end

    if verbose
        nobs_tmp = size(Y, 1);
        fprintf('Mediation analysis\n\nObservations: %3.0f, Replications: %3.0f\n', nobs_tmp, N);
        fprintf('Predictor (X): %s, Outcome (Y): %s: Mediator (M): %s\n', vnames{1}, vnames{2}, vnames{3});

        if ~isempty(additionalM)
            fprintf('Additional mediators: ')
            fprintf(repmat('%s  ', 1, size(additionalM, 2)), vnames{4:end})
            fprintf(1, '\n');
        end

        nms = {'No' 'Yes'};

        fprintf('\nCovariates: %s', nms{~isempty(mediation_covariates) + 1})
        if ~isempty(mediation_covariates)
            fprintf(', %3.0f columns', size(mediation_covariates, 2));
        end
        fprintf('\n')

        if N > 1
            fprintf('\nMulti-level analysis.\n')
            fprintf('Options:\n\tMultilevel weights: %s\n\tPlots: %s\n\tBootstrap: %s\n\tSign perm: %s\n\tRobust: %s\n\tBootstrap 1st level: %s\n', ...
                nms{domultilev+1}, nms{doplots+1}, nms{boottop+1}, nms{dosignperm+1}, nms{dorobust+1}, nms{dobootfirstlevel+1});
        else
            fprintf('\nSingle-level analysis.\n')
            fprintf('Options:\n\tPlots: %s\n\tBootstrap: %s\n\tRobust: %s\n', ...
                nms{doplots+1}, nms{boottop+1}, nms{dorobust+1});
        end

        if arorder, fprintf('\tAR(%d) model for timeseries data at first level.\n', arorder); end

        if boottop, fprintf('\tBootstrap or sign perm samples: %d\n', bootsamples); end

        if dolatent, fprintf('\tLatent HRF model: Yes');
        elseif any(shiftrange), fprintf('\tCross-lagged model with bounds: %3.1f to %3.1f\n', shiftrange(1), shiftrange(2));
        end

        if N > 1, fprintf('\nGetting estimates for replications: %03d', 0); end
    end



end





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
function stats = getstats(bp, additionalM, varargin)
    
    if ~isempty(varargin), stats = varargin{1}; end
    
    stats.names = {'a' 'b' 'c''' 'c' 'ab'};

    for i = 1:size(additionalM, 2)
        stats.names = [stats.names {['a_' num2str(i+1)]} {['b_' num2str(i+1)]} {['ab_' num2str(i+1)]}];
    end

    stats.mean = mean(bp);
    stats.ste = std(bp);
    %     stats.ci95 = [prctile(bp, 97.5); prctile(bp, 2.5)];
    %     stats.ci90 = [prctile(bp, 95); prctile(bp, 5)];
    %     stats.ci99 = [prctile(bp, 99.5); prctile(bp, .5)];
    stats.p = 2.*(min(sum(bp<=0), sum(bp>=0)) ./ size(bp, 1));

    % avoid exactly-zero p-values
    % replace with 1/# bootstrap samples
    stats.p = max(stats.p, 1./size(bp, 1));
end

% -------------------------------------------------------------------------
% Get statistic structure from OLS regression, including p-values and conf. intervals
% -------------------------------------------------------------------------
function stats = get_ols_stats(bp, sterrs, n, additionalM)
    stats.names = {'a' 'b' 'c''' 'c' 'ab'};

    for i = 1:size(additionalM, 2)
        stats.names = [stats.names {['a_' num2str(i+1)]} {['b_' num2str(i+1)]} {['ab_' num2str(i+1)]}];
    end

    stats.mean = bp;
    stats.ste = sterrs;
    stats.t = bp ./ sterrs;
    stats.df = [n-1 n-2 n-2 n-1 n-2];   % not sure about last one ***check***
    stats.p = min(1, (2 .* (1 - tcdf(abs(stats.t), stats.df))));

    % don't do to save time
    % stats.ci95 = [prctile(bp, 97.5); prctile(bp, 2.5)];
    % stats.ci90 = [prctile(bp, 95); prctile(bp, 5)];
    % stats.ci99 = [prctile(bp, 99.5); prctile(bp, .5)];
end


% -------------------------------------------------------------------------
% Collect structure of within-subject statistics
% -------------------------------------------------------------------------
function wistats = collect_within_stats(wistats, stats, x, y, m, i)
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
    if isempty(m) || all(m == m(1)), wistats.sm(i, 1) = Inf; else wistats.sm(i, 1) = std(m); end

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
function [paths] = fast_robust_ab(x, y, m, varargin)
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

    tmp = bootpaths(:,3:4);  % c ones
    tmp = tmp(:);
    nbins = max(10, round(length(tmp) ./ 100));
    nbins = min(nbins, 100);

    [h, xx] = hist(tmp, nbins);

    myperc = mean(xx) .* .2;
	cxlimit = [min(min(xx)-myperc, 0 - myperc) max(xx)+myperc];
    a = bootpaths(:,1);
    b = bootpaths(:,2);
    cp = bootpaths(:,3);
    c = bootpaths(:,4);
    ab = bootpaths(:,5);

    fh = create_figure('Histogram_Plot', 1, 5);
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

    subplot(1, 5, 1);
    shaded_hist(a);
    title(['a: ' vnames{1} '->' vnames{3}]);
    xlimit = [min(min(a), 0 - .2*nanmean(a)) max(a)];
    set(gca, 'XLim', xlimit);

    subplot(1, 5, 2);
    shaded_hist(b);
    title(['b: ' vnames{3} '->' vnames{2}]);
    xlimit = [min(min(b), 0 - .2*nanmean(b)) max(b)];
    set(gca, 'XLim', xlimit);

    subplot(1, 5, 3);
    shaded_hist(cp, xx);
    title(['c'':' vnames{1} '->' vnames{2}]);
    set(gca, 'XLim', cxlimit);

    subplot(1, 5, 4);
    shaded_hist(c, xx);
    title(['c: ' vnames{1} '->' vnames{2}]);
    set(gca, 'XLim', cxlimit);

    subplot(1, 5, 5);
    shaded_hist(ab);
    title(['ab: ' vnames{1} '->' vnames{2}]);
    xlimit = [min(min(ab), 0 - .2*nanmean(ab)) max(ab)];
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

    fignames = {'Path_Diagram' 'Slope_Plot' 'Histogram_Plot' 'HRF_Plot' 'Mediation_Scatterplots', 'Individual Effects'};

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
        dec = num2str(4);
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
