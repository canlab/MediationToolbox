% Example: Run mediation test on clusters
% M = cat(2, allcl{1}.timeseries);
% cl2 = cl(wh);
% cluster_orthviews(cl2, {[1 1 0]})
%
% xyz = cat(2, cl.XYZ);
% xyzmm = cat(2, cl.XYZmm);
% M = cat(2, allcl{1}.all_data);
%
% See mediation_brain
%
% Example:
%
% Do robust regression search over matrix M for mediators
% ---------------------------------------------------------------------
% names = {'rACC antic' 'Neg. emotion report' 'Neg - Neu Stimulation'}
% med_results = mediation_search('M', X, Y, M, 'names', names, 'robust')

% Programmers' notes
% Created by Tor Wager
% 12/4/09: Tor Wager - fixed matrix size compatibility when using multiple
% additional mediators

function med_results = mediation_search(search_var, X, Y, M, varargin)
    % defaults
    names = {'X' 'Y' 'M'};
    thresholds = .05;
    
    robust_opt = 'norobust';
    boot_opt = 'boot1';
    multilev_opt = 'hierarchical';
    arorder = 0;
    verbose = 1;
    verbstr = {'noverbose' 'verbose'};
    plotstr = {'noplots' 'plots'};
    pvals_for_boot = [3 5];
    mediation_covariates = [];
    num_additionalM = 0;
    
    % inputs
    switch(search_var)
        case 'X'
            num_regions = size(X, 2);
            wh_to_analyze = find( find_isok(X) ); %any(~(isnan(X) | X == 0)) );
            
            %%wrongN = size(Y, 2);             % number of replications at 1st level
            
        case 'Y'
            num_regions = size(Y, 2);
            wh_to_analyze = find( find_isok(Y) ); % | Y == 0)) );
             %%           N = size(X, 2);             % number of replications at 1st level

        case 'M'
            num_regions = size(M, 2);
            wh_to_analyze = find( find_isok(M) ); %| M == 0)) );
              %%          N = size(X, 2);             % number of replications at 1st level

        otherwise
            error('Unknown search_var: %s', search_var);
    end

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'names'
                    names = varargin{i+1};
                case 'thresholds'
                    thresholds = varargin{i+1};
                    
                % mediation.m analysis choice options      
                case 'robust', robust_opt = 'robust';
                case 'norobust', robust_opt = 'norobust';
                case {'boot','boot1'}, boot_opt = 'boot1';
                case 'signperm', boot_opt = 'signperm';
                case 'noboot', boot_opt = 'noboot';
                case {'multilevel','hierarchical'}, multilev_opt = 'hierarchical';
                case 'summarystats', multilev_opt = 'summarystats';

                case 'arorder', arorder = varargin{i+1};
                case {'covs','covariates', 'mediation_covariates'}, mediation_covariates = varargin{i+1};
                        
                case 'm'
                    % Number of additional mediators
                    % need num_additionalM only for matrix size compatibility
                    additionalM = varargin{i+1};
                    if iscell(additionalM)
                        num_additionalM = size(additionalM{1}, 2);
                    else
                        num_additionalM = size(additionalM, 2);
                    end
                    
                case 'noverbose', verbose = 0;
            end
        end
    end

    % for matrix size compatibility only
    npaths = 5 + 3 * num_additionalM;

    % initialize summary statistics for regions
    allpaths = zeros(num_regions, npaths);
    stes = zeros(num_regions, npaths);
    pvals = ones(num_regions, npaths);

    if verbose
        fprintf('\nTesting %d valid columns of %d total in matrix - current: %6d', length(wh_to_analyze), num_regions, 1); 
    
        disp(' ')
        disp('The mediation output below (and in the figures) is shown for'); 
        disp('the first valid voxel in the brain. ')
        disp('It is intended to give you an idea of what the ');
        disp('variables and output look like so you can check them.');
        
    end

    
    for region = wh_to_analyze
        if verbose
            if strcmp(boot_opt, 'noboot')
                % print results every 50 to save display time
                if mod(region, 50) == 0, fprintf(1, '\b\b\b\b\b\b%6d', region); end
            else
                % print every iteration
                fprintf(1, '\b\b\b\b\b\b%6d', region);
            end
        end
        
        if region == wh_to_analyze(1)
            extra_args = {plotstr{verbose+1}, verbstr{verbose+1}};
        else
            extra_args = {};
        end

        switch(search_var)
            case 'X'
                [paths, stat] = mediation(X(:,region), Y, M, varargin{:}, 'persistent', extra_args{:}, 'pvals_for_boot', pvals_for_boot, 'arorder', arorder);
                            %[paths, stat] = mediation(X(:,region), Y, M, robust_opt, boot_opt, multilev_opt, 'names', names, 'covs', mediation_covariates, 'persistent', extra_args{:});

            case 'Y'
                [paths, stat] = mediation(X, Y(:,region), M, varargin{:}, 'persistent', extra_args{:}, 'pvals_for_boot', pvals_for_boot, 'arorder', arorder); % Wani added this line 01/24/15
                %[paths, stat] = mediation(X, Y(:,region), M, robust_opt, boot_opt, multilev_opt, 'names', names, 'covs', mediation_covariates, 'persistent', extra_args{:}, pvals_for_boot, 'arorder', arorder);
            case 'M'
                [paths, stat] = mediation(X, Y, M(:,region), varargin{:}, 'persistent', extra_args{:}, 'pvals_for_boot', pvals_for_boot, 'arorder', arorder);
        end

        allpaths(region,:) = paths;
        stes(region,:) = stat.ste;
        pvals(region,:) = stat.p;
    end

    % which regions show sig. X -> M and a*b
    wh_thresh = cell(1, length(thresholds));
    
    for i=1:length(thresholds)
        wh_thresh{i} = pvals(:,1) < thresholds(i) & pvals(:,2) < thresholds(i) & pvals(:,5) < thresholds(i);
    end

    med_results = struct('search_var', search_var, 'paths', allpaths, 'ste', stes, 'pvals', pvals, 'thresholds', thresholds, 'wh_thresholds', {wh_thresh});
end

% Extra stuff for viewing output
% V = spm_vol(p);
% dim = V.dim(1:3);
% ind = sub2ind(dim, xyz(wh, 1), xyz(wh, 2), xyz(wh, 3));
% z = zeros(dim); z(ind) = 1;
% clmed = mask2clusters(z, V.mat);
%cluster_orthviews(clmed, {[1 1 0]})
% nv = cat(1, cl.numVox); clmed(nv < 5) = [];

function isok = find_isok(X)
% returns logical vector of ok (eligible) columns

tol = .000001;

isok = any( ~(isnan(X)) ) & ~all(abs(X) < tol);
end
