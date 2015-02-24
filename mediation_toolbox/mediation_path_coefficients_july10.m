function [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2] = mediation_path_coefficients(x,y,m,domultilev,dorobust,boot1,varargin)
    % [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2] = mediation_path_coefficients(x,y,m,domultilev,dorobust,boot1,arorder)
    %
    % Returns single-level mediation path coefficients and standard errors (in some cases) using
    % subfunctions optimized for computational efficiency, based on the
    % configuration of input options.
    %
    % --------------------------------------------
    % Compute path coefficients (and standard errors in some cases)
    % If bootstrapping is 'on', skips standard errors
    % If bootstrapping is 'off', returns std. errors, intercept, and n
    % needed for OLS stats
    % --------------------------------------------
    %
    % Support function for mediation.m
    % Run mediation.m, not this function
    %
    % Tor Wager, Dec 2006

    % --------------------------------------------
    % Compute path coefficients
    % --------------------------------------------

    % choices are: hierarchical weighting (multilev), robust IRLS, and
    % bootstrapping
    % We want to do something appropriate and computationally efficient for
    % each case

    if length(varargin) > 0, arorder = varargin{1}; else arorder = 0; end
    
    if isempty(arorder), arorder = 0; end
    
    sterrs = NaN;
    n = size(x,1);
    intcpt = ones(n,1);
    
    residm = zeros(n,1);
    residy = zeros(n,1);
    residy2 = zeros(n,1);
                
    % THE code choices below can be simplified once robust multi-level
    % option is implemented.
    
    if domultilev
        % this is done if N > 1 and multi-level weighting is specified
        % we need to save standard errors

        if dorobust
            error('Robust multi-level not implemented yet.')
        else

            [paths, sterrs, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths_withstes(x, y, m,intcpt,n,arorder);
        end
    else
        % this is done if the summary statistics approach is specified
        % or for single-level models

        if boot1
            % bootstrapping; we will get stats later

            if dorobust
                warning('off', 'stats:statrobustfit:IterationLimit');
                [paths, abetas, bbetas, cpbetas, cbetas] = robust_mediation_paths(x, y, m);
            else
                intcpt = ones(size(x));      % faster to define intercept just once
                [paths, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths(x, y, m, intcpt);
            end

        else
            if dorobust
                error('Robust OLS non-bootstrapped results not implemented yet.')
            else
                % no bootstrapping; get stats here
                
                [paths, sterrs, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths_withstes(x, y, m,intcpt,n,arorder);
            end
        end

    end

end




%__________________________________________________________________________
%
%
% *** Path Computation Support ***
%
%__________________________________________________________________________

% -------------------------------------------------------------------------
% Compute path coefficients for OLS regression with St. Errs
% This version has more complete output, including STEs, for multilevel
% -------------------------------------------------------------------------
function [paths, sterrs,abetas,bbetas,cpbetas,cbetas, residm, residy, residy2] = ols_mediation_paths_withstes(x, y, m, intcpt, n, varargin)

    if length(varargin) > 0, arorder = varargin{1}; else arorder = 0; end
    
    % mediation_coefficients are handled in each regression
    % using a call to 'regression' subfunction
    
    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    [abetas, aste, vara, xtxi, px, residm] = regression(x, m, intcpt, n, [], [], arorder);
    a = abetas(1);
    aste = aste(1);
    vara = vara(1);

    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3,  [b,  c'] = Y ~ X + M

    [betas,  myste,  varbeta, xtxi2,  px2, residy] = regression([m x], y, intcpt, n, [], [], arorder);
    bbetas = [betas(1) betas(3)]';
    cpbetas = [betas(2) betas(3)]';
    b = betas(1);
    cp = betas(2);
    bste = myste(1);
    cpste = myste(2);
    varb = varbeta(1);

    % is X related to Y without mediator?
    % Eq. 1,  Y ~ X
    [cbetas,  cste,  varc,  xtxi,  px, residy2] = regression(x, y, intcpt, n, xtxi, px, arorder);
    c = cbetas(1);
    cste = cste(1);

    ab = a .* b;

    % variance of mediation, e.g., kenny, korchmaros, and bolger 2003, eq. 10
    % single level only; for random ab effects, we need est. of cov(a,b) from
    % group model
    % Should check this *****
    abvar = (b^2 * vara) + (a^2 * varb) + (vara^2 * varb^2);
    abste = sqrt(abvar);

    paths = [a b cp c ab];
    sterrs = [aste bste cpste cste abste];
    
    abetas = abetas(1:2);  % truncate to make consistent with pre-covariate version
    cbetas = cbetas(1:2);  % truncate to make consistent with pre-covariate version

end



function [beta, stebeta, varbeta, xtxi, px, resid] = regression(x,y,intcpt,n,xtxi,px,varargin)
    % [beta,stebeta,xtxi,px, resid] = regression(x,y,intcpt,xtxi,px,n)

    global mediation_covariates
    
    if length(varargin) > 0, arorder = varargin{1}; else arorder = 0; end
    X = [x intcpt mediation_covariates];

    p = size(X,2);
    
    if nargin < 5 || isempty(xtxi) || isempty(px)
        xtxi = inv(X'*X);
        px = xtxi * X';          % pinv(X), X' beta-forming matrix
    end

    % Step 1: Find the OLS solution
    beta = px * y;                                % Beta values
    resid = y - X * beta;                         % Residuals
    sigmasq = 1/(n-p) * sum(resid.^2);            % Estimate of Sigma  

    % Add AR(p) model stuff here
    % Many things we might want to save requested as outputs, but nothing
    % is done with them yet.
    if arorder > 0
        [beta,stebeta, varbeta] = ar_iterate_core(X,y,beta,n,arorder);
        %[beta,stebeta, varbeta, i, v, df, Phi] = ar_iterate_core(X,y,beta,n,arorder);
    else
        varbeta = sigmasq.*xtxi;                    % Var(beta)
        varbeta = diag(varbeta);                    % var of coefficients
        stebeta = varbeta.^.5;
    end
end


% -------------------------------------------------------------------------
% Compute path coefficients for OLS regression
% This version has more complete output, to save more info, but is somewhat
% slower than the other version below, which is used for bootstrapping
% -------------------------------------------------------------------------
function [paths, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths(x, y, m,intcpt)
    i = 1;

    %%% *****action item: replace pinv with expression below
    %%% *****action item: could consolidate some code in these 2 functions
    
    global mediation_covariates
    
    xx = [x intcpt mediation_covariates];
%     px = pinv(xx);          % X beta-forming matrix    inv(X'*X)*X'
%     pmx = pinv([m xx]);     % M+X beta-forming matrix

    px = inv(xx'*xx)*xx'; % X beta-forming matrix
    mx = [m xx];
    pmx = inv(mx'*mx)*mx'; % M+X beta-forming matrix
    
    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    tmp = px * m;
    abetas(:,i) = tmp;
    a(i, 1) = tmp(1);
    residm = m - xx * abetas;

    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3, [b, c'] = Y ~ X + M
    tmp = pmx * y;
    bbetas(:,i) = tmp([1 3]');
    cpbetas(:,i) = tmp([2 3]');
    b(i, 1) = tmp(1);
    cp(i, 1) = tmp(2);

    residy = y - [m xx] * tmp;
    
    % is X related to Y without mediator?
    % Eq. 1, Y ~ X
    tmp = px * y;
    cbetas(:,i) = tmp;
    c(i, 1) = tmp(1);

    residy2 = y - xx * cbetas;
    
    ab = a .* b;
    paths = [a b cp c ab];
    
    abetas = abetas(1:2);  % truncate to make consistent with pre-covariate version
    cbetas = cbetas(1:2);  % truncate to make consistent with pre-covariate version
end

% -------------------------------------------------------------------------
% Compute path coefficients for robust IRLS regression
% This version has more complete output, to save more info, but is somewhat
% slower than the other version below, which is used for bootstrapping
% -------------------------------------------------------------------------
function [paths, abetas,bbetas, cpbetas, cbetas] = robust_mediation_paths(x, y, m)
    i = 1;

    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    abetas(:,i) = robustfit([x mediation_covariates], m);
    abetas(1:2,i) = abetas([2 1], i);  % put intercept last
    a(i, 1) = abetas(1, i);

    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3, [b, c'] = Y ~ X + M
    tmp(:,i) = robustfit([m x mediation_covariates], y);
    tmp(1:3,i) = tmp([2 3 1], i);  % put intercept last
    if nargout > 1
        bbetas(:,i) = tmp([1 3], i);
        cpbetas(:,i) = tmp([2 3], i);
    end
    b(i, 1) = tmp(1, i);
    cp(i, 1) = tmp(2, i);

    % is X related to Y without mediator?
    % Eq. 1, Y ~ X
    cbetas(:,i) = robustfit([x mediation_covariates], y);
    cbetas(1:2,i) = cbetas([2 1], i);  % put intercept last
    c(i, 1) = cbetas(1, i);

    ab = a .* b;
    paths = [a b cp c ab];
    
    
    abetas = abetas(1:2);  % truncate to make consistent with pre-covariate version
    cbetas = cbetas(1:2);  % truncate to make consistent with pre-covariate version
end

