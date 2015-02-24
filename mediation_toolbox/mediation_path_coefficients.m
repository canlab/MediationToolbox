function [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2, Vxm, Vmy, Vxy] = mediation_path_coefficients(x,y,m,domultilev,dorobust,boot1,logistic_Y,varargin)
    % [paths, abetas, bbetas, cpbetas, cbetas, sterrs, intcpt, n, residm, residy, residy2, Vxm, Vmy, Vxy] = mediation_path_coefficients(x,y,m,domultilev,dorobust,boot1,arorder)
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
    % Added Var/Cov matrices as output, May 2008
    % Support for returning stats on multiple mediators, Oct 2008
    %
    % Wani Woo, Nov 2012
    % Added logistic regression (when y is a categorical variable)
    
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
    
    % var/cov matrices for each regression equation
    Vxm = {[]};
    Vmy = {[]};
    Vxy = {[]};
    
    % THE code choices below can be simplified once robust multi-level
    % option is implemented.
    
    if domultilev
        % this is done if N > 1 and multi-level weighting is specified
        % we need to save standard errors

        if dorobust
            error('Robust multi-level not implemented yet.')
        else

            [paths, sterrs, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2, Vxm, Vmy, Vxy] = ols_mediation_paths_withstes(x, y, m,intcpt,n,logistic_Y,arorder);
        end
    else
        % this is done if the summary statistics approach is specified
        % or for single-level models

        if boot1
            % bootstrapping; we will get stats later

            if dorobust
                warning('off', 'stats:statrobustfit:IterationLimit');
                [paths, abetas, bbetas, cpbetas, cbetas] = robust_mediation_paths(x, y, m, logistic_Y);
            elseif logistic_Y
                if dorobust
                    error('Robust logistic mediation results not implemented yet.')
                else
                    [paths, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = logistic_mediation_paths(x, y, m);
                end
            else
                intcpt = ones(size(x));      % faster to define intercept just once
                [paths, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths(x, y, m, intcpt);
            end

        else
            if dorobust
                error('Robust OLS non-bootstrapped results not implemented yet.')
            else
                % no bootstrapping; get stats here
                
                [paths, sterrs, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths_withstes(x, y, m,intcpt,n,logistic_Y,arorder);
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
function [paths, sterrs,abetas,bbetas,cpbetas,cbetas, residm, residy, residy2, Vxm, Vmy, Vxy] = ols_mediation_paths_withstes(x, y, m, intcpt, n, logistic_Y, varargin)

    if length(varargin) > 0, arorder = varargin{1}; else arorder = 0; end

    % NOTE: mediation_covariates are added in regression.m
    % so they are included, even though they don't appear here.
    
    % NOTE: 9/25/08 Modified to allow multiple mediators
    % each col of beta_vec, varbeta and stebeta, and each cell of Vxm, etc., refers
    % to a different regression model.
    nmediators = size(m, 2);
    
    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    [abetas, aste, vara, xtxi, px, residm, Vxm] = regression(x, m, intcpt, n, [], [], arorder);
    % abetas should be 2 x nmediators
    a = abetas(1, :);  % the first coefficient for each regression (each mediator); intercept is 2nd., covs at end (ignore covs)
    aste = aste(1, :); % same as above.
    vara = vara(1, :);

    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3,  [b,  c'] = Y ~ X + M

    % **** here we want to include all mediators, to get b's
    % TEMP: CENTER X; should not affect slopes but does affect plots...
    %[betas,  myste,  varbeta, xtxi2,  px2, residy] = regression([scale(m, 1) scale(x, 1)], y, intcpt, n, [], [], arorder);
    
    if ~logistic_Y
        [betas,  myste,  varbeta, xtxi2,  px2, residy, Vmy] = regression([m x], y, intcpt, n, [], [], arorder);
        % with one outcome (y), betas will be col. vector of slope and intercept
        % with 2+ mediators, rows are slope and intercept for each mediator
        % bbetas and cpbetas are used later to reconstruct slopes, and for
        % RB_emp bayes weighting
        % rows: m1 m2 m...x x intcpt.  cols: 1
        % tmp: ones(nmediators, 1) .*
    else
        nregressions = size(y, 2);
        for ii=1:nregressions 
            [betas(:,ii), myste(:,ii), varbeta(:,ii), residy(:,ii), Vmy{ii}] = regression_logistic([m x], y(:,ii), arorder);
        end
    end
    bbetas = [betas(1:nmediators,:); betas(nmediators + 2, :)];  % 1:nmediators are the mediators, then x (c effect), then intercept, then covs
    cpbetas = [betas(nmediators + 1,:); betas(nmediators + 2,:)];  % we want the x effect and intercept here
    
    
    % save key outputs for path model, omitting covariates
    % make each effect be in a row vector
    b = betas(1:nmediators, :)';
    cp = betas(nmediators + 1, :); % x effect, controlling for all mediators
    bste = myste(1:nmediators, :)';
    cpste = myste(nmediators + 1, :)';
    varb = varbeta(1:nmediators, :)';

    % is X related to Y without mediator?
    % Eq. 1,  Y ~ X
    if ~logistic_Y
        [cbetas,  cste,  varc,  xtxi,  px, residy2, Vxy] = regression(x, y, intcpt, n, xtxi, px, arorder);
    else
        nregressions = size(y, 2);
        for ii=1:nregressions 
            [cbetas(:,ii), cste(:,ii), varc(:,ii), residy2(:,ii), Vxy{ii}] = regression_logistic(x, y(:,ii), arorder);
        end
    end
    c = cbetas(1,:);
    cste = cste(1,:);

    % now this should be a row vector of ab effects, one for each mediator
    if ~logistic_Y
        ab = a .* b;
    else
        % Iacobucci(2012) eq(4)
        ab = a./aste .* b./bste;
    end

    % variance of mediation, e.g., kenny, korchmaros, and bolger 2003, eq. 10
    % single level only; for random ab effects, we need est. of cov(a,b) from
    % group model
    if ~logistic_Y
        abvar = (b.^2 .* vara) + (a.^2 .* varb) + (vara .* varb);
        abste = abvar .^ .5;
    else
        % Iacobucci(2012) eq(5)
        abvar = (a./aste).^2 + (b./bste).^2 + 1;
        abste = abvar .^ .5;
    end
    
    % Put everything in standard output format:
    % ------------------------------------------
    % a(1) b(1) cp c ab(1) a(2) b(2) ab(2) a(3)...etc.
    % First 5 elements are paths for first mediator, to be consistent.  
    % Return additional mediators after that, separately, [a b ab] for each
    % mediator
    paths = [a(1) b(1) cp c ab(1)];
    sterrs = [aste(1) bste(1) cpste cste abste(1)];
    
    % get outut for additional mediators
    % ****NOTE: returning these will require additional changes to
    % mediation to handle output.  For now, just ignore them and report the
    % first mediator.
    for i = 2:nmediators
       paths = [paths a(i) b(i) ab(i)];
       sterrs = [sterrs aste(i) bste(i) abste(i)];
    end
    
end


function [beta_vec, stebeta, varbeta, xtxi, px, resid, V] = regression(x,y,intcpt,n,xtxi,px,varargin)
    % [beta_vec,stebeta,xtxi,px, resid] = regression(x,y,intcpt,xtxi,px,n)

    % NOTE: 9/25/08 Modified to allow y to have multiple columns
    % (e.g., for multiple mediators)
    % each col of beta_vec, varbeta and stebeta, and each cell of V, refers
    % to a different regression model.

    global mediation_covariates

    if ~isempty(varargin), arorder = varargin{1}; else arorder = 0; end
    X = [x intcpt mediation_covariates];

    p = size(X,2);

    if nargin < 5 || isempty(xtxi) || isempty(px)
        xtxi = inv(X'*X);
        px = xtxi * X';          % pinv(X), X' beta-forming matrix
    end

    % Step 1: Find the OLS solution for each y col.  : each col of beta_vec
    % is a separate regression.
    nregressions = size(y, 2);

    beta_vec = px * y;                                % Beta values
    resid = y - X * beta_vec;                         % Residuals
    sigmasq = 1/(n-p) * sum(resid.^2);            % Estimate of Sigma

    V = cell(1, nregressions);
    for i = 1:nregressions
        V{i} = xtxi .* sigmasq(i);                % Var/Cov mtx, Precision^-1, used in weighted est. and empirical bayes
    end

    % Add AR(p) model stuff here
    % Many things we might want to save requested as outputs, but nothing
    % is done with them yet.
    if arorder > 0
        if nregressions > 1
            disp('WARNING: AR model not yet tested with multiple regressions/multiple mediators.  This will probably break!');
        end
        [beta_vec, stebeta, varbeta] = ar_iterate_core(X,y,beta_vec,n,arorder);
        %[beta_vec,stebeta, varbeta, i, v, df, Phi] = ar_iterate_core(X,y,beta_vec,n,arorder);
    else
        for i = 1:nregressions
            %varbeta = sigmasq.*xtxi;                    % already got this
            varbeta(:, i) = diag(V{i});                    % var of coefficients
        end
        stebeta = varbeta.^.5;
    end
    
end


function [beta_vec, stebeta, varbeta, resid, V] = regression_logistic(x,y, arorder)

    if arorder, error('error: AR model does not work with logistic regressions yet!'); end
    
    global mediation_covariates

    nmediators = size(x,2);
    
    X = [x mediation_covariates];
    
    k = unique(y);
    if length(k) == 2 && k(1) == 0
        % modified this part to fix a warning message related to "iteration limit reached."
        y1 = y; y1(y1==0)=2;  
        [dummy,dummy,stat] = mnrfit(X,y1);
        stat.resid = stat.resid(:,1);
        % [dummy,dummy,stat] = glmfit(mx,y,'binomial', 'link', 'logit'); %
    else 
        [dummy,dummy,stat] = mnrfit(X,y);
    end
    
    beta_vec = stat.beta([2:(2+nmediators-1) 1],:);
    stebeta = stat.se([2:(2+nmediators-1) 1],:);
    varbeta = diag(stat.covb);
    varbeta = varbeta([2:(2+nmediators-1) 1],:);
    resid = stat.resid;
    if size(mediation_covariates) == 0
        V = stat.covb([2:(2+nmediators-1) 1],[2:(2+nmediators-1) 1]);
    else
        V = stat.covb([2:(2+nmediators-1) 1 2+nmediators:end],[2:(2+nmediators-1) 1 2+nmediators:end]);
    end    
end

function [paths, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = logistic_mediation_paths(x, y, m)
    
    global mediation_covariates
    
    nmediators = size(m, 2);
    
    xx = [x mediation_covariates];
    
    px = inv(xx'*xx)*xx'; % X beta-forming matrix
    mx = [m xx];
    pmx = inv(mx'*mx)*mx'; % M+X beta-forming matrix
    
    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    % do for each mediator
    
    for i = 1:nmediators
        [dummy,dummy,stat] = glmfit(xx,m(:,i));
        all_abetas(:,i) = stat.beta;
        stebeta(:,i) = stat.se;
        residm(:,i) = stat.resid;
    end
    
    abetas = all_abetas([2 1], :);        % slope and intercept, for plotting
    a = all_abetas(2, :);               % a path coefficients, related to x
    aste = stebeta(2, :);
    
    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3, [b, c'] = Y ~ X + M (Here, Y = categorical variables)
    k = unique(y);
    if length(k) == 2 && k(1) == 0
        % modified this part to fix a warning message related to "iteration limit reached."
        y1 = y; y1(y1==0)=2;  
        [dummy,dummy,stat] = mnrfit(mx,y1);
        stat.resid = stat.resid(:,1);
        % [dummy,dummy,stat] = glmfit(mx,y,'binomial', 'link', 'logit');
    else 
        [dummy,dummy,stat] = mnrfit(mx,y);
    end
    
    bbetas(1, :) = stat.beta(2:(2+nmediators-1))'; % slopes for each b effect
    bbetas(2, :) = stat.beta(1,:);                  % intercept for each b effect
    cpbetas = stat.beta([nmediators + 2 1],:);  % slope and intercept for x direct effect
    b(1, :) = bbetas(1, :);
    bste = stat.se(2:(2+nmediators-1))';
        
    cp = cpbetas(1);
    residy = stat.resid;
    
    % is X related to Y without mediator?
    % Eq. 1, Y ~ X
    if length(k) == 2 && k(1) == 0
        % modified this part to fix a warning message related to "iteration limit reached."
        y1 = y; y1(y1==0)=2;  
        [dummy,dummy,stat] = mnrfit(xx,y1);
        stat.resid = stat.resid(:,1);
    else 
        [dummy,dummy,stat] = mnrfit(xx,y);
    end
    
    cbetas = stat.beta;
    c = cbetas(2);              % c total path coefficients
    residy2 = stat.resid;
    cbetas = cbetas([2 1]);       % slope and intercept for c
    
    ab = a./aste .* b./bste;            % vector of ab for each mediator
    paths = [a(1) b(1) cp c ab(1)];
    
    for i = 2:nmediators
        paths = [paths a(i) b(i) ab(i)];
    end
    

end

% -------------------------------------------------------------------------
% Compute path coefficients for OLS regression
% This version has more complete output, to save more info, but is somewhat
% slower than the other version below, which is used for bootstrapping
% -------------------------------------------------------------------------
function [paths, abetas, bbetas, cpbetas, cbetas, residm, residy, residy2] = ols_mediation_paths(x, y, m, intcpt)

    %%% *****action item: Need to update for multiple mediator output
    
    % ***  m can have multiple columns
    
    global mediation_covariates
    
    nmediators = size(m, 2);
    
    xx = [x intcpt mediation_covariates];
    %     px = pinv(xx);          % X beta-forming matrix    inv(X'*X)*X'
    %     pmx = pinv([m xx]);     % M+X beta-forming matrix
    
    px = inv(xx'*xx)*xx'; % X beta-forming matrix
    mx = [m xx];
    pmx = inv(mx'*mx)*mx'; % M+X beta-forming matrix
    
    % Is X related to the mediator?
    % Eq. 2, a = M ~ X
    % do for each mediator
    
    all_abetas = px * m;                % [x intcpt no. covs] x no. mediators
    abetas = all_abetas(1:2, :);        % slope and intercept, for plotting
    a = all_abetas(1, :);               % a path coefficients, related to x
    residm = m - xx * all_abetas;
    
    % Are M (b path) and X (c' path) independently related to Y?
    % Eq. 3, [b, c'] = Y ~ X + M
    all_bbetas = pmx * y;           % [no. mediators x intcpt no. covs] x 1
    
    bbetas(1, :) = all_bbetas(1:nmediators)';       % slopes for each b effect
    bbetas(2,:) = all_bbetas(nmediators + 2, :);    % intercept for each b effect
    
    cpbetas = all_bbetas([nmediators + 1 nmediators + 2]);  % slope and intercept for x direct effect
    b(1, :) = bbetas(1, :) ;
    cp = cpbetas(1);
    
    residy = y - [m xx] * all_bbetas;
    
    % is X related to Y without mediator?
    % Eq. 1, Y ~ X
    cbetas = px * y;
    
    c = cbetas(1);              % c total path coefficients
    
    residy2 = y - xx * cbetas;
    cbetas = cbetas(1:2);       % slope and intercept for c
    
    ab = a .* b;            % vector of ab for each mediator
    paths = [a(1) b(1) cp c ab(1)];
    
    for i = 2:nmediators
        paths = [paths a(i) b(i) ab(i)];
    end

end

% -------------------------------------------------------------------------
% Compute path coefficients for robust IRLS regression
% This version has more complete output, to save more info, but is somewhat
% slower than the other version below, which is used for bootstrapping
% -------------------------------------------------------------------------
function [paths, abetas,bbetas, cpbetas, cbetas] = robust_mediation_paths(x, y, m, logistic_Y)
    
    if logistic_Y, error('Robust logistic mediation results not implemented yet.'); end

    global mediation_covariates

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

