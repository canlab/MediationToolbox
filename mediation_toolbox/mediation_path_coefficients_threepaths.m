function [paths, b1betas, b2betas, b3betas, cpbetas, cbetas, sterrs, intcpt, n, residm1, residm2, residy, residy2, Vxm1,Vm1m2,Vm2y,Vxy] = mediation_path_coefficients_threepaths(x,y,m1,m2,domultilev, logistic_Y)

% function [paths, b1betas, b2betas, b3betas, cpbetas, cbetas, sterrs, intcpt, n, residm1, residm2, residy, residy2, Vxm1,Vm1m2,Vm2y,Vxy] = mediation_path_coefficients_threepaths(x,y,m1,m2,domultilev)
%
% This function is a modified version of mediation_path_coefficients.m
% 
% Support function for mediation_threepaths.m
% Run mediation_threepaths.m, not this function
%
% Wani Woo, 11/04/2012
%
% added logistic option, Wani Woo, 01/13/2013

sterrs = NaN;
n = size(x,1);
intcpt = ones(n,1);

residm1 = zeros(n,1);
residm2 = zeros(n,1);
residy = zeros(n,1);
residy2 = zeros(n,1);

% THE code choices below can be simplified once robust multi-level
% option is implemented.

if domultilev
    % this is done if N > 1 and multi-level weighting is specified
    % we need to save standard errors
    [paths, sterrs, b1betas, b2betas, b3betas, cpbetas, cbetas, residm1, residm2, residy, residy2, Vxm1, Vm1m2, Vm2y, Vxy] = ols_mediation_paths_withstes(x, y, m1, m2, intcpt, n, logistic_Y);
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
function [paths,sterrs,b1betas,b2betas,b3betas,cpbetas,cbetas,residm1,residm2,residy,residy2,Vxm1,Vm1m2,Vm2y,Vxy] = ols_mediation_paths_withstes(x, y, m1, m2, intcpt, n, logistic_Y)


% Eq. 1, b1 = M1 ~ X
[b1betas, b1ste, varb1, xtxib1, pxb1, residm1, Vxm1] = regression(x, m1, intcpt, n, [], []);
% abetas should be 2 x nmediators
b1 = b1betas(1, :);  % the first coefficient for each regression (each mediator); intercept is 2nd., covs at end (ignore covs)
b1ste = b1ste(1, :); % same as above.
varb1 = varb1(1, :);

% Eq. 2,  b2 = M2 ~ X + M1
[b2betas, b2ste, varb2, xtxib2, pxb2, residm2, Vm1m2] = regression([m1 x], m2, intcpt, n, [], []);
% abetas should be 2 x nmediators
b2 = b2betas(1, :);  % the first coefficient for each regression (each mediator); intercept is 2nd., covs at end (ignore covs)
%     b2betas = b2betas_all([1 3], :);  % the first coefficient for each regression (each mediator); intercept is 2nd., covs at end (ignore covs)
b2ste = b2ste(1, :); % same as above.
varb2 = varb2(1, :);

% Eq. 3, [b3, c'] = Y ~ X, M1, M2
if ~logistic_Y
    [b3betas, b3ste_all, varb3, xtxib3, pxb3, residy, Vm2y] = regression([m2 m1 x], y, intcpt, n, [], []);
else
    nregressions = size(y, 2);
    for ii=1:nregressions
        [b3betas(:,ii), b3ste_all(:,ii), varb3(:,ii), residy(:,ii), Vm2y{ii}] = regression_logistic([m2 m1 x], y(:,ii));
    end
end

    % abetas should be 2 x nmediators
b3 = b3betas(1, :);  % the first coefficient for each regression (each mediator); intercept is 2nd., covs at end (ignore covs)
%     b3betas = b3betas_all([1 4], :);  % the first coefficient for each regression (each mediator); intercept is 2nd., covs at end (ignore covs)
b3ste = b3ste_all(1, :); % same as above.
varb3 = varb3(1, :);
cpbetas = b3betas([3 4], :);
cpste = b3ste_all(3, :);

% is X related to Y without mediator?
% Eq. 1,  Y ~ X
if ~logistic_Y
    [cbetas,  cste,  varc,  xtxi,  pxc, residy2, Vxy] = regression(x, y, intcpt, n, xtxib1, pxb1);
else
    nregressions = size(y, 2);
    for ii=1:nregressions
        [cbetas(:,ii), cste(:,ii), varc(:,ii), residy2(:,ii), Vxy{ii}] = regression_logistic(x, y(:,ii));
    end
end
        
c = cbetas(1);
cste = cste(1);


% save key outputs for path model, omitting covariates
% make each effect be in a row vector
%
%     cp = cpbetas; % x effect, controlling for all mediators
%     bste = myste(1:nmediators, :)';
%     cpste = myste(nmediators + 1, :)';
%     varb = varbeta(1:nmediators, :)';


% now this should be a row vector of ab effects, one for each mediator
if ~logistic_Y
    b1b2b3 = b1 .* b2 .* b3;
    % Multivariate delta standard error from Taylor et al.(2007)
    b1b2b3_var = b1.^2 .* b2.^2 .* varb3 + b1.^2 .* b3.^2 .* varb2 + b2.^2 .* b3.^2 .* varb1;
    b1b2b3_ste = b1b2b3_var .^ .5;

else
    % three-path version of Iacobucci(2012) eq(4), to standardize betas
    b1b2b3 = (b1/b1ste) .* (b2/b2ste) .* (b3/b3ste);
    b1b2b3_var = (b1.^2 .* b2.^2 .* varb3 + b1.^2 .* b3.^2 .* varb2 + b2.^2 .* b3.^2 .* varb1)./(varb1.*varb2.*varb3);
    b1b2b3_ste = b1b2b3_var .^ .5;
end


%     % variance of mediation, e.g., kenny, korchmaros, and bolger 2003, eq. 10
%     % single level only; for random ab effects, we need est. of cov(a,b) from
%     % group model
%     abvar = (b.^2 .* vara) + (a.^2 .* varb) + (vara.^2 .* varb.^2);
%     abste = abvar .^ .5;

% Put everything in standard output format:
% ------------------------------------------
% a(1) b(1) cp c ab(1) a(2) b(2) ab(2) a(3)...etc.
% First 5 elements are paths for first mediator, to be consistent.
% Return additional mediators after that, separately, [a b ab] for each
% mediator
paths = [b1 b2 b3 cpbetas(1,:) c b1b2b3];
sterrs = [b1ste b2ste b3ste cpste cste b1b2b3_ste];

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

function [beta_vec, stebeta, varbeta, resid, V] = regression_logistic(x,y)

    % if arorder, error('error: AR model does not work with logistic regressions yet!'); end
    
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