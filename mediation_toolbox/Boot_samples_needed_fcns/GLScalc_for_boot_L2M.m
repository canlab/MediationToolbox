function [b, s2between_ols, stats] = GLScalc_for_boot_L2M(Y, W, X, varargin)

% function [b, s2between_ols, stats] = GLScalc_for_boot_L2M(Y, W, X, varargin)
%
% To implement bootstrapping of the 2nd level moderation, A sub-function (GLScals) 
% of "scn_stats_helper_functions.m" was separately saved as
% "GLScals_for_boot_L2M.m". 
%
% -------------------------------------------------------------------------
% GLS estimation of betas, ols variances, and full stats if 3rd output
% is requested
% -------------------------------------------------------------------------
%
% for more details, ">> help scn_stats_helper_functions". 
%
arorder = 0;
if ~isempty(varargin), arorder = varargin{1}; end

[n, k] = size(X);
nvars = size(Y, 2);

b = zeros(k, nvars);


if k == 1 && all(X == X(1)) && (nargout == 1)
    
    % intercept only; fast computation
    % Weighted mean function:
    % Very efficient when there are no predictors other than the intercept
    % Faster than looping, and faster than using mean() for equal weights
    
    b = diag(W'*Y)';
    
    
else
    % predictors; need full computation
    % OR we are getting variances, and need full computation for them
    % anyway
    % Full GLS function, needed if there are predictors
    
    
    % setup stuff
    
    Wi = cell(1, k);
    invxvx = cell(1, nvars);
    bforming = cell(1, nvars);
    equal_weights = false(1, nvars);
    
    for i = 1:nvars
        
        if all(W(:, i) == W(1, i))
            % weights are equal
            equal_weights(i) = 1;
            
        end
        
    end
    
    isweighted = 0;
    
    % get matrices for outcome vars in OLS case
    if any(equal_weights)
        
        tmp = inv(X' * X);       % Save these for later, for speed; don't need to re-use
        invxvx(equal_weights) = {tmp};
        bforming(equal_weights) = {tmp * X'};
    end
    
    % get betas (coefficients) and other necessary matrices for weighted columns
    for i = 1:nvars
        
        if ~equal_weights(i)
            isweighted = 1;
            
            Wi{i} = diag(W(:, i));              % Wi = V^-1, inverse of cov.matrix
            
            invxvx{i} = inv(X' * Wi{i} * X);       % Save these for later, for speed; don't need to re-use
            bforming{i} = invxvx{i} * X' * Wi{i};
            
        end
        
        b(:, i) = bforming{i} * Y(:, i);
        %b(:, i) = inv(X' * Wi * X) * X' * Wi * Y(:, i);
        
        
        % Add AR(p) model stuff here
        if arorder > 0
            % Warning: Tested March08...not same as fit_gls, check.****
            disp('Warning: Check AR code');
            [b(:, i),stebetatmp, varbetatmp, tmpi, tmpv, dfe_ols(i), Phi] = ...
                ar_iterate_core(X, Y(:, i), b(:, i), n, arorder);
            
        end
        
    end
end


% Optional additional computations: optional in case we want to return just b
% and go really fast (e.g., for bootstrapping)

if nargout > 1
    
    % --------------------------------------
    % * Residuals
    % --------------------------------------
    
    e = Y - X * b;         % residuals
    
    %
    % OLS residual variance
    % If bootstrapping or permuting, use this to get weights
    
    if ~(arorder > 0)
        dfe_ols = (n - k) .* ones(1, nvars);
    end
    
    % --------------------------------------
    % * Residual variance
    % --------------------------------------
    %s2between = diag(e' * e)' ./ dfe;               % var for each col of Y, 1 x nvars
    s2between_ols = (1 ./ dfe_ols .* sum(e .^ 2));  % Estimates of Sigma^2 for each col. of Y
    
end

if nargout > 2
    % Weighted variance (s2) and full stats
    
    % --------------------------------------
    % * Degrees of freedom for each column of Y
    % --------------------------------------
    dfe = dfe_ols;
    
    
    % Get design-related component of STE (includes n)
    % replace dfe with Sattherwaite approx. if necessary
    
    Xste = zeros(k, nvars);
    
    
    % --------------------------------------
    % * Residual variances, std. errors for each column of Y
    % --------------------------------------
    for i = 1:nvars
        
        Xste(:, i) = diag(invxvx{i});                  % design-related contribution to variance, k x nvars
        
        if equal_weights(i)
            % weights are equal
            
            s2between(i) = s2between_ols(i);            % Residual variance for this outcome variable
            
        else
            % all weights are not equal
            % Update dfe and s2between
            
            % R = Wi.^.5 * (eye(n) - X * invxvx * X' * Wi{i});          % Residual inducing matrix
            R = Wi{i} .^ .5 * (eye(n) - X * bforming{i});               % Residual inducing matrix
            
            Q = R * inv(Wi{i});                                         % Q = RV
            dfe(i) = (trace(Q).^2)./trace(Q * Q);                       % Satterthwaite approximation for degrees of freedom
            
            % --------------------------------------
            % * Residual variance
            % --------------------------------------
            e = R * Y(:, i);                                            % weighted residuals
            s2between(i) = diag(e' * e)' ./ dfe(i);                     % var for each col of Y, 1 x nvars
        end
        
    end
    
    % --------------------------------------
    % * Standard errors of coefficients
    % --------------------------------------
    sterrs =  ( Xste .* repmat(s2between, k, 1) ) .^ .5;
    
    % -------------------------------------------------------------------------
    % Get statistic structure from OLS regression, including p-values and conf. intervals
    % -------------------------------------------------------------------------
    
    
    stats.mean = b(1, :);           % intercept;  mean response
    stats.mean_descrip = 'Intercept of each col. of Y; (mean response if predictors are centered)';
    
    stats.beta = b;
    stats.beta_descrip = 'betas (regression coefficients), k predictors x nvars';
    
    stats.var = s2between;
    stats.var_descrip = 'Residual variance of each col. of Y';
    
    stats.ste = sterrs;
    stats.ste_descrip = 'Std. error of each beta for each col. of Y, k predictors x nvars';
    
    stats.t = b ./ sterrs;
    
    stats.dfe = dfe;
    stats.dfe_descrip = 'error DF for each col. of Y, Satterthwaite corrected if necessary;  1 x nvars';
    
    
    
    stats.e = e;
    if ~isweighted
        stats.e_descrip = 'unweighted (OLS) residuals';
    else
        stats.e_descrip = 'weighted residuals (resid. from weighted GLS model.)';
    end
    
    for i = 1:k
        
        stats.p(i, :) = min(1, (2 .* (1 - tcdf(abs(stats.t(i, :)), stats.dfe))));
        
        stats.p(i, :) = max(stats.p(i, :), eps);
        
    end
    
    stats.p_descrip = 'Two-tailed p-values';
    
end

end

