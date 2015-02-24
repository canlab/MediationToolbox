function [b_star, T_star, V_new, otherstats] = EM_step(Y, b, T, V, Vgam_hat, X, W)
% Setup things we'll need

otherstats = [];

Nsubj = length(V);

Tinv = inv(T);

% Handle empty matrices (if no cov data for some bs)
for i = 1:Nsubj, if isempty(V{i}), V{i} = Inf .* eye(Npred); end, end
  
I = eye(size(T));


% ============================================================
% Given: T, V, gamma, var(gamma)
% Find: beta*, u*, V*

% Get precision-weighted shrinkage matrix (Lambda, L)
% Get Empirical Bayes' weighted betas for each random level (i.e., subject)
% Get V*, updated individual (first-level) covariance matrix

for i = 1:Nsubj

    % Lambda, Eq. 3.57
    L{i} = T * Precis{i};   % T(T + V{i})^-1
    
    % beta-star, Eq. 3.56
    b_star(:, i) = L{i} * b(:, i) + (I - L{i}) * W * gam_hat;  % weighted linear combo  of indiv and fixed fx est.
    
    % residuals of random-effects levels
    % also est. of mean of random effect level i
    u_star(:, i)  = b_star(:, i) - gam_hat;
    
    % Variance of u, conditional on other params and data
    % Var of random effect level i
    % This is totally different from V
    Vinv{i} = inv(V{i});
    I_minus_Li = I - L{i};
    
    V_star{i} = inv(Vinv{i} + Tinv) + I_minus_Li * Vgam_hat * I_minus_Li';
end
% ============================================================


% ============================================================
% E - step and M-step together: calc E(sumsq) and update
% Given: u*, v*, 
% Get new gamma, T, sigma^2
% ------------------------------------------------------------
% Get new Gamma : Eq 3.74

for i = 1:Nsubj
    adjY{i} = Y{i} - X{i} * u_star(:, i);
end

ally = cat(1, adjY{:});

bigX = blkdiag(X{:});
pibigX = inv(bigX' * bigX);
gam_star = pibigX * ally;
% ------------------------------------------------------------


% Eq. 3.76, new T
% ------------------------------------------------------------
T_star = zeros(size(u_star, 1));
for i = 1:Nsubj, T_star = T_star + nancov(u_star') + V_star{i}; end
T_star = T_star ./ Nsubj;

% Eq. 3.77, new sigma^2
% ------------------------------------------------------------
for i = 1:Nsubj
    adjY{i} = adjY{i} - W * X{i} * gam_star;
    
    sigma2_hat(i) = (adjY{i}' * adjY{i} + V_star{i}) ./ length(adjY{i});
    
    V_new{i} = inv(X{i}' * X{i}) * sigma2_hat(i);
    
end

otherstats.sigma2_hat = sigma2_hat;

end
