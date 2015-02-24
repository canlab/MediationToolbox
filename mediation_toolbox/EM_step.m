function [b_star, T_star, V_new, gam_star, otherstats] = EM_step(Y, b, gam_hat, T, V, Vgam_hat, X, W, b_orig)
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

% could pass Precis in...
for i = 1:Nsubj
    Dind{i} = T + V{i}; % rand cov + ind. inv(X'*X)*sigma^2
    Precis{i} = inv(Dind{i});  
end

for i = 1:Nsubj

    % Lambda, Eq. 3.57
    L{i} = T * Precis{i};   % T(T + V{i})^-1
    
    % beta-star, Eq. 3.56
    b_star(:, i) = L{i} * b(:, i) + (I - L{i}) * W * gam_hat;  % weighted linear combo  of indiv and fixed fx est.
    %b_star(:, i) = b(:, i);
    
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

bigX = cat(1, X{:});
pibigX = inv(bigX' * bigX);
gam_star = pibigX *  bigX' * ally;
% ------------------------------------------------------------

% % % % OR....
% % % % ***********************************************************
% % % for i = 1:Nsubj
% % %     b(:, i) = pinv(X{i}) * adjY{i};
% % % end
% % % 
% % % % This is WLS estimate, whereas naive est. would be mean(b, 2)
% % % % eq 3.31, use blk diag, simpler...
% % % % Note: if W = I, inv(blkdiag(Precis{:})) == blkdiag(Dind{:})
% % % Npred = 2;
% % % G1 = zeros(size(W)); for i = 1:Nsubj, G1 = G1 + Precis{i}; end
% % % G2 = zeros(Npred, 1); for i = 1:Nsubj, G2 = G2 + Precis{i} * b(:, i); end
% % % gam_star = inv(G1) * G2;
% % % % ------------------------------------------------------------
% % % 
% % % % Inference on fixed effects (gammas)
% % % % ------------------------------------------------------------
% % % % inv(G1) is Var(gam-hat)
% % % Vgam_hat = inv(G1);
% % % gam_t = gam_hat ./ (diag(Vgam_hat) .^ .5);
% % % % ***********************************************************

% Eq. 3.76, new T
% ------------------------------------------------------------
T_star = zeros(size(u_star, 1));
for i = 1:Nsubj, T_star = T_star + nancov(u_star') + V_star{i}; end
T_star = T_star ./ Nsubj;

% Eq. 3.77, new sigma^2
% ------------------------------------------------------------
for i = 1:Nsubj
    adjY{i} = adjY{i} - X{i} * W * gam_star;
    
    sigma2_hat(i) = (adjY{i}' * adjY{i} + trace(V_star{i})) ./ length(adjY{i});
    
    V_new{i} = inv(X{i}' * X{i}) * sigma2_hat(i);
    
end

otherstats.sigma2_hat = sigma2_hat;
otherstats.T_star = T_star;
otherstats.T_star_descrip = 'Var/Cov matrix of random effects';
otherstats.sigma2_btwn = diag(T_star);

end
