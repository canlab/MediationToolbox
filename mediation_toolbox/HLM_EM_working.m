%% Create sim data

clear X Y b T V

len = 50; Nsubj = 20;    
randslopevar = 1; randintvar = 1; withinerr = 1;

fixedslope = 2; fixedint = 1;

x = zeros(len, Nsubj); x(1:2:end, :) = 1;  % fixed-effect signal (same for all subs)

c = normrnd(fixedslope,randslopevar,Nsubj,1);       % slope between-subjects variations
d = normrnd(fixedint,randintvar,Nsubj,1);  % intercept between-subjects variations

clear y

for i = 1:Nsubj
    
    y(:,i) = d(i) + c(i).*x(:,i) + normrnd(0,withinerr,len,1);

    X{i} = [ones(size(x, 1), 1) x(:, i)];

    Y{i} = y(:, i); 
end

%%
out = igls(y, x, 'verbose');  % for igls

%%
stats = glmfit_multilevel(Y, X, [], 'verbose', 'weighted');

b = stats.first_level.beta;

for i = 1:Nsubj
        
    r = Y{i} - X{i} * b(:, i);
    sigma2_hat(i)  = var(r);
    
    V{i} = inv(X{i}' * X{i}) * sigma2_hat(i);

end
%[b_star, T_star, V_new, otherstats] = EM_step(Y, b, T, V, Vgam_hat, X, W)

%

% betas, V, data (Y), X (design) 

Nsubj = length(V);

% Get dimensionality of predictor cov matrix
i = 1; Npred = 0; 
while Npred == 0, Npred = size(V{i}, 2); end
 
% Handle empty matrices (if no cov data for some bs)
for i = 1:Nsubj, if isempty(V{i}), V{i} = Inf .* eye(Npred); end, end
 
% ------------------------------------------------------------
%
% Get initial estimates
%
% ------------------------------------------------------------
 
T = nancov(b');  % between-subjects cov matrix; TOR's naive estimate (not R & B); could at least subtract expected w/i ss contrib?
 
W = eye(Npred);  % 2nd-level design matrix; intercept only = simple

% ============================================================
% Given: T, V, betas
% Estimate fixed effects (gammas), Eq. 3.31 of R & B w/o block diag
% formulation, and their variance
% ------------------------------------------------------------
for i = 1:Nsubj
    Dind{i} = T + V{i}; % rand cov + ind. inv(X'*X)*sigma^2
    
    Precis{i} = inv(Dind{i});  
end

G1 = zeros(size(W)); for i = 1:Nsubj, G1 = G1 + Precis{i}; end
G2 = zeros(Npred, 1); for i = 1:Nsubj, G2 = G2 + Precis{i} * b(:, i); end
gam_hat = inv(G1) * G2;

Vgam_hat = inv(G1);
% ============================================================

% Estimate individual slopes (betas, fixed + random)
% ------------------------------------------------------------
% Now we have 3 choices for betas: the individual b estimates, the fixed
% effect estimates, or an Empirical Bayes combination of the two
% L is "optimal" weighting matrix for each subject
% R & B eq. 3.56
I = eye(size(Precis{1}));
for i = 1:Nsubj
    L{i} = T * Precis{i};   % T(T + V{i})^-1
    b_star(:, i) = L{i} * b(:, i) + (I - L{i}) * W * gam_hat;  % weighted linear combo  of indiv and fixed fx est.
end

%%
create_figure('Gamma');

plot([0 0], gam_hat, 'r^', 'MarkerFaceColor', 'r','MarkerSize', 10);
plot([0 0], diag(T), 'b^', 'MarkerFaceColor', 'g','MarkerSize', 10);

plot([0 0], out.beta, 'rs', 'MarkerSize', 16);
plot([0 0], out.betastar, 'bs', 'MarkerSize', 16);

%gam_hat = gam_hat .* 1.5;

b_orig = b;

%gam_hat = [0 0]';

for step = 1:50

    [b, T, V, gam_hat, otherstats] = EM_step(Y, b, gam_hat, T, V, Vgam_hat, X, W, b_orig);

    for i = 1:Nsubj
        Dind{i} = T + V{i}; % rand cov + ind. inv(X'*X)*sigma^2

        Precis{i} = inv(Dind{i});
    end

    G1 = zeros(size(W)); for i = 1:Nsubj, G1 = G1 + Precis{i}; end
    Vgam_hat = inv(G1);

    create_figure('Gamma', 1, 1, 1);
    plot([step step], gam_hat, 'ro', 'MarkerFaceColor', 'r');

    plot([step step], otherstats.sigma2_btwn, 'bo', 'MarkerFaceColor', 'g');

    
end

plot_horizontal_line(fixedslope);
plot_horizontal_line(fixedint);

plot_horizontal_line(randslopevar);
plot_horizontal_line(randintvar);

