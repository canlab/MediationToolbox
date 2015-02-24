function [beta,stebeta, varbeta, i, v, df, Phi] = ar_iterate_core(X,y,beta,n,p)
% [beta,stebeta, varbeta, i, v, df, Phi] = ar_iterate_core(X,y,beta,n,p)
%
% Martin Lindquist & Tor Wager
% Tested same as fit_gls on feb 6

% Stuff needed for future iterations
iV = eye(n);
A = eye(n);
Phi = 0;
betaold = 0;
 
% Steps 2-4: Find the GLS solution by iteration 
% If p=0, skip this step. Appropriate solution already calculated above.
% Continue iteration until convergence or for at most 10 loops.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note to Jack and Tor: Keith Worsley uses a similar algorithm when fitting
% a GLM with an AR(p) noise model. However, he skips the iterative step
% and only goes through the loop one time. He claims that this is enough. I
% am not entirely convinced, therefore it is probably better to go through
% a few times if needed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
maxIter = 10;

% tolerance: 1% of min beta
tol = min(abs(beta) * .01);

i=1;                % Set counter
 
while (i < maxIter & sum((beta-betaold).^2) > tol & p > 0)
    
    resid = y - X * beta;             % Calculate residuals of current model fit
    
    % Estimate AR parameters using residuals
% % %     [a,e] = aryule(resid,p);  
% % %     
% % %     Phi = -a(2:end)';
% % %     %Same as above: 
% % % %     Phi =zeros(length(a)-1,1);      
% % % %     Phi(1:p) = -a(2:(p+1));
% % %     sigma = sqrt(e);
% % %  
% % %     % Find the inverse of the covariance matrix
% % %     A = diag(ones(n,1));
% % %  
% % %     for j=1:p,
% % %        A = A + diag(-Phi(j)*ones((n-j),1),-j);
% % %     end
% % %  
% % %     iV = A*A';  % The inverse of the covariance matrix
% % %  

% This code produces the same result: Check for accuracy ****
    [Phi,v] = aryule(resid,p);   
    Phi = Phi(2:end);
 
    % Find the inverse of the covariance matrix (iV)
    % This code essentially puts the AR params on the diagonals of iV
    A = diag(ones(n,1));
 
    for j = 1:p
       A = A + diag(Phi(j)*ones((n-j),1),-j);
    end
 
    iV = A*A'; 
    
    % intermediate calculations we can use later if converged
    xtxi = inv(X' * iV * X);
    betaform = xtxi * X' * iV;
    
    betaold = beta;                 % Set old solution to be betaold
    beta = betaform * y;    % Calculate new solution ***CHECK
 
    i = i+1;                        % Add one to counter
    
end
 
% now make Phi neg, that we're done iterating
Phi = -Phi';
 
V = inv(iV);                                % Covariance matrix
R = (eye(n) - X * betaform);                  % Residual inducing matrix
 
W = R*A*V*A';
df = (trace(W).^2)./trace(W*W);       % Satterthwaite approximation for degrees of freedom
 
varbeta = diag(v .* xtxi);                  % Var(beta)
stebeta = varbeta.^.5;           % tor added as output
 
end
