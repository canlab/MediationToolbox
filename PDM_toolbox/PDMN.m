function [w_N,theta,flag,WMi]= PDMN(x,y,m, W, varargin)

% Compute the Nth Principal Direction of Mediation
%
% This code can be used iteratively to compute each principal direction of mediation
%
% INPUT:
%
% x     - treatment (N X 1 vector)
% y     - outcome   (N X 1 vector)
% m     - mediator  (N X p matrix)
% W     - weights of previous directions (cell structure). Use [] if first direction.
%
% OPTIONAL INPUT:
%
% 'initialvalues'
% followed by a cell array of initial values used for the original PDM
% estimation (returned in WMi). Pass this argument for bootstrapping.
%
%
% OUTPUT:
%
% w_N     - weights for the Nth direction of mediation
% theta_N - parameters for Nth direction of mediation [c, c', a, b]
% lambda  - lambda value
% flag    - exit flag returned by fmincon indicating convergence
% WMi     - initial values for the minimization algorithm
%
%
% Examples:
%
% First 3 directions:
%
% [w_1, theta_1, flag] = PDMN(x,y,M_tilde,[]);
% W{1} = w_1;
% Theta{1} = theta_1;
% [w_2, theta_2, flag] = PDMN(x,y,M_tilde,W);
% W{2} = w_2;
% Theta{2} = theta_2;
% [w_3, theta_3, flag] = PDMN(x,y,M_tilde,W);
%


% ..
%     Author and copyright information:
%
%     Copyright (C) 2017 Martin Lindquist
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set-up
rng;

len = size(m,2);                % length of PDM
N1 = length(W);                 % Number of directions previously computed.
K = 50;                         % Number of initial values used
initWM = 1;                     % flag for randomization of initial weights

J = ones(size(m,1),1);          % Create design matrix for first regression
X1= [J x];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse vargin to check for initial values
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'initialvalues'
                WM  = varargin{i+1};
                K = 1;  % only use these initial values for WM
                initWM = 0; % re-set flag
        end
    end
end


%WM = ones(len,1)./sqrt(len);   % Initial value
%WM = WM./sqrt(sum(WM.^2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create orthogonality constraints

A = [];
b = [];

if (N1 >0)
    A = zeros(len,N1);          % Create matrix of weight values
    for i=1:N1,
        A(:,i) = W{i};
    end
    A  = A';
    b = zeros(N1,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find optimal values


crtmax = 0;

for i=1:K

    % initialize with different values when no values are given
    if initWM == 1
        WM = normrnd(0,1,len,1);        % Initial value for weights
        WM = WM./sqrt(sum(WM.^2));
    end
    WMtmp = WM;    
    
    options = optimset('Algorithm','interior-point','MaxFunEvals',1e6,'Maxiter',1e6,'TolX',1e-6,'TolFun',1e-6,'Display','off');
    [w_Nc, fval, flag] = fmincon(@objfun,WM,[],[],A,b,[],[],@nlconstraints,options,m,y,X1);    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate final theta values

    gamma = pinv(X1)*y;             % Compute gamma values

    mw_N = m*w_Nc;        
    alpha = pinv(X1)*mw_N;          % Compute alpha values

    X2 = [X1 mw_N];                 
    beta = pinv(X2)*y;              % Compute beta values
    
    ab = abs(alpha(2)*beta(3));
    
    if (ab > crtmax)
        crtmax = ab;
        theta = [gamma(2); beta(2); alpha(2); beta(3)];
        w_N = w_Nc;
        WMi = WMtmp;
    end
end





% finally, flip path a to be positive
if sign(theta(3))==-1
    theta(3:4) = -1 * theta(3:4);
    w_N = -1 * w_N;
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function

function mval = objfun(WM,m,y,X1)
            
    mw_N = m*WM;                
    alpha = pinv(X1)*mw_N;
    r = mw_N-X1*alpha;
    beta = r'*y/(r'*r);
%    mval = -alpha(2)*beta;
    mval = -abs(alpha(2)*beta);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrain weights to have norm 1.

 function [c,ceq] = nlconstraints(WM,m,y,X1)
 
     c = [];
     ceq = WM'*WM - 1;

 end
