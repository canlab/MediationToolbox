function [pv,Wboot, Tboot] = BootPDMJoint(x, y, M_tilde, W, worig, Dt, Bsamp, ntrials,WMi,torig)
% Perform bootstrap on Degrees of Mediation (DMs)
%
%
% INPUT:
%
% x      - treatment (N X 1 vector)
% y      - outcome (N X 1 vector)
% Mtilde - mediator after PVD (N X B matrix)
% Dt     - inverse weight projection matrix (B x voxels)
% W      - weights from previous directions
% worig  - weights from current direction
% Theta  - parametersfrom previous directions
% Bsamp  - number of bootstrap samples
% ntrials - number of trials per subject
% WMi     - initial values from original DM estimation
% torig   - original Theta parameters from individual PDMs
%
% OUTPUT:
%
% Wboot  - Bootstrapped weights (Voxels X Bsamp X numdir)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bsamp = 200;
% numdir = 1;
p = size(Dt,2);
len = length(y);
Wboot = zeros(p,Bsamp);
Tboot = zeros(4,Bsamp,numel(W));
nsub = length(ntrials);

stop = cumsum(ntrials);
start = [0; stop(1:(end-1))] + 1;

pinvDt = pinv(Dt);
clear Dt;
fprintf('              ');

for i=1:Bsamp,
    
    fprintf([repmat('\b',1,14) '%5d / %5d\n'],i,Bsamp);
    
    % Bootstrap samples of the data
    subind = ceil(unifrnd(0,nsub,nsub,1));
    

    ind = zeros(len,1);
    for j=1:nsub,    
        ind(start(j):stop(j)) = start(subind(j)) - 1 + ceil(unifrnd(0,ntrials(subind(j)),ntrials(j),1));
    end
               
    xB = x(ind);
    yB = y(ind);
    MB = M_tilde(ind,:);
    
    % loop 10 PDMs for each sample
    for k=1:numel(W)
        
        % Estimate the nth DM
        if k==1
            [w_n(:,k), theta_n(:,k)]= PDMN(xB,yB,MB,[],'initialvalues',WMi{k});
        else
            [w_n(:,k), theta_n(:,k)]= PDMN(xB,yB,MB,W(1:k-1),'initialvalues',WMi{k});
        end
        
%         % make path a positive
%         w_n(:,k) = w_n(:,k) * sign(theta_n(3,k));
%         theta_n(3:4,k) = theta_n(3:4,k)*sign(theta_n(3,k));
    end
    
    % compute joint PDM for one bootstrap sample
    a_n = theta_n(3,:);
    b_n = theta_n(4,:);
    Wboot(:,i) = abs(pinvDt*(w_n*(a_n.*b_n)'));
    Tboot(:,i,:) = theta_n;
end

% Compute pseudo-null

% Sort voxels by median
Q = sort(Wboot,2);
[~,ind] = sort(Q(:,ceil(Bsamp/2)));

% Choose half the voxels based on which have smallest medians 
num1 = ceil(p/4);
num2 = ceil(p/2) + num1;
QQ = Q(ind(num1:num2),:);
clear Q;

% Randomly sample voxels from this cohort and compute threshold
[~,ind2] = sort(normrnd(0,1,2*num1,1));
nsamp = 2000;  % "Sample size" for estimating null

v0 = sort(reshape(QQ(ind2(1:nsamp),:),nsamp*Bsamp,1));
clear QQ ind ind2;

% compute original joint PDM from individual PDMs
if iscell(torig)
    torig = [torig{:}];
end
worig = worig.*repmat(sign(torig(3,:)),p,1);
borig = torig(4,:).*sign(torig(3,:));
aorig = torig(3,:).*sign(torig(3,:));
worig = (worig*(aorig.*borig)');

if verLessThan('matlab','9.0')
    para = fitdist([v0; -v0],'Normal');
    pd = makedist('Normal','mu',0,'sigma',para.sigma);
    pv = (1-cdf(pd,abs(worig))) * 2;
else
    para = fitdist(v0,'HalfNormal');
    pd = makedist('HalfNormal','mu',0,'sigma',para.sigma);
    pv = 1-cdf(pd,abs(worig));
end


% phat0 = gamfit(v0);
% pv = 1 - gamcdf(worig,phat0(1), phat0(2));

% pv = zeros(p,1);
% for i=1:p, pv(i) = mean(worig(i) > v0); end; 

% alpha = 0.99;   % 1-alpha
% m = v0(ceil(length(v0)*alpha));
% H = worig > m;
% pv = mean(Wboot > 0,2);

% % Compute p-values
% t = zeros(p,1);
% pv = zeros(p,1);
% 
% opt = statset('MaxIter',1000);
% 
% for i=1:p,
%     try
%         res = fitgmdist(Wboot(i,:)',2, 'Options', opt) ;
%         t(i) = min(abs(res.mu(1)/sqrt(res.Sigma(1))) ,  abs(res.mu(2)/sqrt(res.Sigma(2))));     
%     catch
%         t(i) = 0;
%     end
%     pv(i) = 2*(1 - tcdf(t(i),len-1));
% end


