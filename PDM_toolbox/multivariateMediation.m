function out = multivariateMediation(varargin)

% Multivariate mediation analysis (Principial Directions of Mediation [PDM]
% Chén et al. (2017) Biostatistics; Geuter et al.(under review)).
%
% :Usage:
% ::
%       out = multivariateMediation(X,Y,M, optionalInputs); 
% 
%           ---OR---
% 
%       out = multivariateMediation(out, optionalInputs); 
%
% :Inputs:
%
%   **X**
%       a cell array (Nsubj x 1) containing a treatment vector (Ntrial x 1) 
%       for each subject
%
%   **Y**
%       a cell array (Nsubj x 1) containing an outcome vector (Ntrial x 1)
%       for each subject
%
%   **M**
%      a cell array (Nsubj x 1) containing the mediator matrix (Nvoxel X
%      Ntrial) for each subject
%
% 
%         ---OR---
%
%         
%   **out**
%       a structure with results from previous calls to
%       multivariateMedition. Minimum fields required to compute PDMs
%             out.dat.X
%             out.dat.Y
%             out.dat.M_tilde
%             out.dat.Dt
%             out.dat.B
%             out.dat.nImgs
%       These results are returned by the PVD dimension reduction step
%       e.g., out = multivariateMediation(X,Y,M,'noPDMestimation');
%       Here, the data fields X,Y,M_tilde are numeric arrays with Nobs rows
%       
% 
%  
% 
% :Optional inputs:
%
%   **'B'**
%       followed by the number of PVD components to keep. Defaults to the
%       minimum number of trials across subjects.
%
%   **noPDMestimation**
%       skip PDM estimation, e.g. just reduce the dimensionality of the
%       data with PVD.
% 
%     **nPDM**
%         followed by the number of PDMs to estimate [Default=5]. Max(nPDM)=B.
%                 
%     **jointPDM**
%          followed by [0 1] controls whether the joint PDM will
%          be computed [default = true]
%          
%     **bootPDM**
%          followed by a vector with the indices of the PDMs for which voxel 
%          weights will be bootstrapped. [default = 1:nPDM]
%          
%     **bootJointPDM**
%          followed by [0 1] controls whether the joint PDM will
%          be bootstrapped [default = false]
%                 
%     **Bsamp**
%         followed by the number of bootstrap samples for each PDM [default=5000]
%         
%     **notable**
%         don't print a table with the path coefficients
%
%     **plots**
%         plot coefficients
%
%     **save2file**
%         followed by a .mat-filename to which the results in the output  
%         structure will be saved.
%         
%     **returnbootsamples**
%         return the bootstrap samples for Wfull and Theta (large matrices)
%         
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Martin Lindquist & Stephan Geuter
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
%         
        

% Programmer's notes
%
% 02/23/2018 - Stephan Geuter
% created the file
%
% 4/3/2018 - Stephan Geuter
% added table and figure options
%
% 4/5/2018 - Stephan Geuter
% added significance thresholding
%
% 5/4/2018 - Stephan Geuter
% added SVD dimension reduction
%

%% determine input mode
%%% standard input mode - X,Y,and M as cell arrays %%%
if nargin>=3 && all([iscell(varargin{1}) iscell(varargin{2}) iscell(varargin{3})])

    % get raw mediation data
    X = varargin{1};
    Y = varargin{2};
    M = varargin{3};
    
    optInStart = 4;
    
    % init minimal output struct
    out = struct('dat',struct('X',[],'Y',[],'M_tilde',[],'Dt',[],'B',[]));

    % specific defaults
    dimReduction= 1;                        % reduce dimensionality
    if any(strcmpi(varargin,'SVD'))         % empty for SVD, default 90% var explained
        out.dat.B = [];
    else
        out.dat.B = min(cellfun(@numel,X));     % PVD dimensions 
    end
    
%%% alternative input mode - output structure with previous results to skip some steps %%%   
elseif isstruct(varargin{1})
    
    % get data structure
    out = varargin{1};
    
    optInStart    = 2;

    % input structure must at least contain fields dat.X, dat.Y,
    % dat.M_tilde, dat.Dt, dat.B
    if checkDataStruct(out)     
            dimReduction  = 0;
            fprintf('\nLow dimension data found in input-structure. Skipping dimension reduction.\n');
    else
        error('Check input structure');
    end
    
    
end



% set defaults
doPDM       = 1;                  % compute PDMs
nPDM        = min(5,out.dat.B);           % how many PDMs to compute
doJointPDM  = 1;                  % compute joint PDM
doBootPDM   = 0;                  % bootstrap individual PDMs
doBootJPDM  = 0;                  % bootstrap joint PDM
bootAllPDM  = 1;                  % bootstrap all PDMs 
Bsamp       = 5000;               % number of bootstrap samples
alpha       = 0.05;               % significance threshold for bootstrapped voxels
sigMethod   = 'fdr';              % multiple comparison correction method
doSVD       = 0;                  % use SVD instead of PVD for dimension reduction
normflag    = 'nonorm';           % z-score mediator data before PVD
printtable  = 1;                  % print a table with path coefficients
doplots     = 0;                  % plot coefficients
save2file   = 0;                  % save results to file
returnBootsamples = 0;            % return Bootstrap samples (large array)





%%% parse optional input %%%
for j=optInStart:numel(varargin)
    
    if ischar(varargin{j})
        switch lower(varargin{j}) % all keyword comparisons in lower case
                
            % dimension reduction
            case {'b'}
                if dimReduction==1
                    out.dat.B  = varargin{j+1}; varargin{j+1} = [];
                else
                    warning('No dimension reduction requested, ignoring ''B'' input');
                end
                
            case {'svd'} 
                doSVD = 1; 
                if isempty(strcmpi('b',varargin))
                    out.dat.B = [];
                end    
            
            case {'zscore'}, normflag = 'normalize';    
                     
            % PDM estimation    
            case {'nopdmestimation'}, doPDM=0;
                
            case {'npdm'}, nPDM = varargin{j+1}; varargin{j+1} = [];
                
            case {'jpdm','jointpdm'}, doJointPDM = varargin{j+1}; varargin{j+1} = [];
                
            % bootstrapping    
            case {'bootpdm'}
                doBootPDM = 1;
                if numel(varargin)>j && isnumeric(varargin{j+1})
                    whPDMBoot = varargin{j+1};
                    bootAllPDM = 0;
                end
                
            case {'bootjpdm','bootjointpdm'}, doBootJPDM = 1;
                
            case {'bsamp'}, Bsamp = varargin{j+1}; varargin{j+1} = [];
                
            case {'returnbootsamples'}, returnBootsamples = 1;
    
            case {'alpha'}, alpha = varargin{j+1}; varargin{j+1} = []; 
                
            case {'fdr'}, sigMethod = 'fdr';
                
            case {'bonf','bonferroni'}, sigMethod = 'bonf';    
            
            case {'unc','uncorrected'}, sigMethod = 'unc';    
                
            % output   
            case {'notable','notables'}, printtable=0;     
                
            case {'plot','plots'}, doplots = 1;    
                
            case {'save2file','savetofile'}, save2file = 1; outFn = varargin{j+1}; varargin{j+1} = [];
                
                
            otherwise, warning(['Unknown input string option: ' varargin{j}]);
        end
    end
end


if bootAllPDM==1
    whPDMBoot = 1:nPDM;
end



%% start computing

% project M to orthogonal, lower dimensional space
if dimReduction
    
    if doSVD
        [out.dat.X, out.dat.Y, out.dat.M_tilde, out.dat.Dt, out.dat.nImgs, out.dat.method] = mySVD(X,Y,M, out.dat.B);
        
    else
        [out.dat.X, out.dat.Y, out.dat.M_tilde, out.dat.Dt] = PVD(X,Y,M, 'B', out.dat.B,normflag);
        out.dat.nImgs  = cellfun(@numel,X(:)); out.dat.nImgs = out.dat.nImgs(:);
        out.dat.method = 'PVD';
    end
    
    % save intermediate results to file if requested
    if save2file
        save(outFn,'out'); 
    end
end



% compute PDMs
if doPDM
    
    [out.W, out.Theta, out.Wfull, out.WMinit, out.WfullJoint] = runPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.dat.Dt, out.dat.method, 'nPDM',nPDM,'jointPDM',doJointPDM);
    
    % save intermediate results to file if requested
    if save2file
        save(outFn,'out');
    end
    
    % output table
    if printtable
        printPathCoeff(out.Theta);
    end
end

% coeff figure
if doplots
    plotPathCoeff(out.Theta);
end


% bootstrap individual PDMs
if doBootPDM
    
    nImgs = out.dat.nImgs;
    if all(isfield(out,{'W','Wfull','dat','WMinit'}))
        if returnBootsamples
            [p, Wboot, Tboot] = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nImgs, out.WMinit, Bsamp, whPDMBoot);
        else
            p = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nImgs, out.WMinit, Bsamp, whPDMBoot);
        end
    else
        error('Missing PDM data bootstrapping');
    end
end



% bootstrap joint PDM
if doBootJPDM
    
   nImgs = out.dat.nImgs;    
   if all(isfield(out,{'W','Wfull','dat','WMinit','Theta'}))
       if returnBootsamples
           [pJ, WbootJ, TbootJ] = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nImgs, out.WMinit, Bsamp, 'JointPDM', out.Theta);
       else
           pJ = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nImgs, out.WMinit, Bsamp, 'JointPDM', out.Theta);
       end
   else
       error('Missing PDM data for bootstrapping');
    end
end


%% collect outputs

if doJointPDM==0
    out = rmfield(out,'WfullJoint');
end

if doBootPDM 
    out.boot.p = p;
    out.boot.Nsamples = Bsamp;
    if returnBootsamples
        out.boot.SamplesW = Wboot;
        out.boot.SamplesTheta = Tboot;
    end
end

if doBootJPDM
    out.boot.pJointPDM = pJ;
    out.boot.Nsamples = Bsamp;
     if returnBootsamples
        out.boot.SamplesWJoint = WbootJ;
        out.boot.SamplesThetaJoint = TbootJ;
    end
end

% threshold bootstrapped PDMs
if doBootPDM || doBootJPDM
    out = thresholdPDM(out,alpha,sigMethod);
end

if save2file
    save(outFn,'out','-v7.3');
end






%% inline functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function datOK = checkDataStruct(DAT)
        
        datOK = 0;
        
        % check all fieldnames
        if isfield(DAT,'dat') && all(isfield(DAT.dat,{'X','Y','M_tilde','Dt','B','nImgs'})) && ...
            isnumeric(DAT.dat.X) && isnumeric(DAT.dat.Y) && isnumeric(DAT.dat.M_tilde) ...
            && isnumeric(DAT.dat.Dt)
            
            datOK = 1;
        else
            error('Missing fields in input structure');
        end
        
        % then compare data dim's
        if all(size(DAT.dat.X)==size(DAT.dat.Y)) 
            datOK = 1;
        else
            error('Data dimension mismatch between X and Y');
        end
        if size(DAT.dat.X,1)==size(DAT.dat.M_tilde,1)
            datOK = 1;
        else
            error('Data dimension mismatch between X and M_tilde');            
        end
        if size(DAT.dat.Dt,1)==size(DAT.dat.M_tilde,2)
            datOK = 1;
        else
            error('Data dimension mismatch between Dt and M_tilde');
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [X, Y, M_tilde, Dt, nTrials, mymeth] = mySVD(treatment,outcome,mediator, nComps)
       
        fprintf('SVD for %d subjects...\n',numel(treatment));
        mymeth = 'SVD';
        
        X = []; Y = []; Mtmp = [];
        for i=1:numel(treatment)
            X    = [X; treatment{i}]; % observations x 1
            Y    = [Y; outcome{i}];   % observations x 1
            Mtmp = [Mtmp, double(mediator{i})]; % voxels x observations
            nTrials(i,1) = size(treatment{i},1);
        end
        
        [U,S,V] = svd(Mtmp,'econ');
        
        % default - explain 90% of the variance
        if isempty(nComps)
            nComps = find(cumsum(diag(S))/sum(diag(S))>=0.9,1,'first');
        elseif nComps>0 && nComps<1  
            nComps = find(cumsum(diag(S))/sum(diag(S))>=nComps,1,'first');
        end
        fprintf('keeping %d components.\n',nComps);
        
        M_tilde = V(:,1:nComps) * S(1:nComps,1:nComps);
        
        % compute back-projection matrix to fit runPDM.m:
        % Wfull{j} = V_star*w_k;
        % Dt = U(:,1:nComps);
        Dt = U(:,1:nComps)';
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function printPathCoeff(theta)
        theta = [theta{:}];
        dashes = '_____________________________________________________________________';
        fprintf('\nPDM path coefficients\n%s\n\tpath a\t\tpath b\t\tpath ab\t\tpath c''\n',dashes);
        for k=1:size(theta,2)
            fprintf('PDM%2d\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',k,theta(3,k),theta(4,k),theta(5,k),theta(2,k));
        end
        fprintf('%s\n',dashes);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plotPathCoeff(theta)
        theta = [theta{:}];
        col = lines(4);
        create_figure('PDM paths',1,3); clf;
        subplot(1,3,1);
        plot(theta(3,:),'-o','color',col(1,:),'linewidth',1.5); 
        title('path a'); xlabel('PDM #'); ylabel('coefficients');
        
        subplot(1,3,2);
        plot(theta(4,:),'-o','color',col(2,:),'linewidth',1.5);
        title('path b'); xlabel('PDM #');
        
        subplot(1,3,3);
        plot(abs(theta(5,:)),'-o','color',col(3,:),'linewidth',1.5);
        title('abs(path ab)'); xlabel('PDM #');

        ax=findobj(gcf,'Type','axes');
        set(ax,'xlim',[0.5 size(theta,2)+0.5],'FontSize',12,'Xtick',1:size(theta,2));
        drawnow;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
