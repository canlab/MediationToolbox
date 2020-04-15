function out = multivariateMediation(varargin)

% Multivariate mediation analysis (Principial Directions of Mediation [PDM]
% Ch�n et al. (2017) Biostatistics; Geuter et al.(under review)).
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
%             out.dat.nTrials
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
%

%% determine input mode
%%% standard input mode - X,Y,and M as cell arrays %%%
if all([iscell(varargin{1}) iscell(varargin{2}) iscell(varargin{3})])

    % get raw mediation data
    X = varargin{1};
    Y = varargin{2};
    M = varargin{3};
    
    optInStart = 4;
    
    % init minimal output struct
    out = struct('dat',struct('X',[],'Y',[],'M_tilde',[],'Dt',[],'B',[]));

    % specific defaults
    dimReduction= 1;                        % reduce dimensionality
    
    dim = min(cellfun(@numel,X));
    if (dim > 1)
        out.dat.B = min(cellfun(@numel,X));     % PVD dimensions
    else
        out.dat.B = 5;                          % Default SVD dimensions
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
save2file   = 0;                  % save results to file
returnBootsamples = 0;            % return Bootstrap samples (large array)





%%% parse optional input %%%
for j=optInStart:numel(varargin)
    
    if ischar(varargin{j})
        switch lower(varargin{j}) % all keyword comparisons in lower case
                            
            case {'nopdmestimation'}, doPDM=0;
                
            case {'b'}
                if dimReduction==1
                    out.dat.B  = varargin{j+1}; varargin{j+1} = [];
                else
                    warning('No dimension reduction requested, ignoring ''B'' input');
                end
                
            case {'npdm'}, nPDM = varargin{j+1}; varargin{j+1} = [];
                
            case {'jpdm','jointpdm'}, doJointPDM = varargin{j+1}; varargin{j+1} = [];
                
            case {'bootpdm'}
                doBootPDM = 1;
                if numel(varargin)>j && isnumeric(varargin{j+1})
                    whPDMBoot = varargin{j+1};
                    bootAllPDM = 0;
                end
                
            case {'bootjdpm','bootjointpdm'}, doBootJPDM = 1;
                
            case {'bsamp'}, Bsamp = varargin{j+1}; varargin{j+1} = [];
                
            case {'save2file','savetofile'}, save2file = 1; outFn = varargin{j+1}; varargin{j+1} = [];
                
            case {'returnbootsamples'}, returnBootsamples = 1;
                
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
    
    if (size(X{1},1)>1)
        [out.dat.X, out.dat.Y, out.dat.M_tilde, out.dat.Dt] = PVD(X,Y,M, 'B', out.dat.B);
    else
        [out.dat.X, out.dat.Y, out.dat.M_tilde, out.dat.Dt] = SingleTrial_SVD(X,Y,M, 'B', out.dat.B);
    end
    
    out.dat.nTrials = cellfun(@numel,X)';
    
    % save intermediate results to file if requested
    if save2file
        save(outFn,'out'); 
    end
end



% compute PDMs
if doPDM
    
    [out.W, out.Theta, out.Wfull, out.WMinit, out.WfullJoint] = runPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.dat.Dt,'nPDM',nPDM,'jointPDM',doJointPDM);
    
    % save intermediate results to file if requested
    if save2file
        save(outFn,'out');
    end
end



% bootstrap individual PDMs
if doBootPDM
    
    nTrials = out.dat.nTrials;
    if all(isfield(out,{'W','Wfull','dat','WMinit'}))
        if returnBootsamples
            [p, Wboot, Tboot] = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nTrials, out.WMinit, Bsamp, whPDMBoot);
        else
            p = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nTrials, out.WMinit, Bsamp, whPDMBoot);
        end
    else
        error('Missing PDM data bootstrapping');
    end
end



% bootstrap joint PDM
if doBootJPDM
    
   nTrials = out.dat.nTrials;    
   if all(isfield(out,{'W','Wfull','dat','WMinit','Theta'}))
       if returnBootsamples
           [pJ, WbootJ, TbootJ] = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nTrials, out.WMinit, Bsamp, 'JointPDM', out.Theta);
       else
           pJ = runBootstrapPDM(out.dat.X, out.dat.Y, out.dat.M_tilde, out.W , out.Wfull, out.dat.Dt, nTrials, out.WMinit, Bsamp, 'JointPDM', out.Theta);
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

if save2file
    save(outFn,'out','-v7.3');
end



%% inline functions

    function datOK = checkDataStruct(DAT)
        
        datOK = 0;
        
        % check all fieldnames
        if isfield(DAT,'dat') && all(isfield(DAT.dat,{'X','Y','M_tilde','Dt','B','nTrials'})) && ...
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
end