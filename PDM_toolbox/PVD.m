function [X, Y, M_tilde, Dt, B, nTrials] = PVD(varargin)

% Perform Population Value Decomposition (PVD) prior to estimating Principal 
% Directions of Mediation (PDM)
%
% :Usage:
% ::
%       [x, y, M_tilde, Dt, B, nTrials] = PVD(files, B, varargin); 
% 
%           ---OR---
% 
%       [x, y, M_tilde, Dt B, nTrials] = PVD(x,y,M, varargin); 
%
% :Inputs:
%
%   **files**
%       a cell array with one cell per subject containing the path to the
%       subject's data file. Each subject's data file contains three 
%       variables named:
%       'treatment' (Ntrial X 1 vector)     (lower case)
%       'outcome' (Ntrial X 1 vector)       (lower case)
%       'mediator' (Nvoxel X Ntrial matrix) (lower case)
%
%   **B**
%       The number of components to compute with B >= min(n_trials).
% 
%         ---OR---
%         
%   **x**
%       a cell array containing a treatment vector (Ntrial x 1) for each
%       subject
%
%   **y**
%       a cell array containing an outcom vector (Ntrial x 1) for each
%       subject
%
%   **m**
%      a cell array with on cell per subject containing the mediator 
%      (Nvoxel X Ntrial matrix)
%
%
% 
% :Optional inputs:
%
%   **'B'**
%       followed by the number of components to compute. Defaults to the
%       minimum number of trials across subjects.
%       B must be specified as 2nd input argument if using the filename
%       input mode.
%       This option will override the initial B-value when providing
%       filenames for each subject (2nd input argument).
%
%   **'normalize'**
%       the mediator matrix will be z-scored for each subject prior to the 
%       PVD  computation.
%
%   **'nonormalize'**
%       [default] The mediator matrix will be used as is for the PVD  computation.
%
%
% :Outputs:
%
%   **x**
%       the treatment concatenated across subjects 
%
%   **y**
%       the outcome concatenated across subjects 
%
%   **Mtilde**
%       mediator after PVD (N X B matrix)
%
%   **Dt**
%       Weight projection matrix (B x voxels)
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


%% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:
%
% N     - Number of subjects
% B     - Number of components (should be lesss than the min number of
%               trials per subject)
% files - cell array of filepaths to each subjects' fMRI data (N x 1
% vector)
%
% OUTPUT:
%
% x      - treatment (N X 1 vector)
% y      - outcome (N X 1 vector)
% Mtilde - mediator after PVD (N X B matrix)
% D      - Weight projection matrix (B x voxels)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% defaults %%%
dozscore  = 0;
B = [];


%%% process input %%%
% check for filenames
if sum(cellfun(@ischar,varargin{1})) == numel(varargin{1})
    % input mode is filenames
    readfiles = 1;
    inpStart = 3;
    files = varargin{1};
    B     = varargin{2};
    
elseif iscell(varargin{1}) && iscell(varargin{2}) && iscell(varargin{3})
    % input mode is data as cell arrays
    readfiles = 0;
    inpStart = 4;
    treatment = varargin{1};
    outcome  = varargin{2};
    mediator = varargin{3};
end


%%% parse true varargin %%%
for j=inpStart:numel(varargin)
    
    if ischar(varargin{j})
        switch lower(varargin{j})
            
            case {'b'}, B = varargin{j+1}; varargin{j+1} = [];
                
            case {'nonormalize','nonorm'}, dozscore = 0;
                
            case {'normalize','norm'}, dozscore = 1;
                
            otherwise, warning(['Unknown input string option: ' varargin{j}]);
        end
    end
end



% get subject number
N = numel(varargin{1});

if isempty(B) && ~readfiles
    % direct data input and B unspecified, get minimum number of trials
    B =  min(cellfun(@numel,treatment));
end

fprintf('PVD for %d subjects with %d components\n',N,B);


% init variables.
X = [];
Y = [];
M = [];
V = [];
Res = []; 


for i=1:N

    % Load data - should be mat-file containing M (voxels-by-trial), x
    % (trial x 1), and y (trial x 1)
    if readfiles
        load(files{i},'treatment','outcome','mediator');
        
        X = [X; treatment];  % Trial-specific treatment for subjects
        Y = [Y; outcome];  % Trial-specific outcome for subjects
        M = double(mediator);
        nTrials(i,1) = size(treatment,1);
        
    else
        % data directly passed to functions
        X = [X; treatment{i}];
        Y = [Y; outcome{i}];
        M = double(mediator{i});
        nTrials(i,1) = size(treatment{i},1);
    end
    
    if B>size(M,2), error('Each subject needs at least B trials'); end
    
    % normalize mediator
    if dozscore
        nz    = M~=0; % zeros are usually masked out in fmri_data
        M(nz) = (M(nz)-nanmean(M(nz))) ./ nanstd(M(nz));
    end
    
    % Compute and store trial-specific SVDs
    [Vi,Si,Ui] =svd(M,0);
    Res{i}.U = Ui;
    Res{i}.S = Si;
    Res{i}.V = Vi;
    disp(i);
   
%    V = [V Vi(:,1:B)];

    V = [V Vi];
end


%%    

% Compute D matrix

VV = V'*V;

% Perfom SVD on VV 

[Ct, Bt, ~] = svd(VV);
Bt_inv= diag(1./sqrt(diag(Bt)));
CB_inv = Ct*Bt_inv;

% Obtain A_k tilde, i.e. A_k tilde = V_k tilde C_tilde B_tilde ^ {-1}

A_tilde = V*CB_inv;

D = A_tilde(:,1:B);
Dt= D'; % return D-transpose, needed for reconstruction of voxelweights later
% Dt = pinv(D'); % return pinv(D-transpose), no further transformations needed later

 
%%

% Compute M_tilde

M_tilde = [];

for i = 1:N,
    Ui = Res{i}.U;
    Si = Res{i}.S;
    Vi = Res{i}.V;
    
 %   tmp = Ui(:,1:B)*Si(1:B,1:B)*Vi(:,1:B)'*D;
    tmp = Ui*Si(:,1:B)*Vi(:,1:B)'*D;
    M_tilde = [M_tilde; tmp];
end




