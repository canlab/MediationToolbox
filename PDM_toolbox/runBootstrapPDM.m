function [p, Wboot, Tboot] = runBootstrapPDM(x, y, M_tilde, W , Wfull, Dt, nTrials, WMi, Bsamp, whPDM, varargin)

% function to bootstrap PDM weights. 
%
% input:
%
% x - vector of treatment variable (Nobs x 1)
% y - vector of outcomes (Nobs x 1)
% M_tilde - transformed mediator from PVD (Nobs x B)
% W - cell array with mediator weights from runPDM (1 x P cell, B x 1 in each cell)
% Wfull - cell array with voxel weights from runPDM (1 x P cell, Voxels x 1 in each cell)
% Dt - transposed weight projection matrix (B x voxels)
% nTrials - vector of number of trials per subject (Nsubj x 1)
% WMi - initial W values from the original PDM estimation (1 x P cell, B x 1 in each cell)
% Bsamp - number of bootstrap samples (default = 5000)
% whPDM - vector of indices which PDM to bootstrap OR 'JointPDM'
% Theta - cell array with original PDM path coefficients Theta 
%             (only when bootstrapping JointPDM)
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



Theta = [];


% defaults
if nargin<9, Bsamp = 5000; end
if nargin<10, whPDM = 1:numel(W); end

% get Theta for Joint PDM
if numel(varargin)==1, Theta = varargin{1}; end


% bootstrap individual PDMs
if isnumeric(whPDM)

    fprintf(['\nPDMs selected for bootstrap: ',repmat('%d ',1,numel(whPDM))],whPDM);
    
    
    % loop PDMs to bootstrap
    for k = whPDM
        
        fprintf('\nBootstrapping PDM %d with %d samples\n',k,Bsamp);
        
        if k == 1
            % first PDM
            [pv,Wb, Tb] = BootPDM(x, y, M_tilde, [], Wfull{1}, Dt, nTrials, Bsamp, WMi{1});
        else
            % later PDMs
            [pv,Wb, Tb] = BootPDM(x, y, M_tilde, W(1:k-1), Wfull{k}, Dt,nTrials, Bsamp, WMi{k});
        end
        
        p{k} = pv;
        Wboot{k} = Wb;
        Tboot{k} = Tb;
    end

    
% Joint PDM
elseif ischar(whPDM) && (strcmpi(whPDM,'jointpdm') || strcmpi(whPDM,'jdpm'))
    
    if isempty(Theta), error('Need Theta to compute Joint PDM'); end
    
    fprintf('\nPDM selected for bootstrap: JointPDM\n');

    % compute JPDM
    paths = [Theta{:}];
    paths = paths(3,:) .* paths(4,:);
    WJoint = [Wfull{:}] * paths';
    
    % bootstrap JPDM
    [p{1},Wboot{1}, Tboot{1}] = BootPDMJoint(x, y, M_tilde, W, WJoint, Dt, Bsamp, nTrials, WMi, Theta);
    
end


