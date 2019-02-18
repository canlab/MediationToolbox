function [W, Theta, Wfull, WMinit, WJoint] = runPDM(x, y, M_tilde, Dt, methodFlag, varargin)


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% defaults
nPDM       = min(5,size(M_tilde,2));
doJointPDM = 1;

%%% parse varargin %%%
for j=1:numel(varargin)
    
    if ischar(varargin{j})
        switch lower(varargin{j})
            
            case {'npdm'}, nPDM = varargin{j+1}; varargin{j+1} = [];
                
            case {'jpdm','jointpdm'}, doJointPDM = varargin{j+1}; varargin{j+1} = [];

            otherwise, warning(['Unknown input string option: ' varargin{j}]);
        end
    end
end



%% compute principal directions of mediation

fprintf(1,'computing %d PDMs...',nPDM);

% init var
W = []; Theta = []; Wfull = []; F = [];
% pinvDt = pinv(Dt); %%% check everywhere
if strcmpi(methodFlag,'PVD')
    pinvDt = pinv(Dt);
elseif strcmpi(methodFlag,'SVD')
    pinvDt =  Dt';
end
clear Dt;

% loop PDMs
for j=1:nPDM
    
    fprintf(' %d',j);
    [w_k, theta_k, flag, wmi]= PDMN(x,y,M_tilde,W);
       
    if sign(theta_k(3))==-1
        theta_k(3:4) = -1 * theta_k(3:4);
        w_k = -1 * w_k;
    end
    
    W{j} = w_k;
    Theta{j} = [theta_k; theta_k(3)*theta_k(4)]; % add a*b to the output
    Wfull{j} = pinvDt*w_k;
    WMinit{j} = wmi;
    F(j) = flag;
    
%     save(fn,'x','y','M_tilde','Dt','F','W','Theta','Wfull','WMinit');
end

% compute joint PDM
if doJointPDM
    
    fprintf(' jointPDM');    
    WAll=[W{:}];
    tf=[Theta{:}];
    a =tf(3,:);
    b =tf(4,:);
    WJoint = pinvDt * ((a.*b)*WAll')';
    
else
    WJoint = [];
end

fprintf(1,' - done.\n');
