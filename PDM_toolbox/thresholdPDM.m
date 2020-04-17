function pdm = thresholdPDM(pdm,alpha,sigMethod)

% Threshold bootstrapped PDMs at significane level alpha while correcting
% for multiple comparisons using FDR, Bonferroni, or no correction.
%
% **Usage:**
%
%   pdm = thresholdPDM(pdm,[alpha,type])
%
% **Input:**
%
%   **pdm*
%       PDM structure with bootstrapped results returned from multivariateMediation.m
%         
%   **alpha**
%       significance threshold alpha [default: 0.05]
%         
%   **sigMethod**
%       method for multiple comparison correction [default: 'fdr']. 
%       Either 'fdr','bonf', or 'unc'
%         
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
% 4/5/2018 - Stephan Geuter
% created the file
%




% defaults
if nargin<3, sigMethod = 'fdr'; end
if nargin<2, alpha= 0.05;  end

if ~isfield(pdm,'boot'), error('no bootstrap results found in PDM struct'); end

% get PDM data
pTAll = [];

if isfield(pdm.boot,'p')
    p = pdm.boot.p;
    Wf= pdm.Wfull;
    WfSig = cell(size(Wf));
   
    
    % threshold PDMs with bootstrap results
    for k=1:numel(p)
        
        switch lower(sigMethod)
            case {'unc','uncorrected'}, pT = alpha;
            case 'fdr', p{k}(p{k}==0) = eps; pT = FDR(p{k},alpha);
            case {'bonf','bonferroni'}, pT = alpha/numel(p{k});
        end
        
        if isempty(pT), pT=0; end
        
        WfSig{k} = Wf{k} .* single(p{k}<pT);
        
        pTAll(end+1) = pT;
    end
end


% threshold joint PDM
if isfield(pdm,'WfullJoint') && isfield(pdm.boot,'pJointPDM')
    
    p = pdm.boot.pJointPDM;
    
    switch lower(sigMethod)
        case {'unc','uncorrected'}, pT = alpha;
        case 'fdr', p{1}(p{1}==0) = eps; pT = FDR(p{1},alpha);
        case {'bonf','bonferroni'}, pT= alpha/numel(p{1});
    end
    
    if isempty(pT), pT=0; end
    
    WJSig = pdm.WfullJoint .* single(p{1}<pT);
    
    pTAll(end+1) = pT;
end

% collect output
pdm.pThreshold = pTAll;
pdm.pType = sigMethod;
if exist('WfSig','var'), pdm.WfullThresh = WfSig; end
if exist('WJSig','var'), pdm.WfullJointThresh = WJSig; end



