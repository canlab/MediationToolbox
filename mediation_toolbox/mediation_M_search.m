% Example: Run mediation test on clusters
% M = cat(2, allcl{1}.timeseries);
% cl2 = cl(wh);
% cluster_orthviews(cl2, {[1 1 0]})
% 
% xyz = cat(2, cl.XYZ);
% xyzmm = cat(2, cl.XYZmm);
% M = cat(2, allcl{1}.all_data);
%
% See mediation_brain
%
% Example:
%
% Do robust regression search over matrix M for mediators
% ---------------------------------------------------------------------
% names = {'rACC antic' 'Neg. emotion report' 'Neg - Neu Stimulation'}
% med_results = mediation_M_search(X, Y, M, 'names', names, 'robust')

function med_results = mediation_M_search(X, Y, M, varargin)
    med_results = mediation_search('M', X, Y, M, varargin{:});
end
