% See mediation_brain
%
% Example:
%
% Do robust regression search over matrix Y for mediators
% ---------------------------------------------------------------------
% names = {'rACC antic' 'Neg. emotion report' 'Neg - Neu Stimulation'}
% med_results = mediation_Y_search(X, Y, M, 'names', names, 'robust')

function med_results = mediation_Y_search(X, Y, M, varargin)
    med_results = mediation_search('Y', X, Y, M, varargin{:});
end
