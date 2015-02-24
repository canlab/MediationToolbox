% Search for indirect effects X mediated by M
%
%
% See mediation_brain
%
% Example:
%
% Do robust regression search over matrix M for mediators
% ---------------------------------------------------------------------
% names = {'rACC antic' 'Neg. emotion report' 'Neg - Neu Stimulation'}
% med_results = mediation_X_search(X, Y, M, 'names', names, 'robust')

function med_results = mediation_X_search(X, Y, M, varargin)
    med_results = mediation_search('X', X, Y, M, varargin{:});
end