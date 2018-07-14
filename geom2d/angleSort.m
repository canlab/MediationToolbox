function varargout = angleSort(pts, varargin)
%ANGLESORT sort points of plane according to their angles
%
%
%   PTS2 = angleSort(PTS);
%   compute angle of points with origin, and sort points with increasing
%   angles in Counter-Clockwise direction.
%
%   PTS2 = angleSort(PTS, PTS0);
%   Compute angles between each point of PTS and PT0, which can be
%   different from origin.
%
%   PTS2 = angleSort(..., THETA0);
%   Specifies the starting angle for sorting.
%
%   [PTS2, I] = angleSort(...);
%   Also return in I the indices of PTS, such that PTS2 = PTS(I, :);
%
%
% ------
% Author: David Legland
% e-mail: david.legland@jouy.inra.fr
% Created: 2005-11-24
% Copyright 2005 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).


%   HISTORY :

% default values
pt0 = [0 0];
theta0 = 0;

if length(varargin)==1
    var = varargin{1};
    if size(var, 2)==1
        % specify angle
        theta0 = var;
    else
        pt0 = var;
    end
elseif length(varargin)==2
    pt0 = varargin{1};
    theta0 = varargin{2};
end


n = size(pts, 1);
pts2 = pts - repmat(pt0, [n 1]);
angle = lineAngle([zeros(n, 2) pts2]);
angle = mod(angle - theta0 + 2*pi, 2*pi);

[dummy, I] = sort(angle);

% format output
if nargout<2
    varargout{1} = pts(I, :);
elseif nargout==2
    varargout{1} = pts(I, :);
    varargout{2} = I;
end

        