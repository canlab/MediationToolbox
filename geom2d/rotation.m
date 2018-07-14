function trans = rotation(varargin)
%ROTATION return 3*3 matrix of a rotation
%
%   usage :
%   TRANS = ROTATION(THETA);
%   return the translation corresponding to angle THETA (in radians)
%   The returned matrix has the form :
%   [cos(theta) -sin(theta)  0]
%   [sin(theta)  cos(theta)  0]
%   [0           0           1]
%
%   TRANS = ROTATION(POINT, THETA);
%   TRANS = ROTATION(X0, Y0, THETA);
%   Also specify origin of rotation. The result is similar as performing
%   translation(-dx, -dy), rotation, and translation(dx, dy).
%
%
%   See also :
%   transformPoint, translation
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 06/04/2004.
%

% default values
dx = 0;
dy = 0;
theta = 0;

% get input values
if length(varargin)==1
    % only angle
    theta = varargin{1};
elseif length(varargin)==2
    % origin point (as array) and angle
    var = varargin{1};
    dx = var(1);
    dy = var(2);
    theta = varargin{2};
elseif length(varargin)==3
    % origin (x and y) and angle
    dx = varargin{1};
    dy = varargin{2};
    theta = varargin{3};
end

% compute coefs
cot = cos(theta);
sit = sin(theta);
tx =  dy*sit - dx*cot + dx;
ty = -dy*cot - dx*sit + dy;

% create transformation
trans = [cot -sit tx; sit cot ty; 0 0 1];
