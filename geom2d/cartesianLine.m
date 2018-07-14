function line = cartesianLine(varargin)
%CARTESIANLINE create a line with cartesian coefficients
%   l = CARTESIANLINE(a, b, c) create a line verifying equation :
%   a*x + b*x + c = 0;
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 25/05/2004.
%

if length(varargin)==1
    var = varargin{1};
    a = var(:,1);
    b = var(:,2);
    c = var(:,3);
elseif length(varargin)==3
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
end

% normalisation factor
d = a.*a + b.*b;

x0 = -a.*c./d;
y0 = -b.*c./d;
theta = atan2(-a, b);
dx = cos(theta);
dy = sin(theta);

line = [x0 y0 dx dy];