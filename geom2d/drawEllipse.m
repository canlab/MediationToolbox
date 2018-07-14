function varargout = drawEllipse(varargin)
%DRAWELLIPSE draw an ellipse on the current axis
%
%   DRAWELLIPSE(XC, YC, A, B) draw ellipse with center (XC, YC), with main
%   axis of half-length A, and second axis of half-length B.
%
%   DRAWELLIPSE(..., THETA) also specifies orientation of ellipse, given in
%   radians. Origin of orientation is (Ox) axis.
%
%   DRAWELLIPSE(PARAM) puts all parameters into one single array.
%
%   H = DRAWELLIPSE(...) also return handles to the created line objects.
%
%   -> Parameters can also be arrays. In this case, all arrays are suposed 
%   to have the same size.
%
%
%   [X Y] = DRAWELLIPSE(...) return only positions of points used to draw
%   ellipse, but does not draw the ellipse on the current axe. This allows
%   to compute intersections of ellipse, or to keep result for a later use.
%   In this case, only one ellipse path is computed. If several parameters
%   are entered, only the first one will be returned.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 11/12/2003.
%

%   HISTORY
%   08/01/2004 : add support for 2 variables (x and y arrays) returns,
%   but only for 1 ellipse
%   08/01/2004 : bug in extraction of input parameters, theta was not
%   initialized in case of array of size 1*5
%   13/08/2005 : uses radians instead of degrees


theta = 0;
if length(varargin)==1
    ellipse = varargin{1};
    x0 = ellipse(1);
    y0 = ellipse(2);
    a  = ellipse(3);
    b  = ellipse(4);
    if length(ellipse)>4
        theta = ellipse(5);
    end
elseif length(varargin)>=4
    x0 = varargin{1};
    y0 = varargin{2};
    a  = varargin{3};
    b  = varargin{4};
    if length(varargin)>4
        theta = varargin{5};
    end
else
    error('drawEllipse : please specify center x, center y and radii a and b');
end

if nargout<2
    % compute position of points to draw each ellipse
    h = zeros(length(x0), 1);
    for i=1:length(x0)
        t = 0:.01:2*pi;
        xt = x0(i) + a(i)*cos(t)*cos(theta(i)) - b(i)*sin(t)*sin(theta(i));
        yt = y0(i) + a(i)*cos(t)*sin(theta(i)) + b(i)*sin(t)*cos(theta(i));
        xt = [xt x0(i)+a(i)*cos(theta(i))];
        yt = [yt y0(i)+a(i)*sin(theta(i))];

        h(i) = line(xt, yt);
    end
    
    % return handles if needed
    if nargout>0
        varargout{1}=h;
    end
else
    % return two arrays : x and y coordinates of points

    % compute position of points used to draw first ellipse
    t = 0:.01:2*pi;
    xt = x0(1) + a(1)*cos(t)*cos(theta(1)) - b(1)*sin(t)*sin(theta(1));
    yt = y0(1) + a(1)*cos(t)*sin(theta(1)) + b(1)*sin(t)*cos(theta(1));
    xt = [xt x0(1)+a(1)*cos(theta(1))];
    yt = [yt y0(1)+a(1)*sin(theta(1))];
    
    varargout{1} = xt;
    varargout{2} = yt;
end