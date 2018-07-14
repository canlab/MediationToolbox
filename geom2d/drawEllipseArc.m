function varargout = drawEllipseArc(varargin)
%DRAWELLIPSEARC draw an ellipse on the current axis
%
%   drawEllipseArc(XC, YC, A, B) draw ellipse with center (XC, YC), with main
%   axis of half-length A, and second axis of half-length B.
%
%   drawEllipseArc(..., THETA) also specifies orientation of ellipse, given in
%   degrees. origin of orientation is (Ox) axis.
%
%   drawEllipseArc(PARAM) puts all parameters into one single array.
%
%   Parameters can also be arrays. In this case, all arrays are suposed to
%   have the same size...
%
%   TODO: precise the definition of angle for an ellipse.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 12/12/2003.
%

theta = 0;
if length(varargin)==1
    ellipse = varargin{1};
    x0 = ellipse(1);
    y0 = ellipse(2);
    a  = ellipse(3);
    b  = ellipse(4);
    if size(ellipse, 1)>6
        theta = ellipse(5)*pi/180;
        t1 = ellipse(6)*pi/180;
        t2 = ellipse(7)*pi/180;
    else
        t1 = ellipse(5)*pi/180;
        t2 = ellipse(6)*pi/180;
    end
elseif length(varargin)>=6
    x0 = varargin{1};
    y0 = varargin{2};
    a  = varargin{3};
    b  = varargin{4};
    if length(varargin)>6
        theta = varargin{5}*pi/180;
        t1 = varargin{6}*pi/180;
        t2 = varargin{7}*pi/180;
    else
        t1 = varargin{5}*pi/180;
        t2 = varargin{6}*pi/180;
    end
else
    error('drawellipse : please specify center x, center y and radii a and b');
end


h = zeros(size(x0));

for i=1:length(x0)
    t = t1:.01:t2;
    % use relation between parametric representation of ellipse, and angle :
    t = mod(atan(a(i)/b(i)*tan(t)).*(cos(t)>0) + atan2(a(i)/b(i)*tan(2*pi-t), -1).*(cos(t)<0), 2*pi);
    
    xt = x0(i) + a(i)*cos(t)*cos(theta(i)) - b(i)*sin(t)*sin(theta(i));
    yt = y0(i) + a(i)*cos(t)*sin(theta(i)) + b(i)*sin(t)*cos(theta(i));
    
    % last segment finishes exactly at t=t2 :
    t2 = t(length(t));
    xt = [xt x0(i)+a(i)*cos(t2)*cos(theta(i))-b(i)*sin(t2)*sin(theta(i))];
    yt = [yt y0(i)+a(i)*cos(t2)*sin(theta(i))+b(i)*sin(t2)*cos(theta(i))];

    h(i) = line(xt, yt);
end

if nargout>0
    varargout{1}=h;
end